import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def relative_day_percent_by_monkey(
    adata, group, day_col, monkey_col, prefix="",
    monkey_select=None, use_log2=False
):
    """
    绘制adata对象中 cluster 的比例随 day 变化的图片，按 monkey 分类。
    每只猴子每种细胞类型的 D0 数据作为基准值（设为1），其余天数显示相对于 D0 的变化比例。
    使用每个猴子在每个时间点的细胞类型占比（而非 count）进行比较，避免文库大小影响。

    参数:
     adata: AnnData对象
     group: 表示cluster的列名 (str)
     day_col: 表示day的列名 (str)
     monkey_col: 表示monkey的列名 (str)
     prefix: 保存图片时的文件名前缀 (str)
     monkey_select: list, 只展示指定的monkey类别，例如["M1","M4"]，默认 None 表示全部展示
     use_log2: bool, 是否使用 log2(fold change) 绘图，默认 False（使用 relative proportion）

     使用示例：
     # 使用 relative proportion
     cluster_day_percent_by_monkey(adata, 'celltype_zbl', 'day', 'monkey', prefix="Test", monkey_select=["M1", "M4"])

     # 使用 log2 fold change
     cluster_day_percent_by_monkey(adata, 'celltype_zbl', 'day', 'monkey', prefix="Test", use_log2=True)
    """

    # 验证输入参数
    if group not in adata.obs.columns:
        raise ValueError(f"'{group}' 不在 adata.obs 中。")
    if day_col not in adata.obs.columns:
        raise ValueError(f"'{day_col}' 不在 adata.obs 中。")
    if monkey_col not in adata.obs.columns:
        raise ValueError(f"'{monkey_col}' 不在 adata.obs 中。")

    # 提取需要的列并确保为字符串类型
    data = adata.obs[[group, day_col, monkey_col]].copy()
    data = data.applymap(str)

    # 将 day 列转换为分类数据类型，并按时间顺序排列
    data[day_col] = data[day_col].astype('category')
    data[day_col] = data[day_col].cat.set_categories(['D0', 'D1', 'D3', 'D7', 'D14', 'D21'], ordered=True)

    # 计算每个 day 和 monkey 的总计数（用于后续计算比例）
    data_summary = data.groupby([day_col,monkey_col,group]).size().reset_index(name='count')
    day_monkey_total = data_summary.groupby([day_col, monkey_col])['count'].transform('sum')
    data_summary['day_monkey_total'] = day_monkey_total  # 合并入 data_summary
    
    # 在比例计算后做筛选 
    if monkey_select is not None:
        data_summary = data_summary[data_summary[monkey_col].isin(monkey_select)]
        if data_summary.empty:
            raise ValueError("筛选后的数据为空，请检查monkey_select的取值。")

    # 计算每个细胞类型在该 day + monkey 中的占比（避免文库大小影响）
    data_summary['proportion'] = data_summary['count'] / data_summary['day_monkey_total'] 
    
    # 提取每个 monkey + group 组合的 D0 占比作为基准
    baseline = (data_summary[data_summary[day_col] == 'D0'].set_index([monkey_col, group])['proportion'])
    print(baseline)
    
    # 映射 D0 的 proportion 回原数据
    data_summary['baseline_prop'] = (
        data_summary.set_index([monkey_col, group]).index.map(baseline.get)
    )

    # 检查是否有缺失值（即某些组合没有 D0）
    if data_summary['baseline_prop'].isna().any():
        raise ValueError("某些 monkey + group 组合中没有 D0 数据，请检查数据完整性。")

    # 计算相对于 D0 的比例变化
    data_summary['relative_proportion'] = data_summary['proportion'] / data_summary['baseline_prop']

    # 强制设置 D0 为 1.0（防止浮点误差）
    data_summary.loc[data_summary[day_col] == 'D0', 'relative_proportion'] = 1.0

    # 新增：计算 log2 fold change
    data_summary['log2_fold_change'] = np.log2(data_summary['relative_proportion'])

    # 添加 label 字段用于 facet 分面
    data_summary['label'] = data_summary[group].astype(str)

    data_summary.to_csv(f"{prefix}_{group}_by_monkey_relative.csv", index=True) 
    print(data_summary)
    
    # 创建颜色映射，使用更明显的调色板
    monkeys = data_summary[monkey_col].unique()
    palette = sns.color_palette("bright", len(monkeys))  # 使用 'bright' 调色板
    color_map = {monkey: palette[i] for i, monkey in enumerate(monkeys)}

    # 设置图像大小和风格
    plt.figure(figsize=(20, 15))
    sns.set(style="whitegrid")

    # 选择绘图的 y 值和 Y 轴标签
    if use_log2:
        y_col = 'log2_fold_change'
        y_label = 'log2(Fold Change)'
        title_suffix = " (log2 Fold Change)"
    else:
        y_col = 'relative_proportion'
        y_label = 'Relative Proportion to D0'
        title_suffix = " (Relative Proportion)"

    # 创建分面网格，按 celltype 分面
    g = sns.FacetGrid(
        data_summary, col="label", col_wrap=4, height=4, aspect=1,
        sharex=False, sharey=False  # y轴不共享，便于观察个体变化趋势
    )

    # 绘制折线图
    g.map_dataframe(
        sns.lineplot, x=day_col, y=y_col, hue=monkey_col, marker="o", palette=color_map
    )

    # 调整每个子图样式
    for ax in g.axes.flat:
        ax.set_xticks(data_summary[day_col].cat.categories)
        ax.set_xticklabels(data_summary[day_col].cat.categories, rotation=0)
        ax.set_xlabel('Day', fontsize=12)
        ax.set_ylabel(y_label, fontsize=12)
        ax.tick_params(axis='both', labelsize=10)
        ax.legend(title='Monkey', fontsize=10)
        if use_log2:
            ax.axhline(0.0, color='gray', linestyle='--', linewidth=1)  # 标注 log2(1) = 0
        else:
            ax.axhline(1.0, color='gray', linestyle='--', linewidth=1)  # 标注 D0 基准线

    # 调整整体布局并添加总标题
    plt.subplots_adjust(top=0.92, hspace=0.3, wspace=0.2)
    g.fig.suptitle(f'Clusters Relative Changes Over Days by Monkey{title_suffix}', fontsize=16)

    # 保存图片
    output_file = f"{prefix}_{group}_by_monkey_relative.png"
    g.savefig(output_file, dpi=300, bbox_inches='tight')

    # 打印图形
    plt.show()