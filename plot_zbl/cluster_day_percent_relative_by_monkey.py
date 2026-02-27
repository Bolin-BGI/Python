import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

def cluster_day_percent_relative_by_monkey(
    adata, group, day_col, monkey_col, prefix="",
    monkey_select=None    # 新增参数
):
    """
    绘制adata对象中 cluster 的比例随 day 变化的图片，按 monkey 分类。
    以 D0 时间点的比例为基线（设为1），展示后续时间点相对于 D0 的变化倍数。
    可指定只展示部分 monkey 的数据。
    
    参数:
     adata: AnnData对象
     group: 表示cluster的列名 (str)
     day_col: 表示day的列名 (str)
     monkey_col: 表示monkey的列名 (str)
     prefix: 保存图片时的文件名前缀 (str)
     monkey_select: list, 只展示指定的monkey类别，例如["M1","M4"]，默认 None 表示全部展示
     
     使用示例：
     cluster_day_percent_relative_by_monkey(adata, 'celltype_zbl', 'day', 'monkey', prefix="Test", monkey_select=["M1", "M4"])
    
    """
    # 验证输入参数
    if group not in adata.obs.columns:
        raise ValueError(f"'{group}' 不在 adata.obs 中。")
    if day_col not in adata.obs.columns:
        raise ValueError(f"'{day_col}' 不在 adata.obs 中。")
    if monkey_col not in adata.obs.columns:
        raise ValueError(f"'{monkey_col}' 不在 adata.obs 中。")
    
    # 提取需要的列并确保为字符串类型 (使用 astype(str) 替代废弃的 applymap)
    data = adata.obs[[group, day_col, monkey_col]].copy()
    data = data.astype(str)

    # 将 day 列转换为分类数据类型，并按时间顺序排列
    data[day_col] = data[day_col].astype('category')
    data[day_col] = data[day_col].cat.set_categories(['D0', 'D1', 'D3', 'D7', 'D14', 'D21'], ordered=True)

    # 计算每个 day、每个 group、每个 monkey 的计数
    data_summary = data.groupby([day_col, group, monkey_col], observed=False).size().reset_index(name='count')

    # 计算每个 day 和 monkey 的总计数
    day_group_total = data_summary.groupby([day_col, monkey_col], observed=False)['count'].transform('sum')

    # 计算绝对比例
    data_summary['proportion'] = data_summary['count'] / day_group_total * 100

    # 在 data_summary 中添加 label 信息
    data_summary['label'] = data_summary[group].astype(str)

    # ================= 修改核心逻辑区 =================
    # 1. 提取 D0 时刻的比例作为每个 monkey 和每个 group 的基线 (Baseline)
    d0_baseline = data_summary[data_summary[day_col] == 'D0'][[group, monkey_col, 'proportion']].copy()
    d0_baseline = d0_baseline.rename(columns={'proportion': 'd0_proportion'})
    
    # 2. 将 D0 的基线比例合并回主数据框
    data_summary = data_summary.merge(d0_baseline, on=[group, monkey_col], how='left')
    
    # 3. 计算相对比例 (当前比例 / D0比例)。如果 D0 没细胞(即为0或NaN)，结果会自动成为 NaN，合理规避除0错误
    # 这里 D0 会变成 1.0
    data_summary['relative_proportion'] = data_summary['proportion'] / data_summary['d0_proportion']
    
    # 将缺失值（比如某些猴子在D0完全没有某个细胞群，导致无法计算相对倍数）替换为 NaN 并在作图时跳过
    data_summary['relative_proportion'] = data_summary['relative_proportion'].replace([np.inf, -np.inf], np.nan)
    # ===================================================

    # ====== 新增：在比例计算后做筛选 ======
    if monkey_select is not None:
        data_summary = data_summary[data_summary[monkey_col].isin(monkey_select)]
        if data_summary.empty:
            raise ValueError("筛选后的数据为空，请检查monkey_select的取值。")
    # ====================================
    
    # 调试输出
    print("Data summary head (Relative):")
    print(data_summary.head())
    
    # 创建颜色映射，使用更明显的调色板
    monkeys = data[monkey_col].unique()
    palette = sns.color_palette("bright", len(monkeys))  # 使用 'bright' 调色板
    color_map = {monkey: palette[i] for i, monkey in enumerate(monkeys)}

    # 绘制分面折线图
    plt.figure(figsize=(20, 15))
    sns.set_theme(style="whitegrid") # 更新为 set_theme

    g = sns.FacetGrid(
        data_summary, col="label", col_wrap=4, height=4, aspect=1, 
        palette=palette, sharex=False, sharey=False # 因为各个细胞变化倍数差异可能很大，取消sharey通常效果更好
    )
    
    # 注意这里将 y 替换成了 'relative_proportion'
    g.map_dataframe(
        sns.lineplot, x=day_col, y="relative_proportion", hue=monkey_col, marker="o", palette=color_map
    )

    # 调整图像布局并在每个子图中添加图例和参考线
    for ax in g.axes.flat:
        ax.set_xticks(data_summary[day_col].cat.categories)
        ax.set_xticklabels(data_summary[day_col].cat.categories, rotation=0)
        ax.set_xlabel('Day', fontsize=12)
        ax.set_ylabel('Relative Proportion (D0 = 1)', fontsize=12) # 坐标轴改名
        ax.tick_params(axis='both', labelsize=10)
        
        # 添加一条 y=1 的参考线，视觉上很容易看出相比于 D0 是上调还是下调
        ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        
        ax.legend(title='monkey', fontsize=10)

    # 调整整体布局并添加总标题
    plt.subplots_adjust(top=0.92, hspace=0.3, wspace=0.3) # 增加了wspace防止不共享Y轴时标签重叠
    g.fig.suptitle('Relative Clusters Proportion Over Days by Monkey (Baseline D0=1)', fontsize=16)

    # 保存图片
    g.savefig(f"{prefix}_{group}_by_monkey_relative_percent.png", dpi=300)

    # 打印图形
    plt.show()
    
    
    