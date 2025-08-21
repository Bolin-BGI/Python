import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

# 0429修改

def cluster_day_percent_by_monkey(
    adata, group, day_col, monkey_col, prefix="",
    monkey_select=None    # 新增参数
):
    """
    绘制adata对象中 cluster 的比例随 day 变化的图片，按 monkey 分类
    可指定只展示部分 monkey 的数据。
    
    参数:
     adata: AnnData对象
     group: 表示cluster的列名 (str)
     day_col: 表示day的列名 (str)
     monkey_col: 表示monkey的列名 (str)
     prefix: 保存图片时的文件名前缀 (str)
     monkey_select: list, 只展示指定的monkey类别，例如["M1","M4"]，默认 None 表示全部展示
     
     使用示例：
     cluster_day_percent_by_monkey(adata, 'celltype_zbl', 'day', 'monkey', prefix="Test", monkey_select=["M1", "M4"])
    
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

    # 计算每个 day、每个 group、每个 monkey 的计数
    data_summary = data.groupby([day_col, group, monkey_col]).size().reset_index(name='count')

    # 计算每个 day 和 group 的总计数
    day_group_total = data_summary.groupby([day_col, monkey_col])['count'].transform('sum')

    # 计算比例
    data_summary['proportion'] = data_summary['count'] / day_group_total * 100

    # 在 data_summary 中添加 label 信息
    data_summary['label'] = data_summary[group].astype(str)

    # ====== 新增：在比例计算后做筛选 ======
    if monkey_select is not None:
        data_summary = data_summary[data_summary[monkey_col].isin(monkey_select)]
        if data_summary.empty:
            raise ValueError("筛选后的数据为空，请检查monkey_select的取值。")
    # ====================================
    
    # 调试输出
    print("Data summary head:")
    print(data_summary.head())
    
    # 确保数据中没有缺失值
    print("Checking for missing values:")
    print(data_summary.isnull().sum())

    # 创建颜色映射，使用更明显的调色板
    monkeys = data[monkey_col].unique()
    palette = sns.color_palette("bright", len(monkeys))  # 使用 'bright' 调色板
    color_map = {monkey: palette[i] for i, monkey in enumerate(monkeys)}

    # 绘制分面折线图
    plt.figure(figsize=(20, 15))
    sns.set(style="whitegrid")

    g = sns.FacetGrid(
        data_summary, col="label", col_wrap=4, height=4, aspect=1, 
        palette=palette, sharex=False, sharey=True
    )
    g.map_dataframe(
        sns.lineplot, x=day_col, y="proportion", hue=monkey_col, marker="o", palette=color_map
    )

    # 调整图像布局并在每个子图中添加图例
    for ax in g.axes.flat:
        ax.set_xticks(data_summary[day_col].cat.categories)
        ax.set_xticklabels(data_summary[day_col].cat.categories, rotation=0)
        ax.set_xlabel('Day', fontsize=12)
        ax.set_ylabel('Proportion (%)', fontsize=12)
        ax.tick_params(axis='both', labelsize=10)
        ax.legend(title='monkey', fontsize=10)

    # 调整整体布局并添加总标题
    plt.subplots_adjust(top=0.92,hspace=0.3, wspace=0.2)
    g.fig.suptitle('Clusters Proportion Over Days by Monkey', fontsize=16)

    # 保存图片
    g.savefig(f"{prefix}_{group}_by_monkey_percent.png", dpi=300)

    # 打印图形
    plt.show()

















# def cluster_day_percent_by_monkey(adata, group, day_col, monkey_col, prefix=""):
#     """
#     绘制adata对象中 cluster 的比例随 day 变化的图片，按 monkey 分类

#     参数:
#     adata: AnnData对象
#     group: 表示cluster的列名 (str)
#     day_col: 表示day的列名 (str)
#     monkey_col: 表示monkey的列名 (str)
#     prefix: 保存图片时的文件名前缀 (str)
#     """
#     # 验证输入参数
#     if group not in adata.obs.columns:
#         raise ValueError(f"'{group}' 不在 adata.obs 中。")
#     if day_col not in adata.obs.columns:
#         raise ValueError(f"'{day_col}' 不在 adata.obs 中。")
#     if monkey_col not in adata.obs.columns:
#         raise ValueError(f"'{monkey_col}' 不在 adata.obs 中。")
    
#     # 提取需要的列并确保为字符串类型
#     data = adata.obs[[group, day_col, monkey_col]].copy()
#     data = data.applymap(str)

#     # 将 day 列转换为分类数据类型，并按时间顺序排列
#     data[day_col] = data[day_col].astype('category')
#     data[day_col] = data[day_col].cat.set_categories(['D0', 'D1', 'D3', 'D7', 'D14', 'D21'], ordered=True)

#     # 计算每个 day、每个 group、每个 monkey 的计数
#     data_summary = data.groupby([day_col, group, monkey_col]).size().reset_index(name='count')

#     # 计算每个 day 和 group 的总计数
#     day_group_total = data_summary.groupby([day_col, group])['count'].transform('sum')

#     # 计算比例
#     data_summary['proportion'] = data_summary['count'] / day_group_total * 100

#     # 在 data_summary 中添加 label 信息
#     data_summary['label'] = data_summary[group].astype(str)

#     # 调试输出
#     print("Data summary head:")
#     print(data_summary.head())
    
#     # 确保数据中没有缺失值
#     print("Checking for missing values:")
#     print(data_summary.isnull().sum())

#     # 创建颜色映射，使用更明显的调色板
#     monkeys = data[monkey_col].unique()
#     palette = sns.color_palette("bright", len(monkeys))  # 使用 'bright' 调色板
#     color_map = {monkey: palette[i] for i, monkey in enumerate(monkeys)}

#     # 绘制分面折线图
#     plt.figure(figsize=(20, 15))
#     sns.set(style="whitegrid")

#     g = sns.FacetGrid(data_summary, col="label", col_wrap=4, height=4, aspect=1, palette=palette,sharex=False,sharey=True)
#     g.map_dataframe(sns.lineplot, x=day_col, y="proportion", hue=monkey_col, marker="o", palette=color_map)

#     # 调整图像布局并在每个子图中添加图例
#     for ax in g.axes.flat:
#         ax.set_xticks(data_summary[day_col].cat.categories)
#         ax.set_xticklabels(data_summary[day_col].cat.categories, rotation=0)
#         ax.set_xlabel('Day', fontsize=12)
#         ax.set_ylabel('Proportion (%)', fontsize=12)
#         ax.tick_params(axis='both', labelsize=10)
#         ax.legend(title='monkey', fontsize=10)

#     # 调整整体布局并添加总标题
#     plt.subplots_adjust(top=0.92)
#     g.fig.suptitle('Clusters Proportion Over Days by Monkey', fontsize=16)

#     # 保存图片
#     g.savefig(f"percent_{prefix}_{group}_by_monkey.png", dpi=300)

#     # 打印图形
#     plt.show()

# 使用示例
# cluster_day_percent_by_monkey(adata, 'celltype_zbl', 'day', 'monkey', prefix="Test")