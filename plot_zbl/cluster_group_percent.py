import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc


def cluster_group_percent(adata, group_y, group_x, prefix="Test"):
    """
    绘制adata对象中，按 group_y 分组后，各 group_x 在每个 group_y 中的比例

    参数:
    adata: AnnData对象
    group_y: 分组依据，如细胞类型，字符串
    group_x: 统计比例的子分组，如天数，字符串
    prefix: 图片保存前缀
    """
    # 提取需要的列
    data = adata.obs[[group_y, group_x]].copy()

    # 将 group_x 和 group_y 都强制转换为 category，保持原有顺序
    data[group_x] = data[group_x].astype('category')
    data[group_y] = data[group_y].astype('category')

    # 统计各 group_y 中 group_x 的数量
    summary = data.groupby([group_y, group_x]).size().reset_index(name='count')

    # 计算每个 group_y 下 group_x 数的比例
    summary['proportion'] = summary.groupby(group_y)['count'].transform(lambda x: x / x.sum() * 100)

    # 为绘图方便
    summary['group_x_label'] = summary[group_x].astype(str)
    summary['group_y_label'] = summary[group_y].astype(str)

    # 调色板
    palette = sns.color_palette("Paired", len(summary[group_x].unique()))

    # facet 各 group_y
    plt.figure(figsize=(15, 10))
    sns.set(style="whitegrid")

    # g = sns.FacetGrid(summary, col="group_y_label", col_wrap=4, height=4, aspect=1, palette=palette, sharex=False, sharey=True)
    # g.map_dataframe(sns.lineplot, x=group_x, y="proportion", hue=group_x, marker="o", palette=palette)
    g = sns.FacetGrid(summary, col="group_y_label", col_wrap=4, height=4, aspect=1, sharex=False, sharey=True)
    g.set_titles("{col_name}")  # 移除 "group_y_label=" 前缀，仅显示值
    g.map_dataframe(sns.lineplot, x=group_x, y="proportion", marker="o", color='b')
    
    g.add_legend(title=group_x, bbox_to_anchor=(1.05, 1) ) # loc='upper left'

    # 细调
    for ax, cname in zip(g.axes.flat, summary['group_y_label'].unique()):
        # 当前 group_y 的所有 group_x（有可能不是全部 group_x，故需筛选）
        xticks = summary.loc[summary['group_y_label']==cname, group_x].unique()
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, rotation=0)
        ax.set_xlabel(group_x, fontsize=12)
        ax.set_ylabel("Proportion (%)", fontsize=12)

    plt.subplots_adjust(top=0.9, hspace=0.5)
    g.fig.suptitle(f'Proportion of {group_x} within {group_y}', fontsize=16)
    g.savefig(f"{prefix}_{group_y}_by_{group_x}_percent.png", dpi=300)

    # 绘制总的折线图
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=summary, x=group_x, y="proportion", hue=group_y, marker="o", palette="tab10")
    plt.title(f'Proportion of {group_x} in each {group_y}', fontsize=16)
    plt.xlabel(group_x, fontsize=14)
    plt.ylabel('Proportion (%)', fontsize=14)
    plt.xticks(summary[group_x].cat.categories)
    plt.legend(title=group_y, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    #plt.savefig(f"{prefix}_{group_y}_by_{group_x}_percent_all.png", dpi=300)

    plt.show()

# 用法示例
# group_x_percent_in_group_y(adata, group_y="celltype", group_x="day")





# old_version
# def cluster_group_percent(adata, group_y, group_x, prefix="Test"):
#     """
#     绘制adata对象中 group_y 的比例随 day 变化的图片

#     参数:
#     adata: AnnData对象
#     group_y: 表示 group_y 的列名 (如 celltype 或其它分组信息,str)
#     group_x: 表示分组变量（如 day 或其他分组信息）的列名 (str)
#     prefix: 保存图片时的文件名前缀 (str)
#     """
#     # 提取需要的列
#     data = adata.obs[[group_y, group_x]].copy()

#     # 将 group_x 列转换为分类数据类型，并指定分类的顺序
#     data[group_x] = data[group_x].astype('category')

#     # 计算每个 group_x 中 每个 group_y 的计数
#     data_summary = data.groupby([group_x, group_y]).size().reset_index(name='count')

#     # 计算每个 group_y 在每个 group_x 的比例
#     data_summary['proportion'] = data_summary.groupby(group_x)['count'].transform(lambda x: x / x.sum() * 100)

#     # 在 data_summary 中添加 group_y_label 信息
#     data_summary['group_y_label'] = data_summary[group_y].astype(str)

#     # 使用 seaborn 默认调色板
#     palette = sns.color_palette("Paired", len(data_summary[group_y].unique()))

#     # 绘制折线图
#     plt.figure(figsize=(15, 10))
#     sns.set(style="whitegrid")

#     g = sns.FacetGrid(data_summary, col="group_y_label", col_wrap=4, height=4, aspect=1, palette=palette,sharex=False,sharey=True)
#     g.map_dataframe(sns.lineplot, x=group_x, y="proportion", hue=group_y, marker="o", palette=palette)
#     g.add_legend()

#     # 调整图像布局
#     for ax in g.axes.flat:
#         ax.set_xticks(data_summary[group_x].cat.categories)
#         ax.set_xticklabels(data_summary[group_x].cat.categories, rotation=0)
#         ax.set_xlabel('', fontsize=16)
#         ax.set_ylabel('Proportion (%)', fontsize=14)
#         ax.tick_params(axis='both', labelsize=12)

#     plt.subplots_adjust(top=0.9)
#     g.fig.suptitle('group_ys Proportion Over Days', fontsize=16)
#     g.savefig(f"{prefix}_{group_y}_percent.png", dpi=300)

#     # 绘制总的折线图
#     plt.figure(figsize=(10, 6))
#     sns.lineplot(data=data_summary, x=group_x, y="proportion", hue=group_y, marker="o", palette=palette)
#     plt.title('Total group_ys Proportion Over Days', fontsize=16)
#     plt.xlabel(group_x, fontsize=14)
#     plt.ylabel('Proportion (%)', fontsize=14)
#     plt.xticks(data_summary[group_x].cat.categories, rotation=0)
#     plt.legend(title='group_y', bbox_to_anchor=(1.05, 1), loc='upper left')
#     plt.tight_layout()
#     plt.savefig(f"{prefix}_{group_y}_percent_all.png", dpi=300)

#     # 打印图形
#     plt.show()
