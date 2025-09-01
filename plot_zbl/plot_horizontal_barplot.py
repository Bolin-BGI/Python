import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_horizontal_barplot(
    adata,
    obs_col="final_annotation",
    figsize=(8, 12),
    fontsize=12,
    save=False,
    save_path="celltype_barplot.png"
):
    """
    绘制细胞比例横向柱状图（由多到少排序），在每个柱子后标注百分比值

    eg:
    # 调整字体和尺寸，并保存
    plot_horizontal_barplot(
        adata_ref,
        obs_col="final_annotation",
        figsize=(10, 14),
        fontsize=14,
        save=False,
        save_path="final_annotation_barplot.png"
    )

    参数
    ----
    adata : AnnData
        输入的 AnnData 对象
    obs_col : str, default="final_annotation"
        选择 obs 中的哪一列进行统计
    figsize : tuple, default=(8, 12)
        图像尺寸 (width, height)
    fontsize : int, default=12
        字体大小
    save : bool, default=False
        是否保存图片
    save_path : str, default="celltype_barplot.png"
        保存路径（仅当 save=True 生效）
    """

    # 统计比例
    cell_counts = adata.obs[obs_col].value_counts(normalize=True) * 100
    cell_counts = cell_counts.sort_values(ascending=False)

    # 转为 DataFrame
    df = cell_counts.reset_index()
    df.columns = ['Celltype', 'Percentage']

    # 绘图
    plt.figure(figsize=figsize)
    ax = sns.barplot(
        y='Celltype',
        x='Percentage',
        data=df,
        order=df['Celltype'],
        color='steelblue'
    )

    plt.xlabel('Percentage (%)', fontsize=fontsize)
    plt.ylabel(obs_col, fontsize=fontsize+1)
    plt.title(f'Cell type composition (sorted by percentage) - {obs_col}', fontsize=fontsize+2)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize-3)

    # 在柱子后添加百分比标签
    for i, (perc, celltype) in enumerate(zip(df['Percentage'], df['Celltype'])):
        ax.text(
            perc + 0.1,  # 在柱子末端稍微偏右
            i,
            f"{perc:.2f}%",
            va='center',
            ha='left',
            fontsize=fontsize-3
        )

    plt.tight_layout()

    # 保存 or 显示
    if save:
        plt.savefig(save_path, dpi=300)
        plt.show()
        plt.close()
        
    else:
        plt.show()
