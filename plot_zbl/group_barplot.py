import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

def group_barplot(adata, group_x, group_y, prefix="plot"):
    """
    绘制不同 group_x 中 group_y 类别分布的柱状图（分面显示所有 group_y 类别）
    
    参数:
    adata: AnnData 对象
    group_x: 分组变量1（如 'monkey'）
    group_y: 分组变量2（如 'bulk_label'）
    prefix: 图片保存前缀
    使用示例 ：plot_group_comparison(adata, group_x='monkey', group_y='bulk_label', prefix="monkey_vs_celltype")
    """
    # 生成交叉表
    cross_table = pd.crosstab(adata.obs[group_x], adata.obs[group_y])
    
    # 重塑数据为长格式（便于 Seaborn 绘图）
    data_summary = cross_table.reset_index().melt(id_vars=group_x, var_name=group_y, value_name='count')
    
    # 设置颜色和样式
    palette = sns.color_palette("Paired", len(data_summary[group_x].unique()))
    sns.set(style="whitegrid", font_scale=1.2)
    
    # 创建分面网格
    g = sns.FacetGrid(
        data_summary, 
        col=group_y,  # 按 group_y 分面
        col_wrap=4,   # 每行4个子图
        height=4,     # 子图高度
        aspect=1.2,   # 子图宽高比
        sharey=False, # 每个子图独立y轴
        sharex=False  # 每个子图独立x轴
    )
    
    # 在每个子图中绘制柱状图
    g.map_dataframe(
        sns.barplot, 
        x=group_x, 
        y='count', 
        palette=palette,
        order=data_summary[group_x].unique()  # 保持x轴顺序
    )
    
    # 在每个柱子上方添加数值标签
    for ax in g.axes.flat:
        for p in ax.patches:
            ax.annotate(
                f"{int(p.get_height())}",  # 显示的数值
                (p.get_x() + p.get_width() / 2., p.get_height()),  # 位置
                ha='center', va='center',  # 对齐方式
                fontsize=10, color='black',  # 字体大小和颜色
                xytext=(0, 5),  # 偏移量
                textcoords='offset points'
            )
    
    # 调整子图标签和标题
    g.set_titles(col_template="{col_name}", size=14)  # 子图标题（group_y类别名）
    g.set_axis_labels(x_var=group_x, y_var="Cell Count")  # 设置x轴和y轴标签
    
    # 调整整体布局
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle(f'{group_y} Distribution Across {group_x}', y=0.98, size=16)
    
    # 保存图片
    g.savefig(f"{prefix}_{group_x}_vs_{group_y}.png", dpi=300, bbox_inches='tight')
    plt.show()