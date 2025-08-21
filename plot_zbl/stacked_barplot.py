# stacked_barplot 绘制不同时间点细胞类型比例变化的堆叠柱状图

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

def get_color_palette(color_palette, num_colors):
    """
    生成颜色方案。
    如果颜色数量超过默认方案，使用 matplotlib 的颜色映射生成更多颜色。
    """
    if color_palette in ["tab20", "tab28"]:
        # 使用 Scanpy 默认颜色方案
        if color_palette == "tab20":
            colors = sc.pl.palettes.default_20
        elif color_palette == "tab28":
            colors = sc.pl.palettes.default_28
        # 如果颜色数量不足，循环使用颜色方案
        return (colors * (num_colors // len(colors) + 1))[:num_colors]
    else:
        # 使用 matplotlib 的颜色映射生成更多颜色
        cmap = plt.get_cmap(color_palette)
        return cmap(np.linspace(0, 1, num_colors))


def stacked_barplot(
    adata,
    group_x='day',
    group_y='celltype_zbl',
    prefix='output',
    figsize=(10, 6),
    color_palette='tab20'
):
    """
    绘制不同时间点细胞类型比例变化的堆叠柱状图，并在图中添加比例标签。

    参数：
    adata: Scanpy对象
    group_x: 时间列名，默认'day'
    group_y: 细胞类型列名，默认'celltype_zbl'
    prefix: 输出图片名前缀，默认'output'
    figsize: 图片尺寸，默认(10,6)
    color_palette: 颜色方案，默认'tab20'; matplotlib 提供了许多内置的颜色映射（colormaps），你可以通过 plt.colormaps() 查看所有可用的颜色映射。
    """
    # 提取数据并计算比例
    df = adata.obs[[group_x, group_y]].copy()
    prop_table = pd.crosstab(df[group_x], df[group_y], normalize='index')

    # 生成颜色映射
    categories = prop_table.columns
    colors = get_color_palette(color_palette, len(categories))

    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)

    # 绘制堆叠柱状图
    bars = prop_table.plot.bar(
        stacked=True,
        ax=ax,
        color=colors,
        width=0.8,
        edgecolor='white'
    )

    # 添加比例标签
    for i, (index, row) in enumerate(prop_table.iterrows()):
        cumulative_height = 0  # 记录当前堆叠高度
        for j, value in enumerate(row):
            if value > 0.03:  # 只标注非零比例
                # 计算标签的位置
                label_position = cumulative_height + value / 2
                # 添加文本标签，保留小数点后一位
                ax.text(
                    i,  # x 坐标
                    label_position,  # y 坐标
                    f"{value:.1%}",  # 标签内容，格式为百分比
                    ha='center',  # 水平居中对齐
                    va='center',  # 垂直居中对齐
                    fontsize=8,  # 字体大小
                    color='black'  # 字体颜色
                )
            cumulative_height += value  # 更新累计高度

    # 设置坐标轴和标签
    ax.set_xlabel(group_x, fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,  # X 轴标签旋转 45 度
        ha='right',
        rotation_mode='anchor'
    )

    # 设置图例
    ax.legend(
        title=group_y,
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=False
    )

    # 设置边框和网格
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(False)

    # 保存图片
    plt.tight_layout()
    plt.savefig(f"{prefix}_stacked_barplot.png", dpi=300, bbox_inches='tight')
    plt.close()