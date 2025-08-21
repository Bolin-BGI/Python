import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def split_group_stacked_barplot(
    adata,
    group_sub: str,          # 分组变量（如 'monkey'）
    group_x: str,            # 时间变量（如 'day'，需是分类列）
    group_y: str,            # 细胞类型列（如 'bulk_label'）
    color: str = 'tab20',    # 高对比度颜色映射名称（如 'tab20', 'Set3', 'viridis'）
    figsize: tuple = (25, 12),
    save_path: str = None
):
    """
    绘制分组变量（group_sub）下不同时间（group_x）的细胞类型比例堆积柱状图
    
    参数:
    adata: AnnData 对象
    group_sub: 子图分组列名（如猴子 'monkey'，决定画几个子图）
    group_x: 行名：时间列名（如 'day'，需是分类变量）
    group_y: 列名：细胞类型（如 'bulk_label'）
    color: plt.cm.get_cmap 的颜色映射名称（默认 'tab20'）, 高对比度颜色映射名称（如 'tab10','tab20','tab20b','tab20c', 'Set1','Set2','Set3', 'viridis'）
    figsize: 图像大小（默认 (25, 12)）
    save_path: 图片保存路径名称（如 "monkey_day_celltype_stacked.png"） 如果 为 None 则只显示图片

    """

    
    # --- 数据准备 ---
    # 检查分类顺序
    if not isinstance(adata.obs[group_x].dtype, pd.CategoricalDtype):
        raise ValueError(f"'{group_x}' 必须是分类变量（Categorical）！")
    group_x_order = adata.obs[group_x].cat.categories.tolist()
    
    # group_sub 顺序
    if isinstance(adata.obs[group_sub].dtype, pd.CategoricalDtype):
        group_sub_values = adata.obs[group_sub].cat.categories.tolist()
    else:
        group_sub_values = sorted(adata.obs[group_sub].unique())
    n_subplots = len(group_sub_values)
    
    # 生成交叉表并计算比例
    cross_table = pd.crosstab(
        index=[adata.obs[group_sub], adata.obs[group_x]],
        columns=adata.obs[group_y]
    )
    prop_table = cross_table.div(cross_table.sum(axis=1), axis=0).fillna(0) * 100
    prop_table = prop_table.reset_index()
    
    # 确保 group_x 顺序正确
    prop_table[group_x] = pd.Categorical(
        prop_table[group_x], 
        categories=group_x_order,
        ordered=True
    )
    
    # 转换为长格式
    prop_long = prop_table.melt(
        id_vars=[group_sub, group_x],
        var_name='cell_type',
        value_name='percentage'
    )
    
    # --- 颜色设置 ---
    cell_types = prop_long['cell_type'].unique()
    n_cell_types = len(cell_types)
    cmap = plt.cm.get_cmap(color)
    if n_cell_types > cmap.N:
        colors = [cmap(i % cmap.N) for i in range(n_cell_types)]
    else:
        colors = cmap(np.linspace(0, 1, n_cell_types))
    color_dict = {ct: colors[i] for i, ct in enumerate(cell_types)}
    
    # --- 子图布局 ---
    max_cols = 3  # 每行最多4个
    ncols = min(n_subplots, max_cols)
    nrows = math.ceil(n_subplots / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()
    
    # --- 绘图 ---
    for i, group_val in enumerate(group_sub_values):
        ax = axes[i]
        df_group = prop_long[prop_long[group_sub] == group_val].pivot(
            index=group_x, columns='cell_type', values='percentage'
        ).reindex(group_x_order)
        
        # 绘制堆积柱状图
        bars = df_group.plot.bar(
            stacked=True,
            ax=ax,
            color=[color_dict[ct] for ct in df_group.columns],
            width=0.8,
            edgecolor='white',
            legend=False  # 先不画，后面单独加
        )
        
        # 添加百分比标签
        for bar in bars.containers:
            ax.bar_label(
                bar,
                labels=[f'{val:.1f}%' if val > 5 else '' for val in bar.datavalues],
                label_type='center',
                fontsize=11,
                color='black'
            )
        
        # 设置子图标题和标签
        ax.set_title(f'{group_sub} {group_val}', fontsize=12)
        ax.set_xlabel(group_x, fontsize=10)
        ax.set_xticklabels(df_group.index, rotation=0)
        if i % ncols == 0:
            ax.set_ylabel('Percentage (%)', fontsize=12)
        else:
            ax.set_ylabel('')
        ax.grid(False)
        ax.set_ylim(0, 100)
        
        # 每个子图右上角加图例
        ax.legend(
            df_group.columns,
            title='Cell Type',
            loc='upper left',
            bbox_to_anchor=(1.01, 1),
            fontsize=9,
            title_fontsize=10,
            frameon=False
        )
    
    # 多余的子图隐藏
    for j in range(i+1, nrows*ncols):
        fig.delaxes(axes[j])
    
    plt.tight_layout(rect=[0, 0, 1, 1])
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

