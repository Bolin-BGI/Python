def cluster_split(adata, group: str, prefix: str, basis: str = 'X_umap',) -> None:
    """
    根据分组变量绘制降维图，并将所有分组的图像保存到一个 PNG 文件中。
    
    参数：
    adata: AnnData 对象
    group: 分组变量（adata.obs 中的列名）
    prefix: 文件名前缀
    basis: 降维数据的表示，默认为 'X_umap'
    
    示例调用: cluster_split(adata, 'leiden_res_0.50', 'SMC', basis='X_umap')
    """
    
    import scanpy as sc
    import matplotlib.pyplot as plt
    import math
    import warnings
    
    # 参数验证
    if group not in adata.obs.columns:
        raise ValueError(f"'{group}' 不在 adata.obs 中。可用的选项有: {adata.obs.columns}")
    if basis not in adata.obsm_keys():
        raise ValueError(f"'{basis}' 不在 adata.obsm 中。可用的选项有: {adata.obsm_keys()}")
    
    # 抑制 FutureWarning 和 UserWarning
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
      
    # 将 group 列转换为分类数据类型，获取所有唯一的分组标签，并排序
    adata.obs[group] = adata.obs[group].astype('category') 
    unique_groups = adata.obs[group].cat.categories
    # 计算行列数以适应所有子图
    n_groups = len(unique_groups)  # 计算行列数以适应所有子图
    ncols = 4  # 每行的子图数量
    nrows = math.ceil(n_groups / ncols)  # 计算需要的行数

    # 创建一个大的图形
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 7, nrows * 5))

    # 对每个分组绘制 UMAP 图
    for idx, i in enumerate(unique_groups):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row, col] if nrows > 1 else axes[col]
        
        # 绘制 UMAP 图，并高亮显示当前分组的细胞
        sc.pl.embedding(adata, basis=basis,color=group, groups=[i], ax=ax, size=3, show=False)
        ax.set_title(f"{group}: {i}")

    # 删除多余的子图
    for idx in range(n_groups, nrows * ncols):
        row = idx // ncols
        col = idx % ncols
        fig.delaxes(axes[row, col] if nrows > 1 else axes[col])

    # 调整布局并保存图形为 PNG 文件
    plt.tight_layout()
    output_path = f"cluster_split_{prefix}_{group}.png"
    fig.savefig(output_path, dpi=100)
    plt.show(fig)
    plt.close(fig)

    print(f"PNG 文件已保存为 {output_path}")


