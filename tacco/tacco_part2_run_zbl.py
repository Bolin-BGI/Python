import argparse
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import tacco as tc
import matplotlib.pyplot as plt


# ==========================================
# 0. 投递命令 (来自 Notebook)
# ==========================================
# # 1. M1_POD3
# !python /data/work/02.result/tacco_251118/tacco_part2_run_zbl.py \
#     --sc_name M1_POD3 \
#     --chip_id M1_D3-Y00852K1 \
#     --outdir /data/work/02.result/tacco_251118/02.result/KC/1126_modification/M1_POD3 \
#     --modifications "Bas-mig=0.05,Spi-mig=0.05"

# ==========================================
# 1. 辅助函数 (来自 Notebook)
# ==========================================

def parse_modifications(modifications_str):
    """
    Parse a string of modifications in the format 'key1=value1,key2=value2' into a dictionary.
    """
    if not modifications_str:
        return {}
    
    modifications = {}
    for item in modifications_str.split(','):
        try:
            key, value = item.split('=')
            modifications[key.strip()] = float(value)
        except ValueError:
            print(f"Warning: Format error in modification '{item}'. Expected key=value.")
    return modifications

def update_celltype_proportions2(df, modifications, ct='celltype'):
    # 计算每个粗注释的细胞比例
    celltype_proportions = df[ct].value_counts(normalize=True)

    if not modifications:
        return celltype_proportions
    
    # 过滤掉不存在于数据中的 key
    valid_mods = {k: v for k, v in modifications.items() if k in celltype_proportions.index}
    if len(valid_mods) != len(modifications):
        missing = set(modifications.keys()) - set(celltype_proportions.index)
        print(f"Warning: Following modification keys not found in data: {missing}")

    # 计算未被修改的细胞类型的当前总比例
    current_remaining_proportion = celltype_proportions.drop(valid_mods.keys()).sum()
    # 计算修改后剩余的比例
    remaining_proportion = 1 - sum(valid_mods.values())
    
    if remaining_proportion < 0:
        print("Error: Sum of modified proportions exceeds 1.0!")
        sys.exit(1)

    # 更新未被修改的细胞类型的比例
    if current_remaining_proportion > 0:
        scale_factor = remaining_proportion / current_remaining_proportion
        celltype_proportions.update(
            celltype_proportions.drop(valid_mods.keys()) * scale_factor
        )
    
    # 更新指定的粗注释细胞比例
    celltype_proportions.update(pd.Series(valid_mods))

    return celltype_proportions

def cluster_small_multiples(adata, clust_key, ncol=5, nrow=None, size=3, frameon=False, coord_x='x_adjust',coord_y='y_adjust', **kwargs):
    tmp = adata.copy()
    categories = adata.obs[clust_key].cat.categories
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(6 * ncol, 4 * nrow))
    axs = axs.flatten()

    for i, clust in enumerate(categories):
        ax = axs[i]
        tmp.obs['__highlight__'] = (adata.obs[clust_key] == clust).astype('category')
        # 注意：这里需要确保 adata.uns 中有颜色定义，如果没有可能会报错，建议加上 try-except 或默认色
        if clust_key + '_colors' in adata.uns:
             base_color = adata.uns[clust_key + '_colors'][i]
        else:
             base_color = 'red' # Fallback color
        
        tmp.uns['__highlight___colors'] = ['#d3d3d3', base_color]
        sc.pl.scatter(tmp, basis='spatial', color='__highlight__', ax=ax, title=clust, frameon=frameon, show=False, size=size, **kwargs)

    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    return fig

def markers_plot(adata, clust_key, df, **kwargs):
    adata = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    markers = {i['celltype']: [j for j in i['marker'].split(',') if j in adata.var_names] for (_, i) in df.iterrows()}
    sc.pl.dotplot(adata, markers, groupby=clust_key, **kwargs)

# ==========================================
# 2. 主逻辑
# ==========================================

def main():
    parser = argparse.ArgumentParser(description='Tacco Analysis Script')
    parser.add_argument('--sc_name', type=str, required=True, help='Single cell sample name (e.g. M1_POD3)')
    parser.add_argument('--chip_id', type=str, required=True, help='Spatial chip ID (e.g. M1_D3-Y00852K1)')
    parser.add_argument('--outdir', type=str, required=True, help='Output directory')
    parser.add_argument('--modifications', type=str, default="", help='Proportion modifications, e.g. "FB_reti_1=0.3,FB_inflam_2=0.15"') # "Key=Value,Key2=Value2"
    
    args = parser.parse_args()

    # 配置路径 (根据你的 Notebook 硬编码路径模式)
    # ---------------------------------------------------------
    # 这里的路径前缀是根据你提供的 Notebook 设定的
    # ---------------------------------------------------------
    base_data_dir = "/data/work/02.result/tacco_251118/00.data"
    
    if 'D21' in args.chip_id:
        # **【修正 D21 样本】**
        # 根据 ll 列表中的链接目标，D21 样本的实际文件在另一个目录，并且文件名不同。
        
        spatial_path = f'/data/users/zengxiaoqi/zengxiaoqi_07ef1e4aa8384c83b1fbfa53cc352045/online/002.Identify_pig_cells/07.Monkey_mapping_Cluster/{args.chip_id}/{args.chip_id}_SCT.h5ad'
        
        # **请务必确认此路径在你的环境中是正确的**
        print(f"D21 DETECTED. Using direct target path: {spatial_path}")
    
    else:
        # **【修正其他样本】**
        spatial_path = os.path.join(base_data_dir, f"sp_h5ad/{args.chip_id}_monkey_cells_SCT.h5ad")
    
    ref_path = os.path.join(base_data_dir, f"sc_ref_h5ad_add_M-WT/{args.sc_name}_anno_1111_with_MWT.h5ad")
    meta_sc_path = os.path.join(base_data_dir, "all_KC_meta_1027.csv")
    meta_sp_path = os.path.join(base_data_dir, f"mian_pred/pred_csv_gz/{args.chip_id}/pred.csv.gz")
    markers_path = os.path.join(base_data_dir, "scRNA_KC_markers_for_tacco_1118.csv")
    
    annotation_key = 'kc_1020'
    celltype_sp_key = 'pred_celltype'
    selected_cells = ['KC']
    min_cell_count = 20

    # 准备输出目录
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    os.chdir(args.outdir) # 切换工作目录，方便保存图片
    
    print(f"Start processing: {args.sc_name} | {args.chip_id}")
    print(f"Modifications: {args.modifications}")
    
    ## --- 【新增/修改】 文件名安全处理逻辑 ---
    modifications_str = args.modifications
    # 1. 判断并清理特殊字符
    if modifications_str:
        # 如果字符串非空，则进行清理
        cleaned_modifications = modifications_str.replace("=", "_").replace(",", "_")
    else:
        # 如果是空字符串，则直接将清理后的结果设为空字符串
        cleaned_modifications = "" 

    # 2. 构建后缀：只有 cleaned_modifications 非空时，safe_mod_suffix 才包含下划线
    safe_mod_suffix = f"_{cleaned_modifications}" if cleaned_modifications else ""
    


    # 1. Load Reference
    print("Loading Reference...")
    reference = ad.read_h5ad(ref_path)
    meta_sc = pd.read_csv(meta_sc_path, index_col=0)
    
    # Intersect cells
    cells = reference.obs_names.intersection(meta_sc.index)
    reference = reference[cells, :]
    reference.obs[annotation_key] = meta_sc.loc[reference.obs.index, annotation_key]
    
    # 2. Load Spatial
    print("Loading Spatial...")
    spatial_raw = ad.read(spatial_path)
    meta_sp = pd.read_csv(meta_sp_path, index_col=0)
    spatial_raw.obs[celltype_sp_key] = meta_sp.loc[spatial_raw.obs.index, celltype_sp_key]
    spatial_raw.X = spatial_raw.layers['counts'] # 使用 counts 层

    # 3. Coordinate logic (D21 check)
    if 'D21' in spatial_path: # Check if D21 is in the path/filename
        if 'spatial' in spatial_raw.obsm:
             spatial_raw.obsm['X_spatial'] = spatial_raw.obsm['spatial']
    else:
        if 'x_adjust' in spatial_raw.obs and 'y_adjust' in spatial_raw.obs:
            new_spatial_coords = spatial_raw.obs[['x_adjust', 'y_adjust']].values
            spatial_raw.obsm['spatial'] = new_spatial_coords
            spatial_raw.obsm['X_spatial'] = spatial_raw.obsm['spatial']
        else:
            print("Warning: x_adjust/y_adjust not found, using existing spatial if available.")

    # 4. Filter Spatial
    spatial = spatial_raw[spatial_raw.obs[celltype_sp_key].isin(selected_cells), :].copy()

    # 5. Filter Reference (min cells)
    table = reference.obs[annotation_key].value_counts()
    removed = table[table < min_cell_count].index.tolist()
    reference = reference[~reference.obs[annotation_key].isin(removed), :]
    print(f"Removed cell types < {min_cell_count} cells: {removed}")

    # 6. Calculate Proportions & Modifications
    raw_prob = reference.obs[annotation_key].value_counts() / reference.obs[annotation_key].size
    raw_prob.to_csv('raw_prob.csv')
    print("Raw SC Proportions:\n", raw_prob.head())
    
    #### ----- Modifications -----
    mod_dict = parse_modifications(args.modifications)
    celltype_proportions = update_celltype_proportions2(reference.obs, mod_dict, ct=annotation_key)
    celltype_proportions.to_csv(f'modified_sc_prop{safe_mod_suffix}.csv')
    print("Modified SC Proportions:\n", celltype_proportions.head())

    # 7. Plot Reference Markers
    markers_df = pd.read_csv(markers_path)
    # 设置绘图分辨率
    plt.rcParams['figure.dpi'] = 200
    plt.rcParams['savefig.dpi'] = 200
    
    try:
        markers_plot(reference, annotation_key, markers_df, show=False, save="reference.png")
    except Exception as e:
        print(f"Error plotting reference markers: {e}")

    # 8. Run Tacco Annotate
    print("Running Tacco Annotation...")
    prob = celltype_proportions.copy()
    tc.tl.annotate(spatial, reference, annotation_key,
                   result_key='pred_celltype', annotation_prior=prob, verbose=False)

    # 9. Process Predictions
    prediction = spatial.obsm['pred_celltype'].copy()
    prediction['pred_celltype'] = prediction.idxmax(1)
    prediction['max_score'] = prediction.iloc[:, :-1].max(1)
    
    pred_prob = prediction['pred_celltype'].value_counts(normalize=True)
    pred_prob.to_csv(f'predict_sp_prob{safe_mod_suffix}.csv')
    
    # MSE Calculation
    # Align indices for calculation
    common_idx = pred_prob.index.intersection(prob.index)
    mse = ((pred_prob[common_idx] - prob[common_idx]) ** 2).mean()
    print(f"MSE: {mse:.5f}")

    # 10. Save Results
    prediction.to_csv(f'pred{safe_mod_suffix}.csv.gz')
    
    # Update spatial objects with prediction
    spatial.obs['pred_celltype'] = prediction['pred_celltype']
    
    # Initialize raw with original types, then update with predictions
    spatial_raw.obs['pred_celltype'] = spatial_raw.obs[celltype_sp_key].astype(object)
    spatial_raw.obs.loc[spatial.obs_names, 'pred_celltype'] = prediction['pred_celltype']
    
    # 11. Final Plots
    print("Generating plots...")
    celltype_col = 'pred_celltype' # Use the name used in obs
    
    # Plot 1: Spatial scatter (subset)
    sc.pl.scatter(spatial, basis='spatial', color=celltype_col, show=False)
    plt.savefig(f"pred_{annotation_key}{safe_mod_suffix}.png")
    plt.close()

    # Plot 2: Split (subset)
    try:
        fig = cluster_small_multiples(spatial, clust_key=celltype_col, frameon=False)
        fig.savefig(f"spatial_pred_{annotation_key}{safe_mod_suffix}_split.png", dpi=200)
        plt.close(fig)
    except Exception as e:
        print(f"Error plotting split subset: {e}")

    # Plot 3: Markers spatial
    try:
        markers_plot(spatial, 'pred_celltype', markers_df, show=False, save="spatial.png")
    except Exception as e:
         print(f"Error plotting spatial markers: {e}")

    # Plot 4: Spatial scatter (all)
    sc.pl.scatter(spatial_raw, basis='spatial', color=celltype_col, show=False)
    plt.savefig(f"pred_{annotation_key}{safe_mod_suffix}_all.png")
    plt.close()

    # Plot 5: Split (all)
    try:
        fig = cluster_small_multiples(spatial_raw, clust_key=celltype_col, frameon=False)
        fig.savefig(f"spatial_pred_{annotation_key}{safe_mod_suffix}_split_all.png", dpi=200)
        plt.close(fig)
    except Exception as e:
        print(f"Error plotting split all: {e}")

    print("Analysis Done!")

if __name__ == "__main__":
    main()