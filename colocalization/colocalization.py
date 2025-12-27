import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- 命令行参数解析部分 ---
parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input', help='输入的 h5ad 文件路径')
parser.add_argument('-M', '--meta', help='包含细胞注释(如 celltype)的 csv 文件')
parser.add_argument('-F', '--save_filename', help='结果保存的 csv 文件名')
parser.add_argument('-K', '--k', type=int, default=20, help='计算空间邻居时的邻居数量 (K-nearest neighbors)')

args = parser.parse_args()
print(args)


def count_adjacent_categories(A, categories):
    """
    计算不同细胞类型之间的空间邻接次数。
    
    参数:
    A: 空间邻接矩阵 (Sparse Matrix, N_cells x N_cells)，表示细胞间的连接关系。
    categories: 细胞类型注释 (Series 或 Array)，长度为 N_cells。
    
    返回:
    category_counts_df: DataFrame，行是源细胞类型，列是邻居细胞类型，值是邻接边的数量。
    """
    # 1. 获取所有唯一的细胞类型，并建立索引映射
    unique_categories = np.unique(categories)
    
    # 2. 将由字符串组成的类别列表转换为 One-hot (独热) 编码的稀疏矩阵
    # 结果矩阵 category_matrix 的形状为 (N_cells, N_categories)
    # 如果第 i 个细胞属于第 j 类，则 (i, j) 为 1，否则为 0
    n_nodes = len(categories)
    n_categories = len(unique_categories)
    rows = np.arange(n_nodes)
    # 找到每个细胞对应的类别索引
    cols = [np.where(unique_categories == category)[0][0] for category in categories]
    data = np.ones(n_nodes)
    category_matrix = csr_matrix((data, (rows, cols)), shape=(n_nodes, n_categories))
    
    # 3. 核心矩阵运算：计算邻接节点的类别分布
    # A (N x N) 点乘 category_matrix (N x C) = adjacent_category_matrix (N x C)
    # 结果矩阵的第 i 行第 j 列表示：细胞 i 有多少个邻居属于类别 j
    adjacent_category_matrix = A.dot(category_matrix)
    
    # 4. 汇总统计：将以"细胞"为单位的数据，聚合为以"细胞类型"为单位
    category_counts_df = pd.DataFrame(0, index=unique_categories, columns=unique_categories)
    
    for category in unique_categories:
        # 找到属于当前 category 的所有细胞的行索引
        category_indices = np.where(categories == category)[0]
        rows_in_category = rows[category_indices] # 这里的 rows 其实就是 range(n_nodes)，可以直接用 indices
        
        # 提取这些细胞的邻居情况，并按列求和
        # 逻辑：对于所有属于类型 A 的细胞，它们共有多少个类型 B 的邻居？
        # sum(axis=0) 将该类别所有细胞的邻居向量相加
        category_counts_df.loc[category] = adjacent_category_matrix[rows_in_category, :].sum(axis=0)
    
    return category_counts_df


def hypergeom_test(category_counts_df):
    """
    使用超几何检验评估细胞类型邻接的显著性。
    
    参数:
    category_counts_df: 观察到的邻接计数矩阵 (由上一个函数生成)。
    
    返回:
    results_df: 包含 P值、Log2FC 等统计结果的 DataFrame。
    """
    # M: 总体大小 (Total Population)
    # 整个空间图中所有边的总数（邻接总次数）
    M = category_counts_df.to_numpy().sum()
    
    results_df = pd.DataFrame(columns=['From','To','Counts','Log2FC','Pvalue', 'Percentage'])
    
    # 遍历每一对 (From, To) 细胞类型
    for i, row in category_counts_df.iterrows():
        # N: 抽样数量 (Sample Size)
        # 也就是类型 'i' (From) 发出的所有边的总数（该类型细胞的所有邻居总数）
        N = row.sum()
        
        for j in category_counts_df.columns:
            # n: 总体中成功的数量 (Successes in Population)
            # 也就是整个图中，指向类型 'j' (To) 的边的总数
            n = category_counts_df[j].sum()
            
            # x: 样本中成功的数量 (Observed Successes)
            # 实际观察到的：类型 'i' 的邻居中，属于类型 'j' 的数量
            x = row[j]
            
            # 计算 Log2 Fold Change (观察频率 / 期望频率)
            # 期望频率 = (类型j的总边数 / 总边数M)
            # 观察频率 = (类型i到j的边数x / 类型i的总边数N)
            foldchange = np.log2((x / N) / (n / M))
            
            # 计算 P-value (单尾检验：观察值及更极端情况的概率)
            # hypergeom.cdf(x-1) 计算的是 P(X < x)，即累积分布函数
            # 1 - cdf(x-1) 计算的是 P(X >= x)，即富集的显著性
            p_value = 1 - hypergeom(M, n, N).cdf(x - 1)
            
            # Percentage: 类型 i 的邻居中，类型 j 占比多少
            results_df = results_df._append({
                'From': i, 
                'To': j, 
                'Counts': x, 
                'Log2FC': foldchange, 
                'Pvalue': p_value, 
                'Percentage': x/N
            }, ignore_index=True)
    
    return results_df

# --- 主程序逻辑 ---

# 1. 读取数据
print(f"Loading h5ad: {args.input}")
adata = sc.read_h5ad(args.input)

print(f"Loading meta: {args.meta}")
meta = pd.read_csv(args.meta, index_col=0)

# 2. 取交集：确保只保留 meta 中存在的细胞
# 注意：先取交集，防止 meta 里有 h5ad 里没有的细胞导致报错
common_cells = adata.obs_names.intersection(meta.index)
print(f"Cells in h5ad: {adata.shape[0]}, Cells in meta: {meta.shape[0]}, Overlap: {len(common_cells)}")

if len(common_cells) == 0:
    raise ValueError("Error: No common cells found between h5ad and meta file!")

adata = adata[common_cells, :]
meta = meta.loc[common_cells, :]

# 3. 覆盖 metadata
adata.obs = meta

# 4. 强制更新 spatial 坐标
# 既然你的 meta 里有 x_adjust 和 y_adjust，我们以此为准
if 'x_adjust' in adata.obs.columns and 'y_adjust' in adata.obs.columns:
    print("Updating adata.obsm['spatial'] using 'x_adjust' and 'y_adjust' from metadata...")
    # 确保转换为 float 类型，防止字符串导致的报错
    spatial_coords = adata.obs[['x_adjust', 'y_adjust']].values.astype(float)
    adata.obsm['spatial'] = spatial_coords
else:
    print("Warning: 'x_adjust' or 'y_adjust' not found in metadata. Using existing adata.obsm['spatial'].")

# 5. 检查并移除坐标为 NaN 的细胞
# 检查每一行（每个细胞）是否有 NaN
if 'spatial' in adata.obsm:
    nan_mask = np.isnan(adata.obsm['spatial']).any(axis=1)
    if nan_mask.sum() > 0:
        print(f"⚠️ Warning: Found {nan_mask.sum()} cells with NaN coordinates. Removing them...")
        adata = adata[~nan_mask, :]
    else:
        print("Check passed: No NaN values in spatial coordinates.")
else:
    raise ValueError("Error: adata.obsm['spatial'] is missing!")

# ===============================================

print(f"Final cell count for analysis: {adata.shape[0]}")




# 3. 计算空间邻居图
# 这会在 adata.obsp['spatial_connectivities'] 生成邻接矩阵 A
sq.gr.spatial_neighbors(adata, n_neighs=args.k)

# 4. 运行自定义函数：统计邻接频次
category_counts = count_adjacent_categories(adata.obsp['spatial_connectivities'], adata.obs['pred_celltype'])

# 5. 运行自定义函数：计算超几何检验统计量
hypergeom_results = hypergeom_test(category_counts)

# 6. 保存结果
hypergeom_results.to_csv(args.save_filename)