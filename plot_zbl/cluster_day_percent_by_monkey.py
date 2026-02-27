import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
from scipy.interpolate import PchipInterpolator

def cluster_day_percent_by_monkey(
    adata, group, day_col, monkey_col, prefix="",
    monkey_select=None,
    smooth=False  # ====== 新增：是否平滑曲线参数 ======
):
    """
    绘制adata对象中 cluster 的绝对比例随 day 变化的图片，按 monkey 分类。
    可指定只展示部分 monkey 的数据，支持曲线平滑。
    
    参数:
     adata: AnnData对象
     group: 表示cluster的列名 (str)
     day_col: 表示day的列名 (str)
     monkey_col: 表示monkey的列名 (str)
     prefix: 保存图片时的文件名前缀 (str)
     monkey_select: list, 只展示指定的monkey类别，例如["M1","M4"]
     smooth: bool, 是否对折线进行平滑处理 (默认 False)
    """
    # 1. 验证输入参数
    if group not in adata.obs.columns:
        raise ValueError(f"'{group}' 不在 adata.obs 中。")
    if day_col not in adata.obs.columns:
        raise ValueError(f"'{day_col}' 不在 adata.obs 中。")
    if monkey_col not in adata.obs.columns:
        raise ValueError(f"'{monkey_col}' 不在 adata.obs 中。")
    
    # 2. 提取数据并转为字符串
    data = adata.obs[[group, day_col, monkey_col]].copy()
    data = data.astype(str)

    # 3. 设置 Day 的分类顺序
    days_order = ['D0', 'D1', 'D3', 'D7', 'D14', 'D21']
    data[day_col] = data[day_col].astype('category')
    data[day_col] = data[day_col].cat.set_categories(days_order, ordered=True)

    # 4. 计算绝对比例 (分子 / 分母)
    data_summary = data.groupby([day_col, group, monkey_col], observed=False).size().reset_index(name='count')
    day_group_total = data_summary.groupby([day_col, monkey_col], observed=False)['count'].transform('sum')
    data_summary['proportion'] = data_summary['count'] / day_group_total * 100
    data_summary['label'] = data_summary[group].astype(str)

    # 5. 筛选 Monkey
    if monkey_select is not None:
        data_summary = data_summary[data_summary[monkey_col].isin(monkey_select)]
        if data_summary.empty:
            raise ValueError("筛选后的数据为空，请检查monkey_select的取值。")

    # 调试输出
    print("Data summary head:")
    print(data_summary.head())  
    
    # 6. 准备绘图颜色
    monkeys = data_summary[monkey_col].unique()
    # 你的附图用的是典型的红蓝配色，这里可以用 'Set1' 贴近原图风格，也可使用 'bright' 调色板
    palette = sns.color_palette("Set1", len(monkeys))   
    color_map = {monkey: palette[i] for i, monkey in enumerate(monkeys)}

    # ====== 核心：自定义平滑绘图函数 ======
    def _plot_lines(data, x, y, hue, color_dict, is_smooth, **kwargs):
        # 建立类别到数字的映射，方便插值
        categories = data[x].cat.categories
        cat_to_num = {cat: i for i, cat in enumerate(categories)}
        
        for hue_name, group_data in data.groupby(hue, observed=False):
            # 确保按时间顺序排列
            group_data = group_data.sort_values(by=x).dropna(subset=[y])
            if group_data.empty:
                continue
                
            x_vals = group_data[x].map(cat_to_num).values
            y_vals = group_data[y].values
            current_color = color_dict.get(hue_name, 'black')

            if is_smooth and len(x_vals) >= 3:
                # 开启平滑且数据点>=3时进行插值
                x_smooth = np.linspace(x_vals.min(), x_vals.max(), 300)
                # 使用 Pchip 防止曲线产生负数比例或剧烈抖动
                interp = PchipInterpolator(x_vals, y_vals)
                y_smooth = interp(x_smooth)
                y_smooth = np.clip(y_smooth, 0, None) # 彻底截断负值
                
                # 画平滑线 + 原始散点
                plt.plot(x_smooth, y_smooth, color=current_color, linewidth=2.5)
                plt.scatter(x_vals, y_vals, color=current_color, s=20)
            else:
                # 不平滑或数据点不够时，画普通折线
                plt.plot(x_vals, y_vals, marker='o', color=current_color, linewidth=2.5, markersize=4)
    # ====================================

    # 7. 开始绘图
    plt.figure(figsize=(20, 15))
    sns.set_theme(style="whitegrid")

    # 注意这里 sharey=False 才能让每个细胞群自适应 Y 轴（和你的附图一致）
    g = sns.FacetGrid(
        data_summary, col="label", col_wrap=4, height=3.5, aspect=1.2, 
        sharex=False, sharey=False
    )
    
    # 映射自定义函数
    g.map_dataframe(
        _plot_lines, x=day_col, y="proportion", hue=monkey_col, 
        color_dict=color_map, is_smooth=smooth
    )

    # 8. 调整刻度与外观
    for ax in g.axes.flat:
        # 重置 X 轴刻度，因为插值使用了数字索引
        ax.set_xticks(range(len(days_order)))
        ax.set_xticklabels(days_order, rotation=45) 
        
        ax.set_xlabel('Day', fontsize=12)
        ax.set_ylabel('Proportion (%)', fontsize=12)
        ax.tick_params(axis='both', labelsize=10)
        
        # 添加外边框（类似附图的方框风格）
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.2)

    # 9. 手动添加统一的图例
    handles = [plt.Line2D([0], [0], color=color_map[m], lw=2.5, label=m) for m in monkeys]
    g.fig.legend(handles=handles, title='Monkey', loc='center right', bbox_to_anchor=(0.98, 0.5))

    # 调整布局留出图例空间
    plt.subplots_adjust(top=0.90, right=0.92, hspace=0.4, wspace=0.3)
    title_suffix = "(Smoothed)" if smooth else ""
    g.fig.suptitle(f'Absolute Clusters Proportion Over Days {title_suffix}', fontsize=16, fontweight='bold')

    # 保存与展示
    g.savefig(f"{prefix}_{group}_percent_smooth_{smooth}.png", dpi=300, bbox_inches='tight')
    plt.show()
