import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
from scipy.interpolate import PchipInterpolator

def cluster_day_percent_by_monkey(
    adata, group, day_col, monkey_col, prefix="",
    monkey_select=None, 
    color=None,
    smooth=False,
    figsize=(20, 15)  # ====== 新增：图片大小参数，默认 (20, 15) ======
):
    """
    绘制adata对象中 cluster 的绝对比例随 day 变化的图片，按 monkey 分类。
    可指定只展示部分 monkey 的数据，支持曲线平滑与自定义颜色、图片大小。
    
    参数:
     adata: AnnData对象
     group: 表示cluster的列名 (str)
     day_col: 表示day的列名 (str)
     monkey_col: 表示monkey的列名 (str)
     prefix: 保存图片时的文件名前缀 (str)
     monkey_select: list, 只展示指定的monkey类别，例如["M1","M4"]
     color: list, 自定义monkey的颜色，顺序需与 monkey_select 一致，例如["red", "blue"]
     smooth: bool, 是否对折线进行平滑处理 (默认 False)
     figsize: tuple, 输出图片的宽和高 (默认 (20, 15))
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
    
    if color is not None:
        if monkey_select is not None:
            if len(color) != len(monkey_select):
                raise ValueError("提供的 color 列表长度必须与 monkey_select 列表长度一致！")
            color_map = {m: c for m, c in zip(monkey_select, color)}
        else:
            if len(color) != len(monkeys):
                raise ValueError("未提供 monkey_select 时，color 列表长度必须与实际存在的 monkey 数量一致！")
            color_map = {m: c for m, c in zip(monkeys, color)}
    else:
        palette = sns.color_palette("Set1", len(monkeys)) 
        color_map = {monkey: palette[i] for i, monkey in enumerate(monkeys)}

    # 自定义平滑绘图函数 
    def _plot_lines(data, x, y, hue, color_dict, is_smooth, **kwargs):
        categories = data[x].cat.categories
        cat_to_num = {cat: i for i, cat in enumerate(categories)}
        
        for hue_name, group_data in data.groupby(hue, observed=False):
            group_data = group_data.sort_values(by=x).dropna(subset=[y])
            if group_data.empty:
                continue
                
            x_vals = group_data[x].map(cat_to_num).values
            y_vals = group_data[y].values
            current_color = color_dict.get(hue_name, 'black')

            if is_smooth and len(x_vals) >= 3:
                x_smooth = np.linspace(x_vals.min(), x_vals.max(), 300)
                interp = PchipInterpolator(x_vals, y_vals)
                y_smooth = interp(x_smooth)
                y_smooth = np.clip(y_smooth, 0, None) 
                
                plt.plot(x_smooth, y_smooth, color=current_color, linewidth=2.5)
                plt.scatter(x_vals, y_vals, color=current_color, s=20)
            else:
                plt.plot(x_vals, y_vals, marker='o', color=current_color, linewidth=2.5, markersize=4)

    # 7. 开始绘图
    # ====== 修改：移除了无效的 plt.figure(figsize=(20, 15)) ======
    sns.set_theme(style="whitegrid")

    g = sns.FacetGrid(
        data_summary, col="label", col_wrap=4, height=3.5, aspect=1.2, 
        sharex=False, sharey=False
    )
    
    # ====== 新增：在这里强制应用你传入的 figsize ======
    g.fig.set_size_inches(figsize[0], figsize[1])
    # ===============================================

    g.map_dataframe(
        _plot_lines, x=day_col, y="proportion", hue=monkey_col, 
        color_dict=color_map, is_smooth=smooth
    )

    # 8. 调整刻度与外观
    for ax in g.axes.flat:
        ax.set_xticks(range(len(days_order)))
        ax.set_xticklabels(days_order, rotation=0) 
        
        ax.set_xlabel('Day', fontsize=12)
        ax.set_ylabel('Proportion (%)', fontsize=12)
        ax.tick_params(axis='both', labelsize=10)
        
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.2)

    # 9. 手动添加统一的图例
    legend_monkeys = monkey_select if monkey_select is not None else monkeys
    handles = [plt.Line2D([0], [0], color=color_map[m], lw=2.5, label=m) for m in legend_monkeys]
    g.fig.legend(handles=handles, title='Monkey', loc='center right', bbox_to_anchor=(0.98, 0.5))

    # 调整布局留出图例空间
    plt.subplots_adjust(top=0.90, right=0.92, hspace=0.4, wspace=0.3)
    title_suffix = "(Smoothed)" if smooth else ""
    g.fig.suptitle(f'Absolute Clusters Proportion Over Days {title_suffix}', fontsize=16, fontweight='bold')

    # 10. 保存与展示
    g.savefig(f"{prefix}_{group}_percent_smooth_{smooth}.png", dpi=300, bbox_inches='tight')
    g.savefig(f"{prefix}_{group}_percent_smooth_{smooth}.pdf", bbox_inches='tight')
    
    plt.show()


