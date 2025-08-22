# Python

## plot\_zbl

### stacked\_barplot

```python
plot_zbl.stacked_barplot(
    adata,
    group_x='group',
    group_y='celltype_51',
    prefix='Group',
    figsize=(15, 6),
    color_palette='tab20',
)
```

<img width="2499" height="992" alt="image" src="https://github.com/user-attachments/assets/73daad08-6f62-433e-b10e-b99281fdc861" />

---

### cluster\_split

```python
cluster_split(adata, 'celltype_1229', 'split')
```

<img width="944" height="636" alt="image" src="https://github.com/user-attachments/assets/afa79649-4b76-42cf-bf31-011c3e4597da" />

---

### cluster\_group\_percent

```python
cluster_group_percent(adata, 'celltype_51', 'monkey', prefix='line')
```

<img width="1719" height="860" alt="image" src="https://github.com/user-attachments/assets/0ae961df-885b-4bd8-9037-ba71f2184ffc" />

---

### cluster\_day\_percent\_by\_monkey

```python
cluster_day_percent_by_monkey(
    adata,
    'celltype_51',
    'day',
    'monkey',
    prefix="M14",
    monkey_select=["M1", "M4"]
)
```

<img width="1740" height="852" alt="image" src="https://github.com/user-attachments/assets/0d9c5da9-2b8c-41b9-a89c-f32847afd979" />

---

### split\_group\_stacked\_barplot

```python
split_group_stacked_barplot(
    adata,
    group_sub='monkey',
    group_x='day',
    group_y='celltype_51',
    save_path='split_monkey_celltype_51_stacked_barplot.png'
)
```

<img width="2541" height="1207" alt="image" src="https://github.com/user-attachments/assets/6271f707-85f0-41d6-abe3-c6127a5518fa" />

---

### group\_barplot

```python
group_barplot(adata, group_x='monkey', group_y='celltype_dmt', prefix="barplot")
```

<img width="1566" height="682" alt="image" src="https://github.com/user-attachments/assets/62860f7b-b4c4-4533-8572-3ee7ce0224ea" />

---
