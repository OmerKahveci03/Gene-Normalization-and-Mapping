#!/usr/bin/env python
"""
- Reads "Statistical Significance" sheet from null_age_indexes.xlsx
- Generates:
    • visuals/min_z_scores.png
    • visuals/max_z_scores.png
    • visuals/min_z_scores_bar.png
    • visuals/max_z_scores_bar.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def create_circle_plot(values, labels, output_path, sort_by_magnitude=True):
    orig = np.array(values, dtype=float)
    if sort_by_magnitude:
        idx = np.argsort(np.abs(orig))
        r = np.abs(orig[idx])
    else:
        idx = np.argsort(orig)
        r = orig[idx]

    lbls = labels[idx]
    N = len(r)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)

    # determine ring positions
    max_r = r.max()
    ring_max = int(np.ceil(max_r / 2)) * 2
    rings = np.arange(2, ring_max + 1, 2)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, polar=True)
    ax.scatter(angles, r, s=50, zorder=2)

    # radial rings
    ax.set_yticks(rings)
    labels_rings = [f"{int(rt)}" for rt in rings]
    labels_rings[-1] = ''  # omit outermost ring label
    ax.set_yticklabels(labels_rings, fontsize=8)
    ax.set_ylim(0, ring_max)
    ax.set_rlabel_position(0)
    ax.set_xticks([])

    # tissue names on outer ring
    label_radius = ring_max * 1.05
    for angle, lbl in zip(angles, lbls):
        deg = np.rad2deg(angle)
        rot = deg
        if 90 < deg < 270:
            rot += 180
        if rot > 180:
            rot -= 360
        ha = 'left' if (deg <= 90 or deg >= 270) else 'right'

        ax.text(
            angle, label_radius, lbl,
            rotation=rot,
            rotation_mode='anchor',
            ha=ha,
            va='center',
            fontsize=7,
            zorder=3
        )

    plt.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def create_bar_graph(values, labels, title, ylabel, output_path):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(labels, values, zorder=2)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_xlabel('Tissue')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def main():
    # set up paths
    here = os.path.dirname(__file__)
    project = os.path.abspath(os.path.join(here, '..', '..'))
    excel_path = os.path.join(
        project, 'data', 'age_indexes', 'excel_files',
        'null_age_indexes.xlsx'
    )
    out = os.path.join(project, 'visuals')
    os.makedirs(out, exist_ok=True)

    # load
    df = pd.read_excel(excel_path, sheet_name='Statistical Significance')
    tissues = df['Tissue'].astype(str).values
    min_z = df['Min Data Z Score'].astype(float).values
    max_z = df['Max Data Z Score'].astype(float).values

    # circle plots
    create_circle_plot(
        min_z, tissues,
        os.path.join(out, 'min_z_scores.png'),
        sort_by_magnitude=True
    )
    create_circle_plot(
        max_z, tissues,
        os.path.join(out, 'max_z_scores.png'),
        sort_by_magnitude=False
    )

    # bar graphs
    create_bar_graph(
        min_z, tissues,
        "Min Data Z Score by Tissue", "Min Data Z Score",
        os.path.join(out, 'min_z_scores_bar.png')
    )
    create_bar_graph(
        max_z, tissues,
        "Max Data Z Score by Tissue", "Max Data Z Score",
        os.path.join(out, 'max_z_scores_bar.png')
    )

    print(f"Z-score visuals created in {out}")


if __name__ == '__main__':
    main()
