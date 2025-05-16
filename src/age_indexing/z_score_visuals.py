#!/usr/bin/env python
"""
create_z_score_visuals.py

- Located in src/age_indexing
- Reads "Statistical Significance" sheet from null_age_indexes.xlsx
  (in data/age_indexes/excel_files)
- Generates circle plots:
    * min_z_scores.png
    * max_z_scores.png
- Generates bar charts:
    * min_z_scores_bar.png
    * max_z_scores_bar.png
- Outputs saved to visuals/ at project root
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Paths setup
    script_dir = os.path.dirname(__file__)
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
    excel_path = os.path.join(
        project_root, 'data', 'age_indexes', 'excel_files',
        'null_age_indexes.xlsx'
    )
    visuals_dir = os.path.join(project_root, 'visuals')
    os.makedirs(visuals_dir, exist_ok=True)

    # Read statistical significance sheet
    df = pd.read_excel(excel_path, sheet_name='Statistical Significance')
    tissues = df['Tissue'].astype(str)
    min_z = df['Min Data Z Score'].astype(float).values
    max_z = df['Max Data Z Score'].astype(float).values

    # Compute angles for circle plots
    N = len(tissues)
    angles = np.linspace(0, 2*np.pi, N, endpoint=False)

    # Circle plot: Min Z scores
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, polar=True)
    ax.scatter(angles, min_z, s=100)
    ax.set_xticks(angles)
    ax.set_xticklabels(tissues, fontsize=8)
    ax.set_title('Min Data Z Scores')
    plt.tight_layout()
    fig.savefig(os.path.join(visuals_dir, 'min_z_scores.png'))
    plt.close(fig)

    # Circle plot: Max Z scores
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, polar=True)
    ax.scatter(angles, max_z, s=100)
    ax.set_xticks(angles)
    ax.set_xticklabels(tissues, fontsize=8)
    ax.set_title('Max Data Z Scores')
    plt.tight_layout()
    fig.savefig(os.path.join(visuals_dir, 'max_z_scores.png'))
    plt.close(fig)

    # Bar chart: Min Z scores
    fig, ax = plt.subplots(figsize=(10,6))
    ax.bar(tissues, min_z)
    ax.set_xticklabels(tissues, rotation=45, ha='right')
    ax.set_xlabel('Tissue')
    ax.set_ylabel('Min Data Z Score')
    ax.set_title('Min Data Z Score by Tissue')
    plt.tight_layout()
    fig.savefig(os.path.join(visuals_dir, 'min_z_scores_bar.png'))
    plt.close(fig)

    # Bar chart: Max Z scores
    fig, ax = plt.subplots(figsize=(10,6))
    ax.bar(tissues, max_z)
    ax.set_xticklabels(tissues, rotation=45, ha='right')
    ax.set_xlabel('Tissue')
    ax.set_ylabel('Max Data Z Score')
    ax.set_title('Max Data Z Score by Tissue')
    plt.tight_layout()
    fig.savefig(os.path.join(visuals_dir, 'max_z_scores_bar.png'))
    plt.close(fig)

    print(f"Z-score visuals created in {visuals_dir}")

if __name__ == '__main__':
    main()
