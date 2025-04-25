import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def load_best_df(path: Path) -> pd.DataFrame:
    """
    Load the 'best' sheet from the given Excel file.
    """
    if not path.exists():
        raise FileNotFoundError(f"Excel file not found: {path}")
    df = pd.read_excel(path, sheet_name='best')
    # Ensure correct column names
    df = df.rename(columns={df.columns[1]: 'Metric'})
    return df[['Tissue', 'Metric']]


def prepare_output_dir(root: Path) -> Path:
    """
    Ensure the visuals directory exists at the project root.
    """
    visuals_dir = root / 'visuals'
    visuals_dir.mkdir(parents=True, exist_ok=True)
    return visuals_dir


def plot_bar(df: pd.DataFrame, output_path: Path, title: str) -> None:
    """
    Plot a bar chart of Metric by Tissue, sorted by descending absolute Metric.
    """
    # Sort by absolute Metric descending
    df = df.copy()
    df['abs_metric'] = df['Metric'].abs()
    df = df.sort_values(by='abs_metric', ascending=False)

    # Plot
    fig, ax = plt.subplots()
    ax.bar(df['Tissue'], df['Metric'])
    ax.set_xlabel('Tissue')
    ax.set_ylabel('Age Index')
    ax.set_title(title)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def main():
    # Determine project root (two levels up from this script)
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]

    # Paths to Excel files
    excel_dir = root / 'data' / 'age_indexes' / 'excel_files'
    max_excel = excel_dir / 'a1_max_age_indexes.xlsx'
    min_excel = excel_dir / 'a1_min_age_indexes.xlsx'

    # Load dataframes
    df_max = load_best_df(max_excel)
    df_min = load_best_df(min_excel)

    # Prepare visuals directory
    visuals_dir = prepare_output_dir(root)

    # Plot and save
    plot_bar(df_max, visuals_dir / 'a1_max_bargraph.png', 'Max Age Index by Tissue (Best Run)')
    plot_bar(df_min, visuals_dir / 'a1_min_bargraph.png', 'Min Age Index by Tissue (Best Run)')

    print(f"Bar graphs saved to {visuals_dir}")


if __name__ == '__main__':
    main()
