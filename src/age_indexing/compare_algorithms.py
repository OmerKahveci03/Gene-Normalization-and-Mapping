#!/usr/bin/env python
import os
import pandas as pd

def load_and_prepare(file_path, sheet_name=None):
    """
    Load an Excel file without a header and assign
    column names: 'Tissue' and 'AgeIndex'. Optionally load a specific sheet.
    """
    df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)
    # Assuming the first two columns represent the tissue name and age index.
    df = df.iloc[:, :2]
    df.columns = ['Tissue', 'AgeIndex']
    return df

def compare_mode(df_a0, df_a1, mode):
    """
    Merge algorithm 0 and algorithm 1 dataframes by tissue,
    then compute the difference (algorithm 1 - algorithm 0).
    """
    merged = pd.merge(df_a0, df_a1, on='Tissue', suffixes=('_0', '_1'))
    merged['Difference'] = merged['AgeIndex_1'] - merged['AgeIndex_0']
    return merged

def print_comparison(df, mode):
    """
    For max mode, higher age index is better. For min mode, lower age index is better.
    Print the number of tissues where algorithm 1 outperforms algorithm 0, vice versa, and ties.
    """
    if mode == "max":
        a1_wins = (df['Difference'] > 0).sum()
        a0_wins = (df['Difference'] < 0).sum()
    else:
        # In min mode, lower age index is better:
        a1_wins = (df['Difference'] < 0).sum()
        a0_wins = (df['Difference'] > 0).sum()
    ties = (df['Difference'] == 0).sum()
    
    print(f"\n{mode.upper()} Age Index Comparison:")
    print(f"Algorithm 1 wins in {a1_wins} tissues.")
    print(f"Algorithm 0 wins in {a0_wins} tissues.")
    print(f"Ties in {ties} tissues.")

def main():
    # Define folder paths
    base_path = os.path.join("data", "age_indexes", "excel_files")
    
    # File paths for algorithm 0 and algorithm 1 outputs.
    # Algorithm 1 outputs use the "best" sheet.
    a0_max_file = os.path.join(base_path, "a0_max_age_indexes.xlsx")
    a1_max_file = os.path.join(base_path, "a1_max_age_indexes.xlsx")
    a0_min_file = os.path.join(base_path, "a0_min_age_indexes.xlsx")
    a1_min_file = os.path.join(base_path, "a1_min_age_indexes.xlsx")
    
    # Load and prepare each file.
    df_a0_max = load_and_prepare(a0_max_file)
    df_a1_max = load_and_prepare(a1_max_file, sheet_name="best")
    df_a0_min = load_and_prepare(a0_min_file)
    df_a1_min = load_and_prepare(a1_min_file, sheet_name="best")
    
    # Create merged dataframes for max and min modes
    df_max = compare_mode(df_a0_max, df_a1_max, mode="max")
    df_min = compare_mode(df_a0_min, df_a1_min, mode="min")
    
    # Save the comparison as an Excel file with two sheets.
    comparison_file = os.path.join(base_path, "comparison.xlsx")
    with pd.ExcelWriter(comparison_file) as writer:
        df_min.to_excel(writer, sheet_name="min", index=False)
        df_max.to_excel(writer, sheet_name="max", index=False)
    print(f"Comparison file saved to: {comparison_file}")
    
    # Print out a summary of the results.
    print_comparison(df_max, mode="max")
    print_comparison(df_min, mode="min")

if __name__ == "__main__":
    main()
