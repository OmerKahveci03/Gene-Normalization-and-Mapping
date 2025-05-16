#!/usr/bin/env python
"""
webscrape.py

- Located in src/age_indexing
- Reads 'best' sheets from Excel files in data/age_indexes/excel_files:
    * a1_max_age_indexes.xlsx
    * a1_min_age_indexes.xlsx
- Writes PubMed summary Excel to data/age_indexes/pubmed_genes.xlsx
- Generates bar charts in visuals/:
    * max_pubmed_age.png
    * min_pubmed_age.png
"""
import os
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote_plus
import matplotlib.pyplot as plt


def get_gene_set_from_best(files):
    """
    Reads the 'best' sheet of each file, collects unique genes from the 'Gene_Set' column.
    """
    genes = set()
    for fp in files:
        print(f"Reading 'best' sheet: {fp}")
        df = pd.read_excel(fp, sheet_name='best', usecols=['Gene_Set'])
        for cell in df['Gene_Set'].dropna():
            for gene in str(cell).split(','):
                g = gene.strip()
                if g:
                    genes.add(g)
    return genes


def search_pubmed(gene):
    """
    Searches PubMed for "{gene} aging" (review filter).
    Returns count and top article link.
    """
    query = quote_plus(f"{gene} aging")
    url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}&filter=pubt.review"
    print(f"Searching PubMed for gene: {gene}")
    headers = {"User-Agent": "Mozilla/5.0"}
    try:
        resp = requests.get(url, headers=headers, timeout=10)
        resp.raise_for_status()
    except requests.RequestException as e:
        print(f"Request failed for {gene}: {e}")
        return 0, ""
    soup = BeautifulSoup(resp.text, "html.parser")
    span = soup.find("span", class_="value")
    if span and span.text:
        try:
            count = int(span.text.replace(",", ""))
        except ValueError:
            count = 0
    else:
        count = 0
    link = ""
    if count > 0:
        a = soup.find("a", class_="docsum-title")
        if a and a.get("href"):
            link = "https://pubmed.ncbi.nlm.nih.gov" + a["href"]
    time.sleep(1)
    return count, link


def plot_pubmed_summary(df, kind, output_path):
    """
    Creates a bar chart: green bars = total genes, red bars = non-age-associated.
    """
    tissues = df['Tissue Name']
    age_assoc = df['# of Age Associated Genes']
    non_age = df['# of Non-Age Associated Genes']
    total = age_assoc + non_age

    x = range(len(tissues))
    plt.figure(figsize=(10, 6))
    plt.bar(x, total, color='green', label='Total Genes')
    plt.bar(x, non_age, color='red', label='Non-Age-Associated')
    plt.xticks(x, tissues, rotation=45, ha='right')
    plt.xlabel('Tissue')
    plt.ylabel('# of Genes')
    plt.title(f'{kind.title()} PubMed Aging Gene Summary')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    # Base paths
    script_dir = os.path.dirname(__file__)
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
    excel_dir = os.path.join(project_root, 'data', 'age_indexes', 'excel_files')
    files = [
        os.path.join(excel_dir, 'a1_max_age_indexes.xlsx'),
        os.path.join(excel_dir, 'a1_min_age_indexes.xlsx')
    ]

    # 1. Collect genes from 'best' sheets only
    gene_set = get_gene_set_from_best(files)
    print(f"Total unique genes to query: {len(gene_set)}")

    # 2. PubMed search per gene
    gene_results, gene_info = [], {}
    for gene in sorted(gene_set):
        cnt, link = search_pubmed(gene)
        gene_results.append({'Gene Name': gene, 'number of results': cnt, 'link to article': link})
        gene_info[gene] = cnt

    gene_df = pd.DataFrame(gene_results).sort_values('Gene Name')

    # 3. Summarize per-tissue from 'best' sheets
    max_rows, min_rows = [], []
    for fp in files:
        print(f"Summarizing tissue-level for: {fp}")
        df = pd.read_excel(fp, sheet_name='best', usecols=['Tissue', 'Gene_Set'])
        for _, row in df.iterrows():
            tissue = row['Tissue']
            genes = [g.strip() for g in str(row['Gene_Set']).split(',') if g.strip()]
            age_assoc = [g for g in genes if gene_info.get(g,0)>0]
            non_age = [g for g in genes if gene_info.get(g,0)==0]
            rec = {
                'Tissue Name': tissue,
                '# of Age Associated Genes': len(age_assoc),
                '# of Non-Age Associated Genes': len(non_age),
                'Gene Set from Pubmed': ', '.join(age_assoc),
                'Gene Set Not in Pubmed': ', '.join(non_age)
            }
            if os.path.basename(fp).startswith('a1_max'): max_rows.append(rec)
            else: min_rows.append(rec)

    max_df, min_df = pd.DataFrame(max_rows), pd.DataFrame(min_rows)

    # 4. Write results
    out_excel = os.path.join(project_root, 'data', 'age_indexes', 'pubmed_genes.xlsx')
    with pd.ExcelWriter(out_excel) as w:
        gene_df.to_excel(w, sheet_name='Sheet1', index=False)
        min_df.to_excel(w, sheet_name='min_tissues', index=False)
        max_df.to_excel(w, sheet_name='max_tissues', index=False)
    print(f"PubMed results written to {out_excel}")

    # 5. Create bar charts
    visuals = os.path.join(project_root, 'visuals')
    os.makedirs(visuals, exist_ok=True)
    plot_pubmed_summary(max_df, 'max', os.path.join(visuals, 'max_pubmed_age.png'))
    plot_pubmed_summary(min_df, 'min', os.path.join(visuals, 'min_pubmed_age.png'))
    print(f"Charts saved in {visuals}")

if __name__=='__main__':
    main()
