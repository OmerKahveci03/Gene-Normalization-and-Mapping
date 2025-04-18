#!/usr/bin/env python
"""
webscrape.py

- Current directory: src/
- Age index data directory: data/age_indexes
- Input Excel files: max_age_indexes.xlsx & min_age_indexes.xlsx (located in data/age_indexes)
- Output Excel file: pubmed_genes.xlsx (written to data/age_indexes)

For each gene (from columns 3 onward in each row), this script:
  1. Collects gene names (duplicates removed).
  2. Searches PubMed using the query "{gene name} aging" with a review filter.
  3. Extracts the number of results and, if any are found, the link to the top article.
  4. Writes the gene, result count, and top article link to the first sheet (Sheet1) 
     of the output Excel file.
  5. Additionally, creates separate sheets for tissue-level summary for each input file:
     - "max_tissues" for data from max_age_indexes.xlsx
     - "min_tissues" for data from min_age_indexes.xlsx
     
     Each tissue-level sheet includes the following columns:
       * Tissue Name
       * # of Age Associated Genes
       * # of Non-Age Associated Genes
       * Gene Set from Pubmed
       * Gene Set Not in Pubmed

Note: It is assumed that:
  - Each row in the input files represents a tissue.
  - The tissue name is in the first column.
  - Gene symbols start from the third column (i.e. index 2 onward).
"""

import os
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote_plus

def get_gene_set(file_paths):
    """
    Reads all input files and extracts the set of unique gene symbols.
    """
    gene_set = set()
    for file_path in file_paths:
        print(f"Reading file: {file_path}")
        df = pd.read_excel(file_path, header=None)
        for _, row in df.iterrows():
            # Columns 0 and 1 are not gene symbols; genes start at index 2.
            for gene in row[2:]:
                if pd.notnull(gene):
                    gene_set.add(str(gene).strip())
    return gene_set

def search_pubmed(gene):
    """
    Searches PubMed for the gene with the query "{gene} aging".
    Returns the number of results and, if found, the link to the top article.
    """
    query = quote_plus(f"{gene} aging")
    url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}&filter=pubt.review"
    print(f"Searching PubMed for gene: {gene}")
    
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0 Safari/537.36"
    }
    
    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Request failed for {gene}: {e}")
        return 0, ""
    
    soup = BeautifulSoup(response.text, "html.parser")

    # Extract the number of results from the element with class "value".
    count_element = soup.find("span", class_="value")
    if count_element and count_element.text:
        try:
            count = int(count_element.text.replace(",", ""))
        except ValueError:
            count = 0
    else:
        count = 0

    # If results are found, extract the link to the top article.
    link = ""
    if count > 0:
        a_tag = soup.find("a", class_="docsum-title")
        if a_tag:
            href = a_tag.get("href")
            if href:
                link = "https://pubmed.ncbi.nlm.nih.gov" + href

    # Pause to be polite to the server.
    time.sleep(1)
    return count, link

def main():
    base_path = os.path.join(os.path.dirname(__file__), "..", "data", "age_indexes")
    input_files = [
        os.path.join(base_path, "max_age_indexes.xlsx"),
        os.path.join(base_path, "min_age_indexes.xlsx")
    ]
    
    # 1. Get the unique gene set from all input files.
    gene_set = get_gene_set(input_files)
    print(f"Total unique genes found: {len(gene_set)}")
    
    gene_results = []
    # Dictionary for quick lookup when processing tissue-level info.
    gene_info = {}
    
    # 2. Search PubMed for each gene.
    for gene in gene_set:
        count, top_link = search_pubmed(gene)
        gene_results.append({
            "Gene Name": gene,
            "number of results": count,
            "link to article": top_link if count > 0 else ""
        })
        gene_info[gene] = {"count": count, "link": top_link}
        
    # Create DataFrame for gene-level results (this will be Sheet1).
    gene_df = pd.DataFrame(gene_results)
    gene_df.sort_values(by="Gene Name", inplace=True)
    
    # 3. Process each input file separately to create tissue-level summaries.
    max_tissue_rows = []
    min_tissue_rows = []
    
    for file_path in input_files:
        print(f"Processing tissue data from: {file_path}")
        df = pd.read_excel(file_path, header=None)
        for index, row in df.iterrows():
            # Tissue name is assumed to be in the first column.
            tissue_name = row[0] if pd.notnull(row[0]) else f"Tissue_{index}"
            # Genes are from the third column onward (index 2).
            genes = [str(gene).strip() for gene in row[2:] if pd.notnull(gene)]
            age_associated_genes = []
            non_age_associated_genes = []
            for gene in genes:
                if gene in gene_info and gene_info[gene]["count"] > 0:
                    age_associated_genes.append(gene)
                else:
                    non_age_associated_genes.append(gene)
                    
            summary = {
                "Tissue Name": tissue_name,
                "# of Age Associated Genes": len(age_associated_genes),
                "# of Non-Age Associated Genes": len(non_age_associated_genes),
                "Gene Set from Pubmed": ", ".join(age_associated_genes),
                "Gene Set Not in Pubmed": ", ".join(non_age_associated_genes)
            }
            
            # Determine which file this row came from.
            if "min_age_indexes" in os.path.basename(file_path):
                min_tissue_rows.append(summary)
            elif "max_age_indexes" in os.path.basename(file_path):
                max_tissue_rows.append(summary)
    
    # Create DataFrames for each tissue-level summary.
    min_tissue_df = pd.DataFrame(min_tissue_rows)
    max_tissue_df = pd.DataFrame(max_tissue_rows)
    
    # 4. Write all DataFrames into separate sheets in the same Excel file.
    output_path = os.path.join(base_path, "pubmed_genes.xlsx")
    with pd.ExcelWriter(output_path) as writer:
        gene_df.to_excel(writer, sheet_name="Sheet1", index=False)
        min_tissue_df.to_excel(writer, sheet_name="min_tissues", index=False)
        max_tissue_df.to_excel(writer, sheet_name="max_tissues", index=False)
    
    print(f"Done. Results saved to {output_path}")

if __name__ == "__main__":
    main()
