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
  4. Writes the gene, result count, and top article link to an output Excel file.
"""

import os
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote_plus

def get_gene_set():
    base_path = os.path.join(os.path.dirname(__file__), "..", "data", "age_indexes")
    files = ["max_age_indexes.xlsx", "min_age_indexes.xlsx"]
    gene_set = set()

    for file_name in files:
        file_path = os.path.join(base_path, file_name)
        print(f"Reading file: {file_path}")
        df = pd.read_excel(file_path, header=None)
        for _, row in df.iterrows():
            for gene in row[2:]:
                if pd.notnull(gene):
                    gene_set.add(str(gene).strip())
    return gene_set

def search_pubmed(gene):
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

    # Extract the number of results.
    count_element = soup.find("span", class_="value")
    if count_element and count_element.text:
        try:
            count = int(count_element.text.replace(",", ""))
        except ValueError:
            count = 0
    else:
        count = 0

    # If results were found, try to extract the top article link.
    link = ""
    if count > 0:
        a_tag = soup.find("a", class_="docsum-title")
        if a_tag:
            href = a_tag.get("href")
            if href:
                link = "https://pubmed.ncbi.nlm.nih.gov" + href

    # Be polite and pause between requests to avoid overloading the server.
    time.sleep(1)
    return count, link

def main():
    gene_set = get_gene_set()
    print(f"Total unique genes found: {len(gene_set)}")
    
    results = []
    
    for gene in gene_set:
        count, top_link = search_pubmed(gene)
        results.append({
            "Gene Name": gene,
            "number of results": count,
            "link to article": top_link if count > 0 else ""
        })
    output_df = pd.DataFrame(results)

    # Sort alphabetically by 'Gene Name'
    output_df.sort_values(by="Gene Name", inplace=True)
    output_path = os.path.join(os.path.dirname(__file__), "..", "data", "age_indexes", "pubmed_genes.xlsx")
    output_df.to_excel(output_path, index=False)
    print(f"Done. Results saved to {output_path}")

if __name__ == "__main__":
    main()
