import json
from mtbtk import DATA_DIR
def load_genes():
    genes_file = DATA_DIR / 'genes.json'
    with genes_file.open('r') as f:
        return json.load(f)

def load_lineage_info():
    lineage_info_file = DATA_DIR / 'lineage_info.json'
    with lineage_info_file.open('r') as f:
        return json.load(f)

def load_domains():
    domains_file = DATA_DIR / 'domains.json'
    with domains_file.open('r') as f:
        return json.load(f)

def load_epitopes():
    epitopes_file = DATA_DIR / 'epitopes.json'
    with epitopes_file.open('r') as f:
        return json.load(f)

genes = load_genes()
lineage_info = load_lineage_info()
domains = load_domains()
epitopes = load_epitopes()
