import json
from mtbtk import DATA_DIR
def load_GOs():
    genes_file = DATA_DIR / 'GO.json'
    with genes_file.open('r') as f:
        return json.load(f)
    
def load_KEGG_pathways():
    genes_file = DATA_DIR / 'KEGG_pathways.json'
    with genes_file.open('r') as f:
        return json.load(f)

def load_TF_genes():
    genes_file = DATA_DIR / 'TF_genes.json'
    with genes_file.open('r') as f:
        return json.load(f)
    
GOs = load_GOs()
KEGG_pathways = load_KEGG_pathways()
TF_genes = load_TF_genes()
