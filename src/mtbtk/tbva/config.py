import json
from mtbtk import DATA_DIR
def load_genes():
    genes_file = DATA_DIR / 'genes.json'
    with genes_file.open('r') as f:
        return json.load(f)

genes = load_genes()