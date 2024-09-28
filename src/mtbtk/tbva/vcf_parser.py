import pandas as pd
class VCFParser:
    def __init__(self, vcf_file, filter_pass=True):
        self.SNP, self.INDEL = self._categorize_variants(self._parse_vcf(vcf_file, filter_pass))
        self._all_variants = None  # 延迟加载

    def _parse_vcf(self, vcf_file, filter_pass):
        # Read the VCF file, skipping lines starting with ##
        with open(vcf_file, 'r') as file:
            lines = [line.strip() for line in file if not line.startswith('##')]

        # Parse the header and data
        header = lines[0].lstrip('#').split('\t')
        data = [line.split('\t') for line in lines[1:]]

        # Create a DataFrame
        df = pd.DataFrame(data, columns=header)

        # Keep only POS, REF, ALT columns, and FILTER if needed
        columns_to_keep = ['POS', 'REF', 'ALT']
        if filter_pass:
            columns_to_keep.append('FILTER')

        df = df[columns_to_keep]

        # Convert POS to integer
        df['POS'] = df['POS'].astype(int)

        # Filter PASS variants if required
        if filter_pass:
            df = df[df['FILTER'] == 'PASS'].drop('FILTER', axis=1)

        return df

    def _categorize_variants(self, df):
        # Categorize variants into SNPs and INDELs
        df['is_snp'] = df.apply(lambda row: len(row['REF']) == len(row['ALT']), axis=1)
        snp = df[df['is_snp']].drop('is_snp', axis=1)
        indel = df[~df['is_snp']].drop('is_snp', axis=1)
        return snp, indel

    def get_snps(self):
        return self.SNP

    def get_indels(self):
        return self.INDEL

    def get_all_variants(self):
        if self._all_variants is None:
            self._all_variants = pd.concat([self.SNP, self.INDEL]).sort_values('POS').reset_index(drop=True)
        return self._all_variants