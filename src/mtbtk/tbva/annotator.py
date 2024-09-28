import pandas as pd
from intervaltree import IntervalTree
from .SNP_effect_annotator import SnpEffectAnnotator
from .lineage_identifier import LineageIdentifier
from .config import genes

class VariantAnnotator:
    def __init__(self, SNPdf):
        self.SNPdf = SNPdf
        self.genes = genes
        self.annotator = SnpEffectAnnotator()
        self.gene_intervals = self._create_gene_intervals()

    def _create_gene_intervals(self):
        intervals = IntervalTree()
        for gene_name, gene_info in self.genes.items():
            if gene_info['biotype'] in ['protein_coding', 'pseudogene', 'other']:
                intervals[gene_info['start']:gene_info['end']] = gene_name
        return intervals

    def _is_in_same_codon(self, snp1, snp2):
        overlapping_genes = self.gene_intervals[snp1.POS] & self.gene_intervals[snp2.POS]
        for interval in overlapping_genes:
            gene_name = interval.data
            gene_info = self.genes[gene_name]
            start, end = gene_info['start'], gene_info['end']
            strand = gene_info['strand']
            
            if strand == '+':
                codon1 = (snp1.POS - start) // 3
                codon2 = (snp2.POS - start) // 3
            else:  # strand == '-'
                codon1 = (end - snp1.POS) // 3
                codon2 = (end - snp2.POS) // 3
            
            if codon1 == codon2:
                return True
        return False

    def annotate(self):
        AnnotatedSNP = []
        rows = self.SNPdf.itertuples(index=False)
        current_snp = next(rows, None)
        previous_snps = []  # 用于存储之前的SNP信息

        while current_snp is not None:
            next_snp = next(rows, None)
            
            if next_snp and self._is_in_same_codon(current_snp, next_snp):
                # 如果在同一个codon，合并SNP
                pos = current_snp.POS
                ref = current_snp.REF + next_snp.REF
                alt = current_snp.ALT + next_snp.ALT
                
                # 检查是否还有第三个SNP在同一个codon
                third_snp = next(rows, None)
                if third_snp and self._is_in_same_codon(current_snp, third_snp):
                    ref += third_snp.REF
                    alt += third_snp.ALT
                    current_snp = next(rows, None)
                else:
                    current_snp = third_snp
            else:
                pos = current_snp.POS
                ref = current_snp.REF
                alt = current_snp.ALT
                current_snp = next_snp

            # 注释
            try:
                gene_name, variant, annotation, variant_impact = self.annotator.annotate(int(pos), ref, alt)
            except Exception as e:
                print(f"Error occurred at position: {pos}")
                print(f"Current SNP: POS={pos}, REF={ref}, ALT={alt}")
                
                # 打印之前的三个SNP（如果有的话）
                for i, prev_snp in enumerate(previous_snps[-3:], 1):
                    print(f"Previous SNP {i}: POS={prev_snp['POS']}, REF={prev_snp['REF']}, ALT={prev_snp['ALT']}")
                
                print(f"Error details: {str(e)}")
                raise  # 重新抛出异常
            
            # 存储注释结果
            annotated_snp = {
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'GENE': gene_name,
                'VARIANT': variant,
                'ANNOTATION': annotation,
                'IMPACT': variant_impact
            }
            AnnotatedSNP.append(annotated_snp)
            
            # 保存当前SNP信息到previous_snps列表
            previous_snps.append(annotated_snp)
            if len(previous_snps) > 3:
                previous_snps.pop(0)  # 只保留最近的3个SNP

        return pd.DataFrame(AnnotatedSNP)
    # def annotate(self):
    #     AnnotatedSNP = []
    #     rows = self.SNPdf.itertuples(index=False)
    #     current_snp = next(rows, None)

    #     while current_snp is not None:
    #         next_snp = next(rows, None)
            
    #         if next_snp and self._is_in_same_codon(current_snp, next_snp):
    #             # 如果在同一个codon，合并SNP
    #             pos = current_snp.POS
    #             ref = current_snp.REF + next_snp.REF
    #             alt = current_snp.ALT + next_snp.ALT
                
    #             # 检查是否还有第三个SNP在同一个codon
    #             third_snp = next(rows, None)
    #             if third_snp and self._is_in_same_codon(current_snp, third_snp):
    #                 ref += third_snp.REF
    #                 alt += third_snp.ALT
    #                 current_snp = next(rows, None)
    #             else:
    #                 current_snp = third_snp
    #         else:
    #             pos = current_snp.POS
    #             ref = current_snp.REF
    #             alt = current_snp.ALT
    #             current_snp = next_snp

    #         # 注释
    #         gene_name, variant, annotation, variant_impact = self.annotator.annotate(int(pos), ref, alt)
            
    #         # 存储注释结果
    #         AnnotatedSNP.append({
    #             'POS': pos,
    #             'REF': ref,
    #             'ALT': alt,
    #             'GENE': gene_name,
    #             'VARIANT': variant,
    #             'ANNOTATION': annotation,
    #             'IMPACT': variant_impact
    #         })
    #     return pd.DataFrame(AnnotatedSNP)