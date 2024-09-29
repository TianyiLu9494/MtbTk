import pandas as pd
from .config import lineage_info
class LineageIdentifier:
    def __init__(self):        
        self.lineage_barcode = lineage_info.get('lineage_barcode')
        self.lineage_name = lineage_info.get('lineage_name')

    def get_match_lineage(self,SNP:pd.DataFrame):
        match_lineage = []
        if "931123" not in SNP['POS'].values:
            match_lineage.append({"allele":"931123_T","lineage":'lineage4','lineage_name':"Euro-American"})
            if "1759252" not in SNP['POS'].values:
                match_lineage.append({"allele":"1759252_G","lineage":'lineage4.9',"lineage_name":"Euro-American (H37Rv-like)"})
        for _,row in SNP.iterrows():
            if row['POS'] in self.lineage_barcode.keys():
                if row['ALT'] == self.lineage_barcode[row['POS']]['allele']:
                    lineage = self.lineage_barcode[row['POS']]['lineage']
                    lieage_name = self.lineage_name[lineage]
                    match_lineage.append({"allele":f"{row['POS']}_{row['ALT']}","lineage":lineage,";lineage_name":lieage_name})
        return match_lineage

    # @staticmethod
    # def check_lineage(match_lineage:list):
    #     identified_lineage = match_lineage[0]
    #     for lineage in match_lineage:
    #         if len(lineage) > len(identified_lineage):
    #             # 检查是否是子谱系
    #             if lineage.startswith(identified_lineage):
    #                 identified_lineage = lineage
    #             else:
    #                 return 'conflict_lineages'
    #         else:
    #             # 检查是否是父谱系
    #             if identified_lineage.startswith(lineage):
    #                 pass
    #             else:
    #                 return 'conflict_lineages'
    #     return identified_lineage
    
    # def identify_lineage(self, SNP:pd.DataFrame):
    #     match_lineage = self.get_match_lineage(SNP)
    #     if len(match_lineage) == 0:
    #         return 'Unknown','Unknown','Unknown','Unknown'
    #     elif len(match_lineage) == 1:
    #         sub_lineage = match_lineage[0]
    #         maj_lineage = sub_lineage.split('.')[0]
    #         sub_lineage_name = self.lineage_name[sub_lineage]
    #         maj_lineage_name = self.lineage_name[maj_lineage]
    #     else:
    #         sub_lineage =  self.check_lineage(match_lineage)
    #         if sub_lineage == 'conflict_lineages':
    #             return 'Conflict','Conflict','Conflict','Conflict'
    #         maj_lineage = sub_lineage.split('.')[0]
    #         sub_lineage_name = self.lineage_name[sub_lineage]
    #         maj_lineage_name = self.lineage_name[maj_lineage]
    #     return sub_lineage, sub_lineage_name, maj_lineage, maj_lineage_name