import pandas as pd
class LineageIdentifier:
    def __init__(self):        
        self.lineage_barcode = {
        '615938'	: {'allele':'A', 'lineage':'lineage1'},
        '4404247'	: {'allele':'A', 'lineage':'lineage1.1'},
        '3021283'	: {'allele':'A', 'lineage':'lineage1.1.1'},
        '3216553'	: {'allele':'A', 'lineage':'lineage1.1.1.1'},
        '2622402'	: {'allele':'A', 'lineage':'lineage1.1.2'},
        '1491275'	: {'allele':'A', 'lineage':'lineage1.1.3'},
        '3479545'	: {'allele':'A', 'lineage':'lineage1.2.1'},
        '3470377'	: {'allele':'T', 'lineage':'lineage1.2.2'},
        '497491'	: {'allele':'A', 'lineage':'lineage2'},
        '1881090'	: {'allele':'T', 'lineage':'lineage2.1'},
        '2505085'	: {'allele':'A', 'lineage':'lineage2.2'},
        '797736'	: {'allele':'T', 'lineage':'lineage2.2.1'},
        '4248115'	: {'allele':'T', 'lineage':'lineage2.2.1.1'},
        '3836274'	: {'allele':'A', 'lineage':'lineage2.2.1.2'},
        '346693'	: {'allele':'T', 'lineage':'lineage2.2.2'},
        '3273107'	: {'allele':'A', 'lineage':'lineage3'},
        '1084911'	: {'allele':'A', 'lineage':'lineage3.1.1'},
        '3722702'	: {'allele':'C', 'lineage':'lineage3.1.2'},
        '1237818'	: {'allele':'G', 'lineage':'lineage3.1.2.1'},
        '2874344'	: {'allele':'A', 'lineage':'lineage3.1.2.2'},
        '62657' : {'allele':'A', 'lineage':'lineage4.1'},
        '514245'	: {'allele':'T', 'lineage':'lineage4.1.1'},
        '1850119'	: {'allele':'T', 'lineage':'lineage4.1.1.1'},
        '541048'	: {'allele':'G', 'lineage':'lineage4.1.1.2'},
        '4229087'	: {'allele':'T', 'lineage':'lineage4.1.1.3'},
        '891756'	: {'allele':'G', 'lineage':'lineage4.1.2'},
        '107794'	: {'allele':'T', 'lineage':'lineage4.1.2.1'},
        '2411730'	: {'allele':'C', 'lineage':'lineage4.2'},
        '783601'	: {'allele':'C', 'lineage':'lineage4.2.1'},
        '1487796'	: {'allele':'A', 'lineage':'lineage4.2.2'},
        '1455780'	: {'allele':'C', 'lineage':'lineage4.2.2.1'},
        '764995'	: {'allele':'G', 'lineage':'lineage4.3'},
        '615614'	: {'allele':'A', 'lineage':'lineage4.3.1'},
        '4316114'	: {'allele':'A', 'lineage':'lineage4.3.2'},
        '3388166'	: {'allele':'G', 'lineage':'lineage4.3.2.1'},
        '403364'	: {'allele':'A', 'lineage':'lineage4.3.3'},
        '3977226'	: {'allele':'A', 'lineage':'lineage4.3.4'},
        '4398141'	: {'allele':'A', 'lineage':'lineage4.3.4.1'},
        '1132368'	: {'allele':'T', 'lineage':'lineage4.3.4.2'},
        '1502120'	: {'allele':'A', 'lineage':'lineage4.3.4.2.1'},
        '4307886'	: {'allele':'A', 'lineage':'lineage4.4'},
        '4151558'	: {'allele':'A', 'lineage':'lineage4.4.1'},
        '355181'	: {'allele':'A', 'lineage':'lineage4.4.1.1'},
        '2694560'	: {'allele':'C', 'lineage':'lineage4.4.1.2'},
        '4246508'	: {'allele':'A', 'lineage':'lineage4.4.2'},
        '1719757'	: {'allele':'T', 'lineage':'lineage4.5'},
        '3466426'	: {'allele':'A', 'lineage':'lineage4.6'},
        '4260268'	: {'allele':'C', 'lineage':'lineage4.6.1'},
        '874787'	: {'allele':'A', 'lineage':'lineage4.6.1.1'},
        '1501468'	: {'allele':'C', 'lineage':'lineage4.6.1.2'},
        '4125058'	: {'allele':'C', 'lineage':'lineage4.6.2'},
        '3570528'	: {'allele':'G', 'lineage':'lineage4.6.2.1'},
        '2875883'	: {'allele':'T', 'lineage':'lineage4.6.2.2'},
        '4249732'	: {'allele':'G', 'lineage':'lineage4.7'},
        '3836739'	: {'allele':'A', 'lineage':'lineage4.8'},
        '1799921'	: {'allele':'A', 'lineage':'lineage5'},
        '1816587'	: {'allele':'G', 'lineage':'lineage6'},
        '1137518'	: {'allele':'A', 'lineage':'lineage7'},
        '2831482'	: {'allele':'G', 'lineage':'lineageBOV'},
        '1882180'	: {'allele':'T', 'lineage':'lineageBOV_AFRI'}
        }
        self.lineage_name = {
        'lineage1': 'Indo-Oceanic',
        'lineage1.1': 'Indo-Oceanic',
        'lineage1.1.1': 'Indo-Oceanic',
        'lineage1.1.1.1': 'Indo-Oceanic',
        'lineage1.1.2': 'Indo-Oceanic',
        'lineage1.1.3': 'Indo-Oceanic',
        'lineage1.2': 'Indo-Oceanic',
        'lineage1.2.1': 'Indo-Oceanic',
        'lineage1.2.2': 'Indo-Oceanic',
        'lineage2': 'East-Asian',
        'lineage2.1': 'East-Asian (non-Beijing)',
        'lineage2.2': 'East-Asian (Beijing)',
        'lineage2.2.1': 'East-Asian',
        'lineage2.2.1.1': 'East-Asian',
        'lineage2.2.1.2': 'East-Asian',
        'lineage2.2.2': 'East-Asian',
        'lineage3': 'East-African-Indian',
        'lineage3.1': 'East-African-Indian',
        'lineage3.1.1': 'East-African-Indian',
        'lineage3.1.2': 'East-African-Indian',
        'lineage3.1.2.1': 'East-African-Indian',
        'lineage3.1.2.2': 'East-African-Indian',
        'lineage4': 'Euro-American',
        'lineage4.1': 'Euro-American',
        'lineage4.1.1': 'Euro-American (X-type)',
        'lineage4.1.1.1': 'Euro-American (X-type)',
        'lineage4.1.1.2': 'Euro-American (X-type)',
        'lineage4.1.1.3': 'Euro-American (X-type)',
        'lineage4.1.2': 'Euro-American',
        'lineage4.1.2.1': 'Euro-American (Haarlem)',
        'lineage4.2': 'Euro-American',
        'lineage4.2.1': 'Euro-American (Ural)',
        'lineage4.2.2': 'Euro-American',
        'lineage4.2.2.1': 'Euro-American (TUR)',
        'lineage4.3': 'Euro-American (LAM)',
        'lineage4.3.1': 'Euro-American (LAM)',
        'lineage4.3.2': 'Euro-American (LAM)',
        'lineage4.3.2.1': 'Euro-American (LAM)',
        'lineage4.3.3': 'Euro-American (LAM)',
        'lineage4.3.4': 'Euro-American (LAM)',
        'lineage4.3.4.1': 'Euro-American (LAM)',
        'lineage4.3.4.2': 'Euro-American (LAM)',
        'lineage4.3.4.2.1': 'Euro-American (LAM)',
        'lineage4.4': 'Euro-American',
        'lineage4.4.1': 'Euro-American',
        'lineage4.4.1.1': 'Euro-American (S-type)',
        'lineage4.4.1.2': 'Euro-American',
        'lineage4.4.2': 'Euro-American',
        'lineage4.5': 'Euro-American',
        'lineage4.6': 'Euro-American',
        'lineage4.6.1': 'Euro-American (Uganda)',
        'lineage4.6.1.1': 'Euro-American',
        'lineage4.6.1.2': 'Euro-American',
        'lineage4.6.2': 'Euro-American',
        'lineage4.6.2.1': 'Euro-American',
        'lineage4.6.2.2': 'Euro-American (Cameroon)',
        'lineage4.7': 'Euro-American (mainly T)',
        'lineage4.8': 'Euro-American (mainly T)',
        'lineage4.9': 'Euro-American (H37Rv-like)',
        'lineage5': 'West-Africa 1',
        'lineage6': 'West-Africa 2',
        'lineageBOV': 'M. bovis',
        'lineageBOV_AFRI': 'M. bovis and West-Africa 2',
        'lineage7': 'Lineage 7'
        }

    def get_match_lineage(self,SNP:pd.DataFrame):
        match_lineage = []
        if "931123" not in SNP['POS'].values:
            match_lineage.append('lineage4')
            if "1759252" not in SNP['POS'].values:
                match_lineage.append('lineage4.9')
        for _,row in SNP.iterrows():
            if row['POS'] in self.lineage_barcode.keys():
                if row['ALT'] == self.lineage_barcode[row['POS']]['allele']:
                    match_lineage.append(self.lineage_barcode[row['POS']]['lineage'])
        return match_lineage

    @staticmethod
    def check_lineage(match_lineage:list):
        identified_lineage = match_lineage[0]
        for lineage in match_lineage:
            if len(lineage) > len(identified_lineage):
                # 检查是否是子谱系
                if lineage.startswith(identified_lineage):
                    identified_lineage = lineage
                else:
                    return 'conflict_lineages'
            else:
                # 检查是否是父谱系
                if identified_lineage.startswith(lineage):
                    pass
                else:
                    return 'conflict_lineages'
        return identified_lineage
    
    def identify_lineage(self, SNP:pd.DataFrame):
        match_lineage = self.get_match_lineage(SNP)
        if len(match_lineage) == 0:
            return 'Unknown','Unknown','Unknown','Unknown'
        elif len(match_lineage) == 1:
            sub_lineage = match_lineage[0]
            maj_lineage = sub_lineage.split('.')[0]
            sub_lineage_name = self.lineage_name[sub_lineage]
            maj_lineage_name = self.lineage_name[maj_lineage]
        else:
            sub_lineage =  self.check_lineage(match_lineage)
            if sub_lineage == 'conflict_lineages':
                return 'Conflict','Conflict','Conflict','Conflict'
            maj_lineage = sub_lineage.split('.')[0]
            sub_lineage_name = self.lineage_name[sub_lineage]
            maj_lineage_name = self.lineage_name[maj_lineage]
        return sub_lineage, sub_lineage_name, maj_lineage, maj_lineage_name