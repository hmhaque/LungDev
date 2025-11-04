import bone
pd = bone.pd
re = bone.re
hu = bone.hu
np = bone.np

def getJakubzick2017(self, tn=1):
    self.prepareData("T121", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type (ch1)")
    atypes = ['IM', 'AM']
    ahash = {'IM1':0, 'IM2': 0, 'IM3': 0, 'AM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJakubzick2017 = getJakubzick2017

def getStefan2023(self, tn=1):
    self.prepareData("T122", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type (ch1)")
    atypes = ['IM', 'AM']
    ahash = {'IM':0, 'AM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getStefan2023 = getStefan2023

def getSandra2024(self, tn=1):
    self.prepareData("T123", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type (ch1)")
    atypes = ['IM', 'AM']
    ahash = {'IM':0, 'AM':1}
    
    if (tn == 2):
        atypes = ['Monocytes', 'IM', 'AM']
        ahash = {'Monocytes':0, 'IM':1, 'AM':2}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSandra2024 = getSandra2024

def getSandra_IPT_2024(self, tn=1):
    self.prepareData("T124", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['IM', 'AM']
    ahash = {'IM':0, 'AM':1}
    
    if (tn == 2):
        atypes = ['Monocytes', 'AM', 'IM']
        ahash = {'Monocytes':0, 'AM':1, 'IM':2}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSandra_IPT_2024 = getSandra_IPT_2024


# for deep learning prediction output
def getSandra2024_prediction(self, tn=1):
    self.prepareData("T123.11", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type (ch1)")
    atypes = ['IM', 'AM']
    ahash = {'IM':0, 'AM':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSandra2024_prediction = getSandra2024_prediction

def getSandra_IPT_2024_prediction(self, tn=1):
    self.prepareData("T124.11", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['IM', 'AM']
    ahash = {'IM':0, 'AM':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSandra_IPT_2024_prediction = getSandra_IPT_2024_prediction


# previously processed 

def getSajti2020(self, tn=1):
    self.prepareData("MACV96")
    h = self.h
    atype = h.getSurvName('c cell type')
    atypes = [ 'alveolar macrophages (AM)', 'interstitial macrophages (IM)', 'inflammatory monocytes (iMo)']
    ahash = {}
    if (tn == 2):
        atypes = ['IM', 'AM']
        ahash = {'interstitial macrophages (IM)':0, 'alveolar macrophages (AM)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSajti2020 = getSajti2020

def getSajti2020BL6DBA(self, tn=1):
    self.prepareData("MACV96.3")
    h = self.h
    atype = h.getSurvName('c src1')
    atypes = ['alveolar macrophages', 'interstital macrophages', 'inflammatory macrophages', 'patrolling macrophages']
    ahash = {}
    if (tn == 2):
        atypes = ['IM', 'AM']
        ahash = {'interstital macrophages':0, 'alveolar macrophages':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSajti2020BL6DBA = getSajti2020BL6DBA

def getVanneste2023(self, tn=1):
    self.prepareData("BL34")
    h = self.h
    atype = h.getSurvName('c src1')
    atypes = ['Ly6C+ Mo',
            'CD206+ control IM',
            'AM',
            'CD206+ refilled IM',
            'CD206- refilled IM',
            'CD206- control IM']
    ahash = {}
    if (tn == 2):
        atype = h.getSurvName('c cell type 2')
        atypes = ['IM', 'AM']
        ahash = {'IM':0, 'AM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanneste2023 = getVanneste2023

def getQiang2021_sc(self, tn=1):
    self.prepareData("T132", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c patient status (ch1)")
    atypes = ['Healthy', 'COPD']
    ahash = { 'Healthy no smoker': 0, 'Healthy smoker': 0,
             'COPD': 1}
    
    if (tn == 2):
        atypes = ['Healthy', 'COPD']
        ahash = { 'Healthy no smoker': 0,'COPD': 1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQiang2021_sc = getQiang2021_sc

def getEllen2022(self, tn=1):
    self.prepareData("T126", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Non-smoker', 'Smoker']
    ahash = {'Non-smoker':0, 'Smoker':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEllen2022 = getEllen2022


def getQiang2022(self, tn=1):
    self.prepareData("T127", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Smoker', 'Non-smoker']
    ahash = { 'smoker_monocyte': 0, 'smoker_macrophage': 0,
             'non_smoker_monocyte': 1,
             'non_smoker_macrophage': 1}
    
    if (tn == 2):
        atypes = ['Non-smoker', 'Smoker']
        ahash = { 'non_smoker_macrophage': 0, 'smoker_macrophage': 1}
        
    if (tn == 3):
        atypes = ['Non-smoker', 'Smoker']
        ahash = {'non_smoker_monocyte': 0, 'smoker_monocyte': 1}
     
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQiang2022 = getQiang2022

def getThomas2019(self, tn=1):
    self.prepareData("T128", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Non-smoker', 'Smoker']
    ahash = {'non-smoker':0, 'smoker':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getThomas2019 = getThomas2019

def getYael2019(self, tn=1):
    self.prepareData("T129", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Non-smoker', 'Smoker']
    ahash = {'non-smoker':0, 'smoker':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYael2019 = getYael2019

def getYael2016(self, tn=1):
    self.prepareData("T130", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Non-smoker', 'Smoker']
    ahash = {'non-smoker':0, 'smoker':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYael2016 = getYael2016

def getYaelAll2019(self, tn=1):
    self.prepareData("T131", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Non-smoker', 'Smoker']
    ahash = {'non-smoker':0, 'smoker':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYaelAll2019 = getYaelAll2019

def getQiang2022scRNAseq(self, tn=1):
    self.prepareData("T132", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Smoker', 'Non-smoker', 'COPD']
    ahash = { 'smoker': 0, 'non-smoker': 1, 'COPD': 2}
    
    if (tn == 2):
        atypes = ['Smoker', 'Non-smoker']
        ahash = { 'smoker': 0, 'non-smoker': 1}

    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQiang2022scRNAseq = getQiang2022scRNAseq


# previously processed
def getOBeirne2020(self, tn=1):
    self.prepareData("LU32")
    h = self.h
    atype = h.getSurvName('c phenotype')
    atypes = ['S', 'NS', 'COPDs']
    ahash = {}
    if (tn == 2):
        atypes = ['NS', 'S']
        ahash = {'NS':0, 'S':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOBeirne2020 = getOBeirne2020

def getWoodruff2005mac(self, tn=1):
    self.prepareData("MAC12")
    h = self.h
    atype = h.getSurvName('c status')
    atypes = ['S', 'NS', 'Asth']
    ahash = {}
    if (tn == 2):
        atypes = ['NS', 'S']
        ahash = {'Nonsmoker':0, 'Smoker':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWoodruff2005mac = getWoodruff2005mac

# validation 
def getDeepak2020(self, tn=1):
    self.prepareData("T133", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['WT', 'KO']
    ahash = {'WT':0, 'KO':1}
    
    if (tn == 2):
        atypes = ['WT', 'Pathogen_free']
        ahash = {'WT':0, 'KO':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeepak2020 = getDeepak2020

def getSahar2023(self, tn=1):
    self.prepareData("T134", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type (hr)")
    atypes = ['Long', 'Short']
    ahash = {'5h':0, '0h':1}
    
    if (tn == 2):
        atypes = ['Long', 'Short']
        ahash = {'5h':0, '0h':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSahar2023 = getSahar2023

def getReza2020(self, tn=1): # 
    self.prepareData("T135", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type (hr)")
    atypes = ['Long', 'Short']
    ahash = {'24h':0, '0h':1, '1h':1, '2h':1, '4h':1, '8h':1}
    
    if (tn == 2):
        atypes = ['Long', 'Short']
        ahash = {'24h':0, '0h':1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getReza2020 = getReza2020


def getYushi2025(self, tn=1): 
    self.prepareData("T136", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['IFNAR-/-PBS', 'WT']
    ahash = {'IFNAR-/-PBS':0, 'IFNAR-/-IAV':0, 'WT-PBS':1, 'WT-IAV':1}
    
    if (tn == 2):
        atypes = ['IFNAR-/-PBS', 'WT-PBS']
        ahash = {'IFNAR-/-PBS':0, 'WT-PBS':1}
    if (tn == 3):
        atype = self.h.getSurvName("c type new")
        atypes = ['IFNAR', 'WT']
        ahash = {'IFNAR KO_PBS': 0, 'IFNAR KO_IAV': 0, 'IFNAR KO_B16':0,
                    'WT_SP': 1, 'WT_B16': 1}
    if (tn == 4):
        atype = self.h.getSurvName("c genotype (ch1)")
        atypes = ['IFNAR KO', 'WT']
        ahash = {'IFNAR KO':0 , 'WT': 1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYushi2025 = getYushi2025

def getQiang2021_sc(self, tn=1):
    self.prepareData("T132", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c patient status (ch1)")
    atypes = ['Healthy', 'COPD']
    ahash = { 'Healthy no smoker': 0, 'Healthy smoker': 0,
             'COPD': 1}
    
    if (tn == 2):
        atypes = ['Healthy', 'COPD']
        ahash = { 'Healthy no smoker': 0,'COPD': 1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQiang2021_sc = getQiang2021_sc

def getHeo2025(self, tn=1):
    self.prepareData("LU24")
    h = self.h
    atype = h.getSurvName('c group')
    atypes = ['Control', 'COPD']
    ahash = {'Control': 0, 'COPD': 1}

    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHeo2025 = getHeo2025

def getWang2024(self, tn=1):
    self.prepareData("LU26")
    h = self.h
    atype = h.getSurvName('c copd status')
    atypes = ['Control', 'COPD']
    ahash = {'none': 0, 'COPD': 1}

    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2024 = getWang2024

def getAdams2020(self, tn=1):
    self.prepareData("LU27")
    h = self.h
    atype = h.getSurvName('c Disease_Identity')
    atypes = ['Control', 'COPD']
    ahash = {'Control': 0, 'COPD': 1}

    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAdams2020 = getAdams2020


def getMcDonough_sc(self, tn=1):
    self.prepareData("T148", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    atype = self.h.getSurvName("c type")
    atypes = ['Healthy', 'COPD']
    ahash = { 'Control': 0,'COPD': 1}
        
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMcDonough_sc = getMcDonough_sc


def getWang2020CoV2scblk(self, tn=1, tb=1):
    self.prepareData('COV384')
    atype = self.h.getSurvName("c type of death")
    atypes = ['DBD', 'DCD']
    ahash = {}
    if tn == 2:
        atype = self.h.getSurvName("c age")
        atypes = ['N', 'C', 'A']
        ahash = {'3 yr':1, '29 wkGA':0, '31 wkGA':0,
                '31 yrs':2, '29 yrs':2, '33 yrs':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2020CoV2scblk = getWang2020CoV2scblk

def getMoon2024copd(self, tn=1, ta=0):
    self.prepareData("LU25")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub(".* C", "C", str(k)) for k in atype]
    atype = [hu.re.sub(", scRNA", "", str(k)) for k in atype]
    atypes = ['Ct', 'COPDt', 'Co', 'COPDo']
    ahash = {'Control, Tissue':0, 'COPD, Tissue':1, 'Control, Organoid':2,
            'COPD, Organoid':3, 'COPD, Organid':3}
    if tn == 2:
        atypes = ['Ct', 'COPDt']
        ahash = {'Control, Tissue':0, 'COPD, Tissue':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMoon2024copd = getMoon2024copd

def getSountoulidis2023EmbryonicLungHM(self, tn=1):
    self.prepareData("T117", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")

    if tn == 1:
        atype = self.h.getSurvName("c gender (ch1)")
        atypes = ['Female', 'Male']
        ahash = {'female':0, 'male':1}
    if tn ==2:
        atype = self.h.getSurvName("c age (ch1)")
        atypes = ['Early', 'Mid', 'Late']
        ahash = {
            'PCW5': 0, 'PCW5.5': 0, 'PCW6': 0, 'PCW7': 0, 'PCW8': 0, 'PCW8.5': 0,
            'PCW10': 1, 'PCW11.5': 1, 'PCW12': 1,
            'PCW13': 2, 'PCW14': 2
        }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSountoulidis2023EmbryonicLungHM = getSountoulidis2023EmbryonicLungHM

def getNegretti2021mmLung(self, tn=1):
    self.prepareData("T118", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")
    
    if tn ==1:
        atype = self.h.getSurvName("c age (ch1)")
        atypes = ['Embryonic', 'Post Natal']
        ahash = {
            'embryonic day 12': 0, 'embryonic day 15': 0,
            'post natal day 3': 1, 'post natal day 5-0': 1, 'post natal day 5-1': 1
        }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNegretti2021mmLung = getNegretti2021mmLung

def getNegretti2021mmLung_with_epi(self, tn=1):
    self.prepareData("T119", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")
    
    if tn ==1:
        atype = self.h.getSurvName("c age (ch1)")
        atypes = ['Embryonic', 'Post Natal']
        ahash = {
            'embryonic day 12': 0, 'embryonic day 15': 0,
            'post natal day 3': 1, 'post natal day 5-0': 1, 'post natal day 5-1': 1, 'post natal day 7': 1
        }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNegretti2021mmLung_with_epi = getNegretti2021mmLung_with_epi

def getAmyWong2024fetalLung(self, tn=1):
    self.prepareData("T120", cfile="/Users/m2haque/public_html/Hegemon/explore.conf")
    
    if tn ==1:
        atype = self.h.getSurvName("c type (ch1)")
        atypes = ['Early', 'Mid', 'Late']
        ahash = {
            'Early': 0, 'Mid': 1, 'Late': 2
        }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAmyWong2024fetalLung = getAmyWong2024fetalLung