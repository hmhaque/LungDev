import bone
pd = bone.pd
re = bone.re
hu = bone.hu
np = bone.np

def plotTitleBar(cval, atypes, params):
    dpi = 100
    if 'dpi' in params:
        dpi = params['dpi']
    w,h = (5, 0.8)
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
    color_sch1 = ["#3B449C", "#B2509E","#EA4824"]
    color_sch1 = ["#00CC00", "#EFF51A","#EC008C", "#F7941D", "#808285",
            'cyan', 'blue', 'black', 'green', 'red']
    if 'acolor' in params:
        color_sch1 = params['acolor']
    if 'cval' in params:
        cval = params['cval']

    ax = None
    if 'ax' in params:
        ax = params['ax']
    if ax is None:
        fig = bone.plt.figure(figsize=(w,h), dpi=dpi)
        ax = fig.add_subplot(1, 1, 1)
    nAt = len(cval[0])
    extent = [0, nAt, 0, 5]
    ax.axis(extent)
    cmap = bone.colors.ListedColormap(color_sch1)
    boundaries = range(len(color_sch1) + 1)
    norm = bone.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    #ax.imshow(cval, interpolation='none', cmap=cmap, \
    #                  norm=norm, extent=extent, aspect="auto")
    y = [0, 5]
    x = bone.np.arange(nAt + 1)
    ax.pcolormesh(x, y, cval, cmap=cmap, norm=norm, zorder=-1.0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(top=False, left=False, bottom=False, right=False)
    ax.set_xticks(bone.np.arange(0, nAt, 1))
    ax.grid(which='major', alpha=0.2, linestyle='-', linewidth=0.5,
            color='black', zorder=2.0)
    for edge, spine in ax.spines.items():
                spine.set_visible(False)
    divider = bone.make_axes_locatable(ax)
    width = bone.axes_size.AxesX(ax, aspect=1./20)
    spaceAnn = 70
    widthAnn = 3
    tAnn = 1
    if 'spaceAnn' in params:
        spaceAnn = params['spaceAnn']
    if 'widthAnn' in params:
        widthAnn = params['widthAnn']
    if 'tAnn' in params:
        tAnn = params['tAnn']
    pad = bone.axes_size.Fraction(0.1, width)
    lax = divider.append_axes("top", size="100%", pad="20%", frame_on=False)
    lax.axison = False
    lax.axis(extent)
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.grid(False)
    lax.tick_params(top=False, left=False, bottom=False, right=False)
    if 'atypes' in params:
        atypes = params['atypes']
    bone.barTop(lax, atypes, color_sch1, params)
    return ax

bone.plotTitleBar = plotTitleBar

def adjustFigure(ana, fig, xlabel="Score", wfactor=1):
    rocauc = ana.getROCAUC()
    ax = fig.axes[1]
    #ax.set_xlim([-5, 10])
    ax.set_title(f"ROC-AUC = {rocauc}", fontsize=16)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=20)
    #fig.patch.set_facecolor('#FFF7CD')
    #ax.set_facecolor('#FFF7CD')
    for tobj in ax.texts:
        tobj.set_fontsize(16)
        #tobj.set_color('black')
    #for tobj in fig.findobj(plt.Text):
    #    tobj.set_color('black')        
    fw, fh = fig.get_size_inches()
    bbox = fig.axes[2].get_position()
    aw, ah = bbox.width * fw, bbox.height * fh
    xlim = fig.axes[2].get_xlim()
    ylim = fig.axes[2].get_ylim()
    for p in fig.axes[2].patches:
        bbox1 = p.get_bbox()
        pw = bbox1.width/(xlim[1]-xlim[0]) * aw * 72
        ph = bbox1.height/(ylim[1]-ylim[0]) * ah * 72
        sw = ph/72/aw*(xlim[1]-xlim[0])/wfactor
        p.set_width(sw)
    return
bone.adjustFigure = adjustFigure

def getMSigDB(gs):
    url = "https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=" + gs + "&fileType=txt"
    df = pd.read_csv(url, sep="\t")
    df.columns.values[0] = 'ID'
    l1 = [list(df.ID[1:])]
    wt1 = [1]
    return wt1, l1
bone.getMSigDB = getMSigDB

def getYang2021(self, tn=1):
    self.prepareData("MACV412")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'LPS']
    ahash = {'saline':0, 'challenge':1}
    if tn == 2:
        bhash = {'GSM4291101', 'GSM4291098', 'GSM4291105', 'GSM4291100'}
        atype = [atype[i] if self.h.headers[i] not in bhash
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2021 = getYang2021

def getMarwick2018(self, tn=1):
    self.prepareData("MACV413")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub(" repl.*", "", str(k)) for k in atype]
    atype = [hu.re.sub(" pati.*", "", str(k)) for k in atype]
    atypes = ['C', 'LPS', 'AM', 'AM-IPF', 'AN', 'AN-LPS']
    ahash = {'Control MDM':0, 'MDM treated with LPS':1,
            'AM from RB-ILD':2, 'AM from IPF':3, 'MDM co-cultured with AN':4,
            'MDM co-cultured with AN and treated with LPs':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMarwick2018 = getMarwick2018

def getPatel2017Mac(self, tn=1):
    self.prepareData("MACV414")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub(" Don.*", "", str(k)) for k in atype]
    atype = [hu.re.sub("_Lung$", "", str(k)) for k in atype]
    atype = [hu.re.sub("_Blood$", "", str(k)) for k in atype]
    atypes = ['iMDC', 'iMDM', 'MDC', 'MDM', 'AM', 'Lung']
    ahash = {'Blood_Immature Monocyte derived DC':0,
            'Blood_Immature Monocyte derived macrophages':1,
            'Blood_Mature Monocyte derived DC':2,
            'Blood_Mature Monocyte derived macrophages':3,
            'Lung_Alveolar macrophages':4,
            'Lung_BDCA1+CD14+':5, 'Lung_BDCA1+CD14- DC':5,
            'Lung_BDCA1+CD14+ DC':5, 'Lung_CD14- DC':5, 'Lung_CD14+ DC':5,
            'Lung_Langerin+ DC':5}
    if tn == 2:
        atypes = ['MDM', 'AM']
        ahash = {'Blood_Mature Monocyte derived macrophages':0,
                'Lung_Alveolar macrophages':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPatel2017Mac = getPatel2017Mac

def getSMgutAD(self, tn=1):
    self.prepareData("AD55")
    age = self.h.getSurvName("c Age")
    atype = self.h.getSurvName("c Genotype")
    atypes = ['WT', 'PS19']
    ahash = {}
    if (tn == 2):
        atype = [str(atype[i])+"-"+str(age[i])
                for i in range(len(atype))]
        atypes = ['WT-3', 'WT-9', 'PS19-3', 'PS19-9']
        ahash = {'PS19-3 months':2, 'PS19-9 months':3,
                'WT-3 months':0, 'WT-9 months':1}
    if (tn == 3):
        atype = [str(atype[i])+"-"+str(age[i])
                for i in range(len(atype))]
        atypes = ['WT-3', 'PS19-3']
        ahash = {'PS19-3 months':1, 'WT-3 months':0}
    if (tn == 4):
        atype = [str(atype[i])+"-"+str(age[i])
                for i in range(len(atype))]
        atypes = ['WT-9', 'PS19-9']
        ahash = {'PS19-9 months':1, 'WT-9 months':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSMgutAD = getSMgutAD
def getPeters(self, tn=1):
    self.prepareData("PLP7")
    atype = self.h.getSurvName("c clinical condition")
    atypes = ['Normal', 'UC', 'CD']
    ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c gender")
        atypes = ['F', 'M']
        ahash = {'female':0, 'male':1}
        atype = [atype[i] if aval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atypes = ['Normal', 'UC']
        ahash = {"control": 0, "Ulcerative Colitis":1}
    if (tn == 4):
        atypes = ['Normal', 'CD']
        ahash = {"control": 0, "Crohn's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeters = getPeters


def getNoble(self):
    self.prepareData("PLP6")
    atype = self.h.getSurvName("c disease")
    atypes = ['Normal', 'UC']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNoble = getNoble


def getArijs2018(self, tn=1):
    self.prepareData("PLP10")
    atype = self.h.getSurvName("c Response")
    atypes = ['Control', 'UC R', 'UC NR', 'Active UC', 'UC other']
    ahash = {}
    if (tn == 2):
        atypes = ['UC R', 'UC NR']
    if (tn == 3):
        pid = self.h.getSurvName("c study individual number")
        res = self.h.getSurvName("c Response")
        phash = {}
        for i in range(len(atype)):
            if res[i] == 'UC R':
                phash[pid[i]] = 'R'
            if res[i] == 'UC NR':
                phash[pid[i]] = 'NR'
        time = self.h.getSurvName("c week (w)")
        atype = [phash[pid[i]] if pid[i] in phash and time[i] == 'W0'
                else None for i in range(len(atype))]
        atypes = ['R', 'NR']
        ahash = {}
    if (tn == 4):
        pid = self.h.getSurvName("c study individual number")
        res = self.h.getSurvName("c Response")
        phash = {}
        for i in range(len(atype)):
            if res[i] == 'UC R':
                phash[pid[i]] = 'R'
            if res[i] == 'UC NR':
                phash[pid[i]] = 'NR'
        time = self.h.getSurvName("c week (w)")
        therapy = self.h.getSurvName("c induction therapy_maintenance therapy")
        for i in range(len(atype)):
            if pid[i] in phash and time[i] == 'W0':
                atype[i] = "Pre " + phash[pid[i]]
            if pid[i] in phash and time[i] != 'W0':
                atype[i] = "Post " + phash[pid[i]]
        atype = [atype[i] if therapy[i] == 'IFX'
                else None for i in range(len(atype))]
        atypes = ['Pre R', 'Pre NR', 'Post R', 'Post NR']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArijs2018 = getArijs2018


def getWu2007(self):
    self.prepareData("PLP12")
    atype = self.h.getSurvName("c type")
    atypes = ['N', 'UC', 'CD']
    ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 2}
    atypes = ['N', 'UC-un', 'CD-un', 'IC-un', 'INF', 
            'UC-Aff', 'CD-Aff', 'IC-Aff']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWu2007 = getWu2007


def getVancamelbeke(self):
    self.prepareData("PLP16")
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'UC', 'CD']
    ahash = {'Biopsy from inflamed colonic mucosa of active UC patient':1,
            'Biopsy from inflamed colonic mucosa of active CD patient':2,
            'Biopsy from normal colonic mucosa of control individual':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVancamelbeke = getVancamelbeke


def getDePreter(self, t1 = 1):
    self.prepareData("PLP24")
    atype = self.h.getSurvName("c src1")
    atypes = ['UC Rp', 'UC Rb', 'UC p', 'UC b']
    ahash = {'Colonic mucosal biopsy from UC patient in remission before probiotics intake':1,
            'Colonic mucosal biopsy from UC patient in remission before placebo intake':0,
            'Colonic mucosal biopsy from active UC patient before placebo intake':2,
            'Colonic mucosal biopsy from active UC patient before probiotics intake':3}
    if t1 == 2:
        atypes = ['UC R', 'UC']
        ahash = {'Colonic mucosal biopsy from UC patient in remission before placebo intake':0,
                'Colonic mucosal biopsy from active UC patient before placebo intake':1}
        if t1 == 3:
            atypes = ['UC R', 'UC']
        ahash = {'Colonic mucosal biopsy from UC patient in remission before probiotics intake':0,
                'Colonic mucosal biopsy from active UC patient before probiotics intake':1}
        if t1 == 4:
            atypes = ['UC R', 'UC']
        ahash = {'Colonic mucosal biopsy from UC patient in remission before probiotics intake':0,
                'Colonic mucosal biopsy from UC patient in remission before placebo intake':0,
                'Colonic mucosal biopsy from active UC patient before placebo intake':1,
                'Colonic mucosal biopsy from active UC patient before probiotics intake':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDePreter = getDePreter


def getArijs2009(self, tn=1):
    self.prepareData("PLP27")
    atype = self.h.getSurvName("c tissue")
    ahash = {'Colon':0, 'Ileum':1}
    tissue = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c before or after first infliximab treatment")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'Before':1, 'After':2, 'Not':0}
    treatment = [ahash[i] if i in ahash else None for i in atype]
    response = self.h.getSurvName("c response to infliximab")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'UC', 'CD']
    ahash = {}
    if (tn == 2):
        atypes = ["R", "NR"]
        ahash = {"Yes": 0, "No": 1}
        atype = response
        atype = [atype[i] if tissue[i] == 0 and treatment[i] == 1
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ["R", "NR"]
        ahash = {"Yes": 0, "No": 1}
        atype = response
        atype = [atype[i] if tissue[i] == 0 and treatment[i] == 2
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArijs2009 = getArijs2009


def getHaberman2014(self, tn=1):
    self.prepareData("PLP11")
    atype = self.h.getSurvName("c deep ulcer")
    ahash = {'NA':0, 'No':1, 'Yes':2, 'no':1}
    ulcer = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['Control', 'UC', 'CD']
    ahash = {'UC':1, 'Not IBD':0, 'CD':2, 'not IBD':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if ulcer[i] == 0 or ulcer[i] == 1
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [ulcer[i] if aval[i] == 2
                else None for i in range(len(atype))]
        atypes = ['No', 'Yes']
        ahash = {1:0, 2:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHaberman2014 = getHaberman2014


def getHaberman2018(self):
    self.prepareData("PLP14")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['Control', 'CD']
    ahash = {'Non-IBD':0, 'CD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHaberman2018 = getHaberman2018


def getVanhove(self, dtype=0, tn=1):
    self.prepareData("PLP23")
    activity = self.h.getSurvName("c disease activity")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'UC', 'CD', 'I', 'A', 'N']
    ahash = {'ulcerative colitis':1, "Crohn's disease":2, 'control':0,
            'active':4, 'inactive':3, 'normal': 5}
    if (tn == 2):
        atypes = ['I', 'A']
        ahash = {'active':1, 'inactive':0}
        atype = activity
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanhove = getVanhove


def getVanderGoten(self, tn=1):
    self.prepareData("PLP25")
    activity = self.h.getSurvName("c disease activity")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'UC', 'CD', 'I', 'A', 'NA']
    ahash = {'control':0,
            'active':4, 'inactive':3, 'not applicable': 5}
    if (tn == 2):
        atypes = ['I', 'A']
        ahash = {'active':1, 'inactive':0}
        atype = activity
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanderGoten = getVanderGoten


def getPekow(self):
    self.prepareData("PLP59")
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'qUC', 'nUC']
    ahash = {'normal control': 0, 
            'quiescent ulcerative colitis': 1,
            'ulcerative colitis with neoplasia': 2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPekow = getPekow


def getGao(self, tn=1):
    self.prepareData("PLP34")
    atype = self.h.getSurvName("c group")
    atypes = ['C', 'D', 'A', 'A/D']
    ahash = {'control': 0, 'DSS': 1, 'AOM': 2, 'AOM/DSS': 3}
    if (tn == 2):
        atypes = ['C', 'D']
        ahash = {'control': 0, 'DSS': 1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGao = getGao


def getWatanabe(self):
    self.prepareData("PLP57")
    atype = self.h.getSurvName("c desc")
    atypes = ['UC-NonCa', 'UC-Ca']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWatanabe = getWatanabe


def getEColi(self):
    self.prepareData("CRC141")
    series = self.h.getSurvName("c Series")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'K12', 'O157']
    ahash = {'control, 60min':0, 'K-12, 60min':1,
            'O157:H7, 120min':2, 'control, 90min':0,
            'O157:H7, 60min':2, 'O157:H7, 90min':2,
            'control, 120min':0, 'K12, 120min':1, 'K-12, 90min':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEColi = getEColi


def getMatsuki(self):
    self.prepareData("CRC142")
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'Lc', 'Bb']
    ahash = {'Caco-2 cells cultured with B. breve':2,
            'Caco-2 cells alone':0,
            'Caco-2 cells cultured with L. casei': 1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatsuki = getMatsuki


def getArbibe(self):
    self.prepareData("CRC143")
    atype = self.h.getSurvName("c Title")
    atypes = ['C', 'M90T', 'OspF', 'OspF C']
    ahash = {'Caco2_OspF complementation infected_Rep1':3,
            'Caco2_Non infected_Rep2':0,
            'Caco2_OspF mutant infected_Rep1':2,
            'Caco2_OspF complementation infected_Rep2':3,
            'Caco2_M90T wild type infected_Rep1':1,
            'Caco2_Non infected_Rep1':0,
            'Caco2_M90T wild type infected_Rep2':1,
            'Caco2_OspF mutant infected_Rep2':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArbibe = getArbibe


def getPereiraCaro(self):
    self.prepareData("CRC141")
    series = self.h.getSurvName("c Series")
    atype = self.h.getSurvName("c agent")
    atypes = ['C', 'HTy', 'EHTy']
    ahash = {'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPereiraCaro = getPereiraCaro


def getKonnikova(self):
    self.prepareData("PLP68")
    media = self.h.getSurvName("c collection media")
    atype = self.h.getSurvName("c status")
    atypes = ['uninflamed', 'inflamed', 'D', 'R']
    ahash = {'RNA Later':3, 'DMSO':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKonnikova = getKonnikova


def getKarns2019(self):
    self.prepareData("PLP69")
    atype = self.h.getSurvName("c disease subtype")
    atypes = ['Control', 'UC', 'CD', 'iCD']
    ahash = {'ileal CD':3, 'not IBD':0, 'colonic CD':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKarns2019 = getKarns2019


def getPeck2015(self):
    self.prepareData("PLP70")
    inflamed = self.h.getSurvName("c inflamed")
    etype = self.h.getSurvName("c Type")
    atype = self.h.getSurvName("c disease_stage")
    atypes = ['NA', 'B1', 'B2', 'B3']
    ahash = {'B1/non-strictuing, non-penetrating':1,
            'B3/penetrating':3,
            'NA':0,
            'B2/stricturing':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeck2015 = getPeck2015


def getCorraliza2018(self):
    self.prepareData("PLP71")
    response = self.h.getSurvName("c hsct responder")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'CD', "YES", "NO", "C"]
    ahash = {'Healthy non-IBD':0, "Crohn's disease (CD)":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCorraliza2018 = getCorraliza2018


def getArze2019(self):
    self.prepareData("PLP72")
    atype = self.h.getSurvName("c disease status")
    atypes = ['Control', 'UC', "CD"]
    ahash = {'Non IBD':0, 'Ulcerative Colitis':1, "Crohn's Disease":2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArze2019 = getArze2019


def getVerstockt2019(self):
    self.prepareData("PLP73")
    atype = self.h.getSurvName("c clinical history")
    atypes = ['R', 'NR']
    ahash = {'responder':0, 'non-responder':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVerstockt2019 = getVerstockt2019


def getHasler(self):
    self.prepareData("PLP74")
    tissue = self.h.getSurvName("c tissue")
    inflammation = self.h.getSurvName("c inflammation")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['Control', 'UC', "CD", 'non inflamed', 'inflamed']
    ahash = {'disease control':0, 'healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHasler = getHasler


def getZhao2019(self):
    self.prepareData("PLP75")
    atype = self.h.getSurvName("c disease state")
    atype = self.h.getSurvName("c src1")
    atypes = ['Normal', "CD u", "CD i"]
    ahash = {'control':0,
            'Crohn\xe2\x80\x99s disease uninvolved':1,
            'Crohn\xe2\x80\x99s disease involved':2}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhao2019 = getZhao2019


def getKugathasan2008(self):
    self.prepareData("PLP76")
    atype = self.h.getSurvName("c Type")
    atypes = ['Normal', 'UC', "CD", "CD i"]
    ahash = {'Healthy control':0,
            'Colon-only CD':2,
            'Ileo-colonic CD':3,
            'Ulcerative colitis':1,
            'Internal control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKugathasan2008 = getKugathasan2008


def getZhao2015(self):
    self.prepareData("PLP77")
    atype = self.h.getSurvName("c disease state")
    atypes = ['Normal', 'UC I', "UC A"]
    ahash = {'ulcerative colitis inactive':1,
            'healthy control':0,
            'ulcerative colitis active':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhao2015 = getZhao2015


def getTang2017(self):
    self.prepareData("PLP78")
    state = self.h.getSurvName("c inflammation")
    atype = self.h.getSurvName("c disease state")
    atypes = ['Normal', 'UC', "CD", 'Inactive', "Active"]
    ahash = {'non-IBD control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTang2017 = getTang2017


def getCarey2008(self):
    self.prepareData("PLP79")
    atype = self.h.getSurvName("c Type")
    atypes = ['Normal', 'UC', "CD", "CD t"]
    ahash = {'healthy control reference':0,
            'CD':2,
            'treated CD':3,
            'UC':1,
            'Internal Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCarey2008 = getCarey2008


def getDotti2017(self):
    self.prepareData("PLP80")
    culture = self.h.getSurvName("c organoid_culture")
    atype = self.h.getSurvName("c case_phenotype")
    atypes = ['Normal', 'UC'];
    ahash = {'ulcerative colitis (UC) patient':1, 'non-IBD control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDotti2017 = getDotti2017


def getDenson2018(self):
    self.prepareData("PLP86")
    response = self.h.getSurvName("c week 4 remission")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['Control', 'UC', 'Yes', 'No', 'NA'];
    ahash = {'Ulcerative Colitis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDenson2018 = getDenson2018


def getBoyd2018(self):
    self.prepareData("PLP87")
    response = self.h.getSurvName("c condition")
    ahash = {'CD active':4, 'CD inactive':3, 'control':5,
            'UC active':4, 'UC inactive':3}
    rval = [ahash[i] if i in ahash else None for i in response]
    atype = self.h.getSurvName("c condition")
    atypes = ['Control', 'UC', 'CD', 'Inactive', 'Active', 'NA'];
    ahash = {'CD active':2, 'CD inactive':2, 'control':0,
            'UC active':1, 'UC inactive':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBoyd2018 = getBoyd2018


def getBreynaert2013(self, tn=1):
    self.prepareData("PLP38")
    atype = self.h.getSurvName("c colitis group")
    atypes = ['C', 'DA', 'D1', 'D2', 'D3', 'A'];
    ahash = {'2 cycles DSS with additional recovery period':1,
            'control':0,
            '1 cycle DSS':2,
            'acute colitis':5,
            '3 cycles DSS':4,
            '2 cycles DSS':3}
    if (tn == 2):
        atypes = ['C', 'D']
        ahash = {'control':0, '3 cycles DSS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBreynaert2013 = getBreynaert2013


def getGerstgrasser2017(self):
    self.prepareData("PLP36")
    atype = self.h.getSurvName("c Title")
    atypes = ['C', 'DSS'];
    ahash = {'colon wt DSS rep3':1,
            'colon wt DSS rep2':1,
            'colon wt DSS rep4':1,
            'colon wt DSS rep1':1,
            'colon wt water rep3':0,
            'colon wt water rep2':0,
            'colon wt water rep1':0,
            'colon wt water rep4':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGerstgrasser2017 = getGerstgrasser2017


def getTang2012(self, tn=1):
    self.prepareData("PLP40")
    atype = self.h.getSurvName("c disease state")
    atypes = ['N', 'I', 'LD', 'HD', 'C'];
    ahash = {'low grade dysplasia lesion':2,
            'inflamed colorectal mucosa':1,
            'high grade dysplasia':3,
            'normal':0,
            'colorectal adenocarcinoma':4}
    if (tn == 2):
        atypes = ['N', 'D']
        ahash = {'low grade dysplasia lesion':1,
                'inflamed colorectal mucosa':1,
                'high grade dysplasia':1,
                'normal':0,
                'colorectal adenocarcinoma':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTang2012 = getTang2012


def getJensen2017(self):
    self.prepareData("PLP66")
    atype = self.h.getSurvName("c disease")
    atypes = ['normal', 'colitis']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJensen2017 = getJensen2017


def getGkouskou2016(self, tn=1):
    self.prepareData("PLP84")
    tissue = self.h.getSurvName("c tissue")
    ahash = {'proximal colon':3, 'distal colon':4}
    rval = [ahash[i] if i in ahash else None for i in tissue]
    atype = self.h.getSurvName("c src1")
    atypes = ['normal', 'AD2', 'AD4', 'proximal', 'distal']
    ahash = {'UNTREATED':0, 'AOM, 4 DSS CYCLES':2, 'AOM, 2 DSS CYCLES':1}
    if (tn == 2):
        atype = [str(atype[i])+ " " + str(tissue[i]) for i in
                range(len(atype))]
        atypes = ['proximal', 'distal']
        ahash = {'UNTREATED proximal colon':0,
                'UNTREATED distal colon':1}
    if (tn == 3):
        atypes = ['normal', 'colitis']
        ahash = {'UNTREATED':0, 'AOM, 4 DSS CYCLES':1, 'AOM, 2 DSS CYCLES':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGkouskou2016 = getGkouskou2016


def getGkouskou2016ProxDis(self):
    self.prepareData("PLP85")
    atype = self.h.getSurvName("c tissue")
    atypes = ['proximal', 'distal']
    ahash = {'PROXIMAL COLON':0, 'DISTAL COLON':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGkouskou2016ProxDis = getGkouskou2016ProxDis


def getLopezDee2012(self, tn=1):
    self.prepareData("PLP49")
    gtype = self.h.getSurvName("c genotype/variation")
    ahash = {'Wild-type':8, 'TSP-null':9}
    rval = [ahash[i] if i in ahash else None for i in gtype]
    atype = self.h.getSurvName("c src1")
    atypes = ['WC', 'TC', 'WD', 'TD', 'WDS', 'WDST2', 'WDST3', 'WRFK',
            'WT', 'TN']
    ahash = {'Wt, water control':0,
            'TSP-null-water':1,
            'Wt, DSS treated':2,
            'DSS-treated TSP-null':3,
            'Wt, DSS-saline treated':4,
            'Wt, DSS-TSR2 treated':5,
            'Wt, DSS-3TSR treated':6,
            'Wt, DSS-TSR+RFK treated':7}
    if (tn == 2):
        atypes = ['C', 'D']
        ahash = {'Wt, water control':0,
                'Wt, DSS treated':1,
                'Wt, DSS-saline treated':1,
                'Wt, DSS-TSR2 treated':1,
                'Wt, DSS-3TSR treated':1,
                'Wt, DSS-TSR+RFK treated':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLopezDee2012 = getLopezDee2012


def getSarvestani2018(self):
    self.prepareData("ORG16")
    atype = self.h.getSurvName("c Pathology")
    atypes = ['N', 'UC', 'D', 'C']
    ahash = {'Normal':0,
            'Chronic active colitis':1,
            'Colitis':1,
            'Colitis with benign strictures':1,
            'Colitis with fibrosis and treatment effect':1,
            'Low-grade dysplasia':2,
            'T3N1a':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSarvestani2018 = getSarvestani2018


def getPlanell2013(self):
    self.prepareData("PLP88")
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'NI', 'Re', 'I']
    ahash = {
            'Human colon biopsies from UC patient with active disease (involved mucosa)':3,
            'Human colon biopsies from non-inflammatory control':0,
            'Human colon biopsies from UC patient with active disease (non-involved mucosa)':1,
            'Human colon biopsies from UC patient in remission (involved mucosa)':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPlanell2013 = getPlanell2013


def getLyons2018(self):
    self.prepareData("PLP89")
    atype = self.h.getSurvName("c inflammation level")
    atypes = ['NI', 'M', 'S']
    ahash = {'severe':2, 'moderate':1, 'non-inflamed':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLyons2018 = getLyons2018


def getFang2012(self):
    self.prepareData("PLP90")
    atype = self.h.getSurvName("c time")
    atypes = ['W0', 'W2', 'W4', 'W6']
    ahash = {'0 week':0, '4 weeks':2, '6 weeks':3, '2 weeks':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFang2012 = getFang2012


def getSchiering2014(self):
    self.prepareData("PLP91")
    gtype = self.h.getSurvName("c genotype")
    ahash = {'wild type':3, 'Il23r-/-':4, 'Foxp3gfp':5}
    rval = [ahash[i] if i in ahash else None for i in gtype]
    atype = self.h.getSurvName("c cell type")
    atypes = ['TC', 'TP', 'TN', 'WT', 'I23', 'F3']
    ahash = {
            'TCR\xce\xb2+CD4+ T cells from colon':0,
            'TCR\xce\xb2+CD4+Foxp3+ from colon lamina propria (cLP)':1,
            'TCR\xce\xb2+CD4+Foxp3+ from mesenteric lymph node (MLN)':2}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSchiering2014 = getSchiering2014


def getKremer2012(self):
    self.prepareData("PLP92")
    atype = self.h.getSurvName("c disease status")
    atypes = ['N', 'UC']
    ahash = {'TNBS colitis':1, 'Healthy Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKremer2012 = getKremer2012


def getHo2014(self):
    self.prepareData("PLP93")
    gtype = self.h.getSurvName("c tissue")
    ahash = {'Spleen':3, 'Colon':4}
    rval = [ahash[i] if i in ahash else None for i in gtype]
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'UC', 'UCt', 'Sp', 'CO']
    ahash = {
            'Mock':0,
            'EA treatment/TNBS-induced colitis':2,
            'TNBS-induced colitis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHo2014 = getHo2014


def getDohi2014(self):
    self.prepareData("PLP94")
    gtype = self.h.getSurvName("c treated with")
    ahash = {'none (untreated control)':2,
            '10 mg/kg control IgG2a mAb (anti-human CD20)':3,
            '0.3 mg/kg TNFR-Ig':4,
            '10 mg/kg anti-TWEAK mP2D10':5,
            'combination of TNFR-Fc (0.3 mg/kg) and anti-TWEAK mP2D10 (10 mg/kg)':6}
    rval = [ahash[i] if i in ahash else None for i in gtype]
    atype = self.h.getSurvName("c injected with")
    atypes = ['N', 'UC']
    ahash = {'trinitrobenzene sulfonic acid (TNBS)':1,
            'none (na\xc3\xafve control)':0}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDohi2014 = getDohi2014


def getDeBuhr2006(self):
    self.prepareData("PLP95")
    atype = self.h.getSurvName("c Title")
    atypes = ['B6-WT', 'B6-IL10', 'C3-WT', 'C3-IL10']
    ahash = {'C57BL/6J, sample 1':0,
            'C57BL/6J, 4 week old, sample 2':0,
            'C57BL/6J-Il10tm1Cgn, sample 2':1,
            'C57BL/6J-Il10tm1Cgn, sample 1':1,
            'C3H/HeJBir, sample 2':2,
            'C3H/HeJBir, sample-1':2,
            'C3H/HeJBir-Il10tm1Cgn, sample 2':3,
            'C3H/HeJBir-Il10tm1Cgn, sample 1':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeBuhr2006 = getDeBuhr2006


def getRuss2013(self, tval):
    self.prepareData("PLP96")
    tissue = self.h.getSurvName("c tissue")
    ahash = {'colon epithelium':2, 'colon':3}
    rval = [ahash[i] if i in ahash else None for i in tissue]
    atype = self.h.getSurvName("c genotype/variation")
    atypes = ['WT', 'IL10']
    ahash = {'IL10-/-':1, 'wildtype':0}
    atype = [atype[i] if rval[i] == tval else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRuss2013 = getRuss2013


def getPunit2015(self):
    self.prepareData("PLP97")
    atype = self.h.getSurvName("c genotype")
    atypes = ['WT', 'TNFR2']
    ahash = {'Wildtype':0, 'TNFR2-knockout':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPunit2015 = getPunit2015


def getTam2019(self):
    self.prepareData("PLP98")
    atype = self.h.getSurvName("c genotype")
    atypes = ['WT', 'Tnfr1']
    ahash = {'Tnfrsf1a-/-':1, 'WT':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTam2019 = getTam2019


def getLamas2018(self, tn=1):
    self.prepareData("PLP99")
    atype = self.h.getSurvName("c src1")
    atypes = ['D0', 'D4', 'D12', 'D22']
    ahash = {'Cecum_Flore WT_Day0':0,
            'Cecum_Flore KO_Day0':0,
            'Cecum_Flore WT_Day4':1,
            'Cecum_Flore KO_Day4':1,
            'Cecum_Flore WT_Day12':2,
            'Cecum_Flore KO_Day12':2,
            'Cecum_Flore WT_Day22':3,
            'Cecum_Flore KO_Day22':3}
    if (tn == 2):
        atypes = ['WT', 'KO']
        ahash = {'Cecum_Flore WT_Day0':0,
                'Cecum_Flore KO_Day0':1,
                'Cecum_Flore WT_Day4':0,
                'Cecum_Flore KO_Day4':1,
                'Cecum_Flore WT_Day12':0,
                'Cecum_Flore KO_Day12':1,
                'Cecum_Flore WT_Day22':0,
                'Cecum_Flore KO_Day22':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLamas2018 = getLamas2018


def getFuso1(self):
    self.prepareData("CRC112.2")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'FN']
    ahash = {"non-infected": 0, "F. nucleatum":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFuso1 = getFuso1


def getGEOMac(self):
    self.prepareData("GL4")
    atype = self.h.getSurvName("c Source")
    atypes = ['Control', 'Treated']
    aval = [None, None]
    for i in range(2, len(atype)):
        if re.search(r'control', atype[i]) or \
                re.search(r'untreated', atype[i]):
                    aval += [0]
        else:
            aval += [1]
    self.st1 = [ i for i in self.h.aRange() if aval[i] == 0]
    self.st2 = [ i for i in self.h.aRange() if aval[i] == 1]
    self.st3 = []
    self.aval = aval
    self.atype = atype
    self.atypes = atypes
    self.order = self.st1 + self.st2 + self.st3
    self.printInfo()

def getGEOMacAnn(self):
    self.prepareData("G16")
    atype = self.h.getSurvName("c Type")
    atypes = ['M0', 'M1', "M2"]
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGEOMacAnn = getGEOMacAnn


def getBeyer2012(self):
    self.prepareData("MAC1")
    atype = self.h.getSurvName("c cell type")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'M1 macrophages':1, 'M0 macrophages':0, 'M2 macrophages':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBeyer2012 = getBeyer2012


def getSaliba2016(self):
    self.prepareData("MAC5")
    atype = self.h.getSurvName("c group")
    atypes = ['NM', 'M1', 'M2', 'BM']
    ahash = {'Macrophage with non growing bacteria':1,
            'Macrophage with growing bacteria':2,
            'Naive Macrophage':0,
            'Bystanders':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSaliba2016 = getSaliba2016


def getAhrens2013(self):
    self.prepareData("LIV9")
    atype = self.h.getSurvName("c group")
    atypes = ['C', 'HO', 'S', 'N']
    ahash = {'Control':0, 'Healthy obese':1, 'Nash':3, 'Steatosis':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAhrens2013 = getAhrens2013


def getSmith2015(self):
    self.prepareData("PLP82")
    atype = self.h.getSurvName("c disease")
    atypes = ['HC', 'CD']
    ahash = {'Healthy Control':0, "Crohn's Disease":1 }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2015 = getSmith2015


def getDill2012HPC(self):
    self.prepareData("LIV10")
    atype = self.h.getSurvName("c disease state")
    atypes = ['HC', 'AH']
    ahash = {'acute hepatitis':1, 'healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDill2012HPC = getDill2012HPC


def getDill2012IFN(self):
    self.prepareData("LIV11")
    atype = self.h.getSurvName("c treatment")
    atypes = ['Un', 'IFNa', 'TFNg']
    ahash = {'IFNalpha':1, 'IFNgamma':2, 'untreated':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDill2012IFN = getDill2012IFN


def getTrepo2018(self, tn=1, tb=0):
    self.prepareData("LIV12")
    atype = self.h.getSurvName("c cohort")
    ahash = {'Derivation cohort (BRAH)':0, '':2, 'Validation cohort (BRAH)':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'AH', 'AFL', 'AC']
    ahash = {'liver tissue':0,
            'Alcoholic_cirrhosis, liver tissue':3,
            'Mild_acute_alcoholic_hepatitis, liver tissue':1,
            'alcoholic_steatosis, liver tissue':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName('c outcome at 6 months')
        atypes = ['A', 'D']
        ahash = {'Dead or liver transplantation':1, 'Alive':0}
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTrepo2018 = getTrepo2018


def getTrepo2018II(self, tn=1, tb=0):
    self.prepareData("LIV13")
    atype = self.h.getSurvName("c cohort")
    ahash = {'BRAH cohort':0, 'MSMC cohort':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c outcome at 6 months')
    atypes = ['A', 'D']
    ahash = {'Dead or liver transplantation':1, 'Alive':0}
    if (tn == 2):
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTrepo2018II = getTrepo2018II


def getSuppli2019(self, tn=1):
    self.prepareData("LIV14")
    atype = self.h.getSurvName("c disease")
    atypes = ['H', 'O', 'FL', 'SH']
    ahash = {'NAFLD':2, 'NASH':3, 'obese':1, 'healthy':0}
    if (tn == 2):
        atypes = ['H', 'FL']
        ahash = {'NAFLD':1, 'healthy':0}
    if (tn == 3):
        atypes = ['H', 'SH']
        ahash = {'NASH':1, 'healthy':0}
    if (tn == 4):
        atypes = ['H', 'O']
        ahash = {'obese':1, 'healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSuppli2019 = getSuppli2019


def getHoshida2013(self, tn=1):
    self.prepareData("LIV15.2")
    time1 = self.h.getSurvName("c days to death")
    time1 = ["", ""] + [float(time1[k]) if time1[k] != 'NA' else None for k in self.h.aRange()]
    time2 = self.h.getSurvName("c days to decomp")
    time2 = ["", ""] + [float(time2[k]) if time2[k] != 'NA' else None for k in self.h.aRange()]
    time3 = self.h.getSurvName("c days to child")
    time3 = ["", ""] + [float(time3[k]) if time3[k] != 'NA' else None for k in self.h.aRange()]
    time4 = self.h.getSurvName("c days to hcc")
    time4 = ["", ""] + [float(time4[k]) if time4[k] != 'NA' else None for k in self.h.aRange()]
    status1 = self.h.getSurvName("c death")
    status1 = ["", ""] + [int(status1[k]) if status1[k] != 'NA' else None for k in self.h.aRange()]
    status2 = self.h.getSurvName("c decomp")
    status2 = ["", ""] + [int(status2[k]) if status2[k] != 'NA' else None for k in self.h.aRange()]
    status3 = self.h.getSurvName("c child")
    status3 = ["", ""] + [int(status3[k]) if status3[k] != 'NA' else None for k in self.h.aRange()]
    status4 = self.h.getSurvName("c hcc")
    status4 = ["", ""] + [int(status4[k]) if status4[k] != 'NA' else None for k in self.h.aRange()]
    days = 365 * 15
    atype = [None] * len(time1)
    for k in self.h.aRange():
        for g in [[time1, status1], [time2, status2],
                  [time3, status3], [time4, status4]]:
            if (g[0][k] is not None and g[1][k] is not None): 
                if (g[0][k] < days and g[1][k] == 1):
                    atype[k] = 1
                    break
                else:
                    atype[k] = 0
    atypes = [0, 1]
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHoshida2013 = getHoshida2013


def getRamilo2007v1(self):
    self.prepareData("BL9.1")
    atype = self.h.getSurvName("c Pathogen")
    atypes = ['None', 'Ecoli', 'MRSA', 'MSSA', 'InfA', 'pneu']
    ahash = {'S. aureus, MRSA':2,
            'E. coli':1,
            'S. aureus, MSSA':3,
            'None':0,
            'Influenza A':4,
            'S. pneumoniae':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRamilo2007v1 = getRamilo2007v1


def getPG2019lps(self, tn=1):
    self.prepareData("MAC28.3")
    atype = self.h.getSurvName("c Type")
    atypes = ['LPS0', 'LPS5', 'GIV0', 'GIV5', 'U']
    ahash = {'sh2_GIV_LPS_0hr':2, 'Undetermined':4, 'sh1_GIV_LPS_5hr':3,
            'shC_LPS_0hr':0, 'sh2_GIV_LPS_5hr':3, 'sh1_GIV_LPS_0hr':2,
            'shC_LPS_5hr':1}
    if (tn == 2):
        atypes = ['LPS0', 'LPS5', 'GIV0', 'GIV5']
        ahash = {'sh2_GIV_LPS_0hr':2, 'sh1_GIV_LPS_5hr':3, 'shC_LPS_0hr':0,
                'sh2_GIV_LPS_5hr':3, 'sh1_GIV_LPS_0hr':2, 'shC_LPS_5hr':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2019lps = getPG2019lps


def getPG2020lps(self, tn=1):
    self.prepareData('PG7', "/Users/dtv004/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c Protocol")
    atypes = ['C0', 'C5', 'C16', 'GIV0', 'GIV5', 'GIV16']
    ahash = {'shGIV 16h':5, 'shC 16h':2, 'shGIV 5h':4,
            'shC 5h':1, 'shGIV 0h':3, 'shC 0h':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020lps = getPG2020lps


def getHugo2016(self):
    self.prepareData("ML12")
    atype = self.h.getSurvName("c anti-pd-1 response")
    atypes = ['CR', 'PR', 'PD'];
    ahash = {'Progressive Disease':2,
            'Partial Response':1,
            'Complete Response':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHugo2016 = getHugo2016


def getPrat2017(self):
    self.prepareData("ML13")
    atype = self.h.getSurvName("c response")
    atypes = ['RC_RP_SD', 'PD'];
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPrat2017 = getPrat2017


def getCassetta2019(self):
    self.prepareData("MAC29")
    atype = self.h.getSurvName("c condition")
    atypes = ['N', 'NB', 'NE', 'BC', 'EC'];
    ahash = {'Endometrial cancer':4,
            'Normal':0,
            'Breast cancer':3,
            'Normal Breast':1,
            'Normal Endometrial':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCassetta2019 = getCassetta2019


def getPoczobutt2016(self):
    self.prepareData("MAC30")
    atype = self.h.getSurvName("c src1")
    atypes = ['L', 'TL2', 'TL3']
    ahash = {'tumor-bearing lung, 2 week':1,
            'tumor-bearing lung, 3 week':2,
            'lung':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPoczobutt2016 = getPoczobutt2016


def getWoetzel2014(self, tn=1):
    self.prepareData("MAC31")
    atype = self.h.getSurvName("c disease state")
    atype2 = self.h.getSurvName("c clinical status")
    atype = [atype[i] + atype2[i] for i in range(len(atype))]
    atypes = ['HC', 'RA', 'OA']
    ahash = {'healthy control':0,
            'rheumatoid arthritis':1,
            'synovial tissue isolated from osteoarthritic joint':2,
            'osteoarthritis':2,
            'normal control':0}
    if (tn == 2):
        atypes = ['HC', 'OA']
        ahash = {'healthy control':0,
                'synovial tissue isolated from osteoarthritic joint':1,
                'osteoarthritis':1, 'normal control':0}
    if (tn == 3):
        atypes = ['HC', 'RA']
        ahash = {'healthy control':0, 'normal control':0,
                'rheumatoid arthritis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWoetzel2014 = getWoetzel2014


def getXue2014(self, mt=1):
    self.prepareData("MAC2")
    atype = self.h.getSurvName("c Type")
    atypes = ['M0', 'M1', 'M2']
    ahash = {\
            'M_GMCSF_IL4_72h':2,
            'M_GMCSF_IFNg_72h':1,
            'M0_GMCSF_72h':0}
    if mt == 2:
        ahash = {\
                'M_GMCSF_IL4_24h':2,
                'M_GMCSF_IFNg_24h':1,
                'M0_GMCSF_24h':0}
    if mt == 3:
        ahash = {\
                'M_GMCSF_IL4_12h':2,
                'M_GMCSF_IFNg_12h':1,
                'M0_GMCSF_12h':0}
    if mt == 4:
        ahash = {\
                'M0_GMCSF_0h':0,
                'M0_GMCSF_12h':0,
                'M0_GMCSF_24h':0,
                'M0_GMCSF_48h':0,
                'M0_GMCSF_6h':0,
                'M0_GMCSF_72h':0,
                'M0_MCSF_0h':0,
                'M1/2_GMCSF_24h':0,
                'M_GMCSF_IFNg_30min':0,
                'M_GMCSF_IFNg_1h':0,
                'M_GMCSF_IFNg_2h':0,
                'M_GMCSF_IFNg_4h':1,
                'M_GMCSF_IFNg_6h':1,
                'M_GMCSF_IFNg_12h':1,
                'M_GMCSF_IFNg_24h':1,
                'M_GMCSF_IFNg_72h':1,
                'M_GMCSF_IL4_30min':2,
                'M_GMCSF_IL4_1h':2,
                'M_GMCSF_IL4_2h':2,
                'M_GMCSF_IL4_4h':2,
                'M_GMCSF_IL4_6h':2,
                'M_GMCSF_IL4_12h':2,
                'M_GMCSF_IL4_24h':2,
                'M_GMCSF_IL4_72h':2,
                'M_MCSF_IL4_72h':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXue2014 = getXue2014


def getShaykhiev2009(self, tn=1):
    self.prepareData("MAC11.2")
    atype = self.h.getSurvName("c desc")
    atype = [str(i).split("-")[0] for i in atype]
    atypes = ['NS', 'S', 'COPD']
    ahash = {}
    if (tn == 2):
        atypes = ['H', 'COPD']
        ahash = {'NS':0, 'S':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShaykhiev2009 = getShaykhiev2009


def getWoodruff2005(self, tn=1):
    self.prepareData("MAC12")
    atype = self.h.getSurvName("c status")
    atypes = ['NS', 'S', 'A']
    ahash = {'Asthmatic':2, 'Smoker':1, 'Nonsmoker':0}
    if (tn == 2):
        atypes = ['H', 'A']
        ahash = {'Asthmatic':1, 'Smoker':0, 'Nonsmoker':0}
    if (tn == 3):
        atypes = ['NS', 'S']
        ahash = {'Smoker':1, 'Nonsmoker':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWoodruff2005 = getWoodruff2005


def getWS2009(self, tn=1):
    self.prepareData("MAC12.2")
    atype = self.h.getSurvName("c Type")
    atypes = ['NS', 'S', 'A', 'C']
    ahash = {'Nonsmoker':0, 'Asthmatic':2, 'Smoker':1, 'smoker':1,
            'non-smoker':0, 'COPD':3}
    if (tn == 2):
        atypes = ['NS', 'A']
        ahash = {'Nonsmoker':0, 'Asthmatic':1, 'non-smoker':0}
    if (tn == 3):
        atypes = ['NS', 'COPD']
        ahash = {'Nonsmoker':0, 'non-smoker':0, 'COPD':1}
    if (tn == 4):
        atypes = ['NS', 'S']
        ahash = {'Nonsmoker':0, 'non-smoker':0, 'Smoker':1, 'smoker':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWS2009 = getWS2009


def getZhang2015(self, tn=1):
    self.prepareData("MAC3")
    atype = self.h.getSurvName("c Type")
    atypes = ['M', 'iM', 'M1', 'iM1', 'M2', 'iM2', 'i']
    ahash = {'IPSDM M2':5,
            'IPSDM MAC':1,
            'IPSDM M1':3,
            'iPS':6,
            'HMDM M1':2,
            'HMDM MAC':0,
            'HMDM M2':4}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'IPSDM M2':2,
                'IPSDM MAC':0,
                'IPSDM M1':1}
    if (tn == 3):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'HMDM M1':1,
                'HMDM MAC':0,
                'HMDM M2':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2015 = getZhang2015


def getHaribhai2016(self):
    self.prepareData("MAC13")
    atype = self.h.getSurvName("c cell type")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'M2a macrophages':2, 'M0 macrophages':0, 'M1 macrophages':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHaribhai2016 = getHaribhai2016


def getOhradanovaRepic2018(self, tn=1):
    self.prepareData("MAC14")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2', 'IL10']
    ahash = {'mock-activated (medium only; control) for 2d':0,
            'activated with 100 ng/ml LPS + 25ng/ml IFN\xce\xb3 for 2d':1,
            'activated with 20 ng/ml IL-4 for 2d':2,
            'activated with 20 ng/ml IL-10 for 2d':3}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'mock-activated (medium only; control) for 2d':0,
                'activated with 100 ng/ml LPS + 25ng/ml IFN\xce\xb3 for 2d':1,
                'activated with 20 ng/ml IL-4 for 2d':2}
    if (tn == 3):
        atypes = ['M0', 'IL10']
        ahash = {'mock-activated (medium only; control) for 2d':0,
                'activated with 20 ng/ml IL-10 for 2d':1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOhradanovaRepic2018 = getOhradanovaRepic2018


def getGharib2019CF(self, val = 0):
    self.prepareData("MAC15")
    atype = self.h.getSurvName("c patient identification number")
    atype = [str(i).split(" ")[0] for i in atype]
    ahash = {'Non':0, 'CF':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c title")
    atype = [ str(i).split(" ")[-1] for i in atype ]
    atypes = ['None', 'IL4', 'IL10', 'MP', 'Az', 'PMNs']
    ahash = {'alone':0,
            'methylprednisolone':3,
            'azithromycin':4}
    atype = [atype[i] if rval[i] == val else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGharib2019CF = getGharib2019CF


def getGharib2019Alv(self, tn=1, ri=0):
    self.prepareData("MAC16")
    atype = self.h.getSurvName("c time point")
    atypes = ['D1', 'D4', 'D8']
    ahash = {'Day 1':0, 'Day 4':1, 'Day 8':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("n ventilator-free days (vfd)")
    atypes = ['VFD-Extubated/Alive', 'VFD-Intubated/Dead']
    ahash = {'0':1, '7':0, '18':0, '19':0, '21':0,
            '22':0, '23':0, '24':0, '25':0}
    atype = [atype[i] if rval[i] == ri else None for i in range(len(atype))]
    if (tn == 2):
        atypes = ['VFD-Extubated/Alive', 'VFD-Intubated/Dead']
        ahash = {'0':1, '7':1, '18':1, '19':0, '21':0,
                '22':0, '23':0, '24':0, '25':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGharib2019Alv = getGharib2019Alv


def getOrozco2012(self):
    self.prepareData("MAC17")
    atype = self.h.getSurvName("c treatment condition")
    atypes = ['Control', 'LPS', 'OxPAPC']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOrozco2012 = getOrozco2012


def getRegan2018(self):
    self.prepareData("MAC18")
    treatment = self.h.getSurvName("c drug treatment")
    ahash = {'1 UNT':0,
            '2 CTRL':1,
            '3 Escitalopram':2,
            '4 Nortriptyline':3,
            '5 Anti-TNFa':4,
            '6 Indomethacin':5,
            '7 Prednisolone':6}
    rval = [ahash[i] if i in ahash else None for i in treatment]
    atype = self.h.getSurvName("c inflammatory stimulus")
    atypes = ['None', 'IFN7', 'IFN24', 'LPS7', 'LPS24']
    ahash = {'24 h No inflammation':0,
            'IFN 07 h':1,
            'IFN 24 h':2,
            'LPS 07 h':3,
            'LPS 24 h':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRegan2018 = getRegan2018


def getMartinez2013II(self, tn=1):
    self.prepareData("MAC19.5")
    atype = self.h.getSurvName("c src1")
    atypes = ['None', 'Mono', 'IFN', 'IL4', 'IL10', 'MCSF', 'Med']
    ahash = {'monocyte-derived macrophages, IFN-y':2,
            'monocyte-derived macrophages, M-CSF':5,
            'monocyte-derived macrophages':0,
            'monocyte-derived macrophages, IL-4':3,
            'monocyte-derived macrophages, IL-10':4,
            'monocytes':1,
            'monocyte-derived macrophages, Medium':6}
    if (tn == 2):
        atypes = ['Mono', 'Mac']
        ahash = {'monocyte-derived macrophages':1, 'monocytes':0}
    if (tn == 3):
        atypes = ['Mac', 'MCSF']
        ahash = {'monocyte-derived macrophages':0,
                'monocyte-derived macrophages, M-CSF':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinez2013II = getMartinez2013II


def getMartinez2013(self):
    self.prepareData("MAC19.1")
    atype = self.h.getSurvName("c Type")
    atypes = ['M0', 'M1', 'M2', 'T0', 'T3']
    ahash = {'Monocyte at 3 days':4,
            'classical or M1 activated macrophages':1,
            'Macrophage at 7 days':0,
            'Alternative or M2 activated macrophages':2,
            'Monocyte at T0':3,
            'Monocyteat 3 days':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinez2013 = getMartinez2013


def getMartinez2013Mm(self):
    self.prepareData("MAC19.3")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'none':0, '18 hours with 20 ng/ml of mIL-4':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinez2013Mm = getMartinez2013Mm


def getKoziel2009(self):
    self.prepareData("MAC22")
    atype = self.h.getSurvName("c characteristics")
    atypes = ['control hMDMs', 'SA-treated hMDMs']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKoziel2009 = getKoziel2009


def getJiang2017(self):
    self.prepareData("MAC24")
    atype = self.h.getSurvName("c source name")
    atypes = ['M1', 'M2']
    ahash = {'Bone marrow-derived macrophages, M1, LPS+IFN-r stimulated BMDM':0,
            'Bone marrow-derived macrophages, M2, IL-4 stimulated BMDM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJiang2017 = getJiang2017


def getHan2017(self, tn=1):
    self.prepareData("MAC25")
    atype = self.h.getSurvName("c Title")
    atype = [str(i).split(" ")[2] if len(str(i).split(" ")) > 3 else i \
            for i in atype]
    atypes = ['SR1078', 'M1', 'M0', 'Veh', 'M2', 'SR3335']
    ahash = {}
    if tn == 2:
        atypes = ['M0', 'M1', 'M2']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHan2017 = getHan2017


def getFuentesDuculan2010(self):
    self.prepareData("MAC26")
    atype = self.h.getSurvName("c source name")
    ahash = {'Macrophages culture':0,
            'Psoriasis Non-lesional skin':1,
            'Psoriasis Lesional skin':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment group")
    atypes = ['control', 'IL17', 'IFNg', 'IL4', 'TNFa', 'LPS', "M1", "PS" ]
    ahash = {'LPS and IFNg':6, "":7}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFuentesDuculan2010 = getFuentesDuculan2010


def getZhou2009(self):
    self.prepareData("MAC27")
    atype = self.h.getSurvName("c time")
    ahash = {'1h':0, '4h':1, '24h':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c disease state")
    atypes = ['control (lean)', 'obesity']
    ahash = {}
    atype = [atype[i] if rval[i] == 2 else None for i in range(len(atype))]
    self.rval = rval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhou2009 = getZhou2009


def getSurvival(self, dbid = "CRC35.3"):
    self.prepareData(dbid)
    atype = self.h.getSurvName("status")
    atypes = ['Censor', 'Relapse']
    ahash = {"0": 0, "1":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSurvival = getSurvival


def getJablonski2015(self):
    self.prepareData("MAC32")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'received media alone (M0 condition)':0,
            'classically activated (M1 condition) with LPS + IFN-gamma':1,
            'alternatively activated (M2 condition) with IL-4':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJablonski2015 = getJablonski2015


def getZhao2017(self):
    self.prepareData("MAC33")
    atype = self.h.getSurvName("c genotype/variation")
    atype2 = self.h.getSurvName("c Stage")
    atype = [ str(atype[i]) + " " + str(atype2[i]) for i in range(self.h.end+1)]
    atypes = ['W5', 'K5', 'W24', 'K24']
    ahash = {'wildtype week 5':0,
            'Mecp2 knockout week 5':1,
            'Mecp2 knockout week 24':3,
            'wildtype week 24':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhao2017 = getZhao2017


def getChiu2013(self):
    self.prepareData("MAC34")
    atype = self.h.getSurvName("c sample type")
    atypes = ['U', 'L', 'PSC', 'PS', 'ESC', 'ES']
    ahash = {
            'B6 untreated':0,
            'Pre-symptomatic':3,
            'Symptomatic':3,
            'End-stage':5,
            'B6 End-stage':5,
            'Pre-symptomatic control':2,
            'Symptomatic control':2,
            'WTSOD1 symptomatic control':2,
            'End-stage control':4,
            'WTSOD1 end-stage control':4,
            'LPS injected 48 hr timepoint':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChiu2013 = getChiu2013


def getGrabert2016(self):
    self.prepareData("MAC35")
    atype = self.h.getSurvName("c age")
    ahash = {'22month':2, '12month':1, '4month':0}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Title")
    atype = [ " ".join(str(i).split(" ")[0:2]) if \
            len(str(i).split(" ")) > 3 else None for i in atype]
    atypes = ['CeM', 'CoM', 'HiM', 'StM',
            'CeH', 'CoH', 'HiH', 'StH']
    ahash = {'Cerebellar microglia':0,
            'Cortical microglia':1,
            'Hippocampal microglia':2,
            'Striatal microglia':3,
            'Cerebellar homogenates':4,
            'Cortical homogenates':5,
            'Hippocampal homogenates':6,
            'Striatal homogenates':7}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGrabert2016 = getGrabert2016


def getMatcovitchNatan2016(self):
    self.prepareData("MAC36")
    atype = self.h.getSurvName("c developmental stage")
    atypes = ['E', 'N', 'A']
    ahash = {'E11.5':0, 'E12.5':0, 'E13.5':0, 'E14.5':0, 'E16.5':0, 'E10.5':0,
            'newborn':1, 'day3':1, 'day6':1, 'day9':1, '4wk':2, '5wk':2, '6wk':2,
            '8wk':2 }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatcovitchNatan2016 = getMatcovitchNatan2016


def getCho2019(self):
    self.prepareData("MAC37")
    atype = self.h.getSurvName("c Treatment")
    atypes = ['M0', 'M1', 'M2', 'IL10', 'M1+M2']
    ahash = {'LPS_lo + IFNg':1,
            'IL4':2,
            'None': 0,
            'LPS_lo + IL4 + IL10 + IL13':4,
            'IL10':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCho2019 = getCho2019


def getKrasemann2017(self):
    self.prepareData("MAC38")
    atype = self.h.getSurvName("c Title")
    atype = [str(i)[0:-3] for i in atype]
    atypes = ['WP', 'WN', 'AP', 'AN', 'SK', 'SH', 'ACp', 'ACn', 'WCn']
    ahash = {
            'WT_Phagocytic':0,
            'WT_NonPhagocytic':1,
            'Apoe knock-out_Phagocytic':2,
            'Apoe knock-out_NonPhagocytic':3,
            'SOD1:TREM2-KO (Female)':4,
            'SOD1:TREM2-KO_(Male)':4,
            'SOD1:TREM2-Het (Female)':5,
            'SOD1:TREM2-Het (Male)':5,
            'APP-PS1_Clec7apositive':6,
            'APP-PS1_Clec7anegative':7,
            'WT_Clec7anegative':8}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKrasemann2017 = getKrasemann2017


def getZhang2013(self, dbid = "MAC39.1"):
    self.dbid = dbid
    atype = self.h.getSurvName("c disease")
    atypes = ['A', 'N']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2013 = getZhang2013


def getCoates2008(self):
    self.prepareData("MAC40")
    atype = self.h.getSurvName("c Title")
    atype = [ str(i)[0:-2] for i in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'C57BL/6 4Gy':2, 'CBA/Ca 0Gy':0, 'CBA/Ca 4Gy':1,
            'C57BL/6 0Gy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCoates2008 = getCoates2008


def getHutchins2015(self):
    self.prepareData("MAC41")
    atype = self.h.getSurvName("c Replicate")
    ahash = {'1':1, '2':2}
    kval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Title")
    atype = [str(i).split(" ")[0] for i in atype]
    atypes = ['Mast', 'N', 'S', 'E', 'Mac']
    ahash = {'Mast':0,
            'Neutrophil':1,
            'Splenic':2,
            'Eosinophils':3,
            'PEC':4,
            'na\xc3\xafve':5}
    ahash = asciiNorm(ahash)
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['None', 'IL10', 'IL10, LPS', 'LPS']
    ahash = {}
    atype = [atype[i] if kval[i] == 2 else None for i in range(len(atype))]
    self.rval = rval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHutchins2015 = getHutchins2015


def getHutchins2015TPM(self, tn = 1):
    self.prepareData("MAC41.2")
    atype = self.h.getSurvName("c source_name")
    ahash = {'Bone marrow neutrophil':2,
            'spleen-purified dendritic cells':3,
            'Peritoneal exudate cells (adherent cells)':4,
            'Bone marrow-derived Eosinophils':5,
            'Bone marrow-derived mast cell':6,
            'CD4+ na\xc3\xafve T cells':7}
    ahash = asciiNorm(ahash)
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'LPS':1, 'IL10':2,'IL10\\, LPS':1, 'None':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn >= 2):
        #atype = [atype[i] if rval[i] == tn or aval[i] == 0
        #        else None for i in range(len(atype))]
        atype = [atype[i] if rval[i] == tn
                else None for i in range(len(atype))]
    self.rval = rval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHutchins2015TPM = getHutchins2015TPM


def getMo2018(self, tn=1):
    self.prepareData("PLP15")
    atype = self.h.getSurvName("c disease state (diagnosis)")
    atypes = ['Normal', 'UC', 'CD', 'sJ', 'oJ', 'pJ']
    ahash = {'Control':0,
            'Ulcerative Colitis':1,
            "Crohn's Disease":2,
            'Systemic JIA':3,
            'Oligoarticular JIA':4,
            'Polyarticular JIA':5}
    if (tn == 2):
        atypes = ['Normal', 'UC', 'CD']
        ahash = {'Control':0,
                'Ulcerative Colitis':1,
                "Crohn's Disease":2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMo2018 = getMo2018


def getBurczynski2006(self):
    self.prepareData("PLP13")
    atype = self.h.getSurvName("c Disease")
    atypes = ['Normal', 'UC', 'CD']
    ahash = {'Ulcerative Colitis':1, 'Normal':0, "Crohn's Disease":2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurczynski2006 = getBurczynski2006


def getPlanell2017(self):
    self.prepareData("PLP17")
    atype = self.h.getSurvName("c case_phenotype")
    atypes = ['Normal', 'UC', 'CD']
    ahash = {'Crohn':2, 'Colitis':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPlanell2017 = getPlanell2017


def getMehraj2013(self, tn=1):
    self.prepareData("MAC42")
    atype = self.h.getSurvName("c cell type")
    ahash = {'Monocyte':0, 'Monocytes-Derived Macrophages':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'NS':0, 'IFNg':1, 'IL4':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    atype = [atype[i] if rval[i] == 1
            else None for i in range(len(atype))]
    if (tn == 2):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub(" [0-9]*$", "", str(k)) for k in atype]
        ahash = {'Monocyte NS 6h', 'Monocyte IFNG 6h', 'Macrophage IL4 18h',
                'Monocyte IFNG 18h', 'Monocyte IL4 18h', 'Macrophage IFNG 18h',
                'Macrophage NS 18h', 'Monocyte IL4 6h', 'Monocyte NS 18h'}
        atypes = ['Mono', 'Mac']
        ahash = {'Macrophage NS 18h':1, 'Monocyte NS 6h':0, 'Monocyte NS 18h':0}
    if (tn == 3):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub(" [0-9]*$", "", str(k)) for k in atype]
        ahash = {'Monocyte NS 6h':0, 'Monocyte IFNG 6h':1,
                'Monocyte IFNG 18h':2, 'Monocyte IL4 18h':2,
                'Monocyte IL4 6h':2, 'Monocyte NS 18h':0}
        atypes = ['M0', 'M1', 'M2']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMehraj2013 = getMehraj2013


def getReynier2012(self):
    self.prepareData("MAC43")
    atype = self.h.getSurvName("c treatment")
    atypes = ['saline', 'LPS']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getReynier2012 = getReynier2012


def getChamberland2009(self):
    self.prepareData("MAC44")
    atype = self.h.getSurvName("c src1")
    atypes = ['H', 'Asthma']
    ahash = {'control_bronchial biopsies':0,
            'allergic_asthmatic_bronchial biopsies':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChamberland2009 = getChamberland2009


def getNascimento2009(self):
    self.prepareData("MAC45")
    atype = self.h.getSurvName("c src1")
    atypes = ['ND', 'DF', 'DHF']
    ahash = {'PBMCs from DF patient':1,
            'PBMCs from DHF patient':2,
            'PBMCs from ND patient':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNascimento2009 = getNascimento2009


def getWu2013(self):
    self.prepareData("MAC46")
    atype = self.h.getSurvName("c disease status")
    atypes = ['HC', 'R', 'NR']
    ahash = {'virological failure':2,
            'HIV sero-negative healthy controls':0,
            'sustained virus suppression':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWu2013 = getWu2013


def getWong2009(self, tn=1):
    self.prepareData("MAC47")
    ctype = self.h.getSurvName("c disease state")
    ttype = self.h.getSurvName("c outcome")
    atype = [ str(ctype[i]) + " " + str(ttype[i]) for i in
                            range(len(ctype))]
    atypes = ['HC', 'Alive', 'Dead']
    ahash = {'septic shock patient Survivor':1,
            'septic shock patient Nonsurvivor':2,
            'normal control Survivor':0,
            'normal control n/a':0}
    if (tn == 2):
        atypes = ['HC', 'Alive', 'Dead']
        ahash = {'septic shock patient Survivor':1,
                'septic shock patient Nonsurvivor':2,
                'normal control Survivor':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWong2009 = getWong2009


def getYoon2011(self, tn=1):
    self.prepareData("MAC48")
    atype = self.h.getSurvName("c disease state")
    atypes = ['Normal', 'mild', 'severe', 'THP1']
    ahash = {'asthma severe':2, 'asthma mild':1, '':3, 'normal':0}
    if (tn == 2):
        atypes = ['H', 'Asthma']
        ahash = {'asthma severe':1, 'asthma mild':1, 'normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYoon2011 = getYoon2011


def getVoraphani2014(self, tn=1):
    self.prepareData("MAC49")
    atype = self.h.getSurvName("c disease state")
    atypes = ['Control', 'Moderate', 'Severe']
    ahash = {'Severe Asthma':2, 'Moderate Asthma':1, 'Control':0}
    if (tn == 2):
        atypes = ['H', 'Asthma']
        ahash = {'Severe Asthma':1, 'Moderate Asthma':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVoraphani2014 = getVoraphani2014


def getLund2003(self):
    self.prepareData("MAC50.1")
    atype = self.h.getSurvName("c Cell Type")
    atypes = ['NT1', 'NT2', 'IL4', 'IL12', 'IL12+TGFb', 'IL4+TGFb']
    ahash = {'antiCD3+antiCD28+IL4':2,
            'antiCD3+antiCD28':1,
            'antiCD3+antiCD28+IL12':3,
            'no treatment':0,
            'antiCD3+antiCD28+IL12+TGFbeta':4,
            'antiCD3+antiCD28+IL4+TGFbeta':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLund2003 = getLund2003


def getLund2003II(self):
    self.prepareData("MAC50.2")
    atype = self.h.getSurvName("c Cell Type")
    atypes = ['NT1', 'NT2', 'IL4', 'IL12', 'IL12+TGFb', 'IL4+TGFb']
    ahash = {'antiCD3+antiCD28+IL4':2,
            'antiCD3+antiCD28':1,
            'antiCD3+antiCD28+IL12':3,
            'no treatment':0,
            'antiCD3+antiCD28+IL12+TGFbeta':4,
            'antiCD3+antiCD28+IL4+TGFbeta':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLund2003II = getLund2003II


def getYanez2019(self):
    self.prepareData("MAC51")
    atype = self.h.getSurvName("c culture condition")
    atypes = ['NT', 'Th0', 'Th1', 'Th2', 'Act']
    ahash = {
            'Stimulated under Th0 for 20 hours':1,
            'Stimulated under Th0 for 12 hours':1,
            'Stimulated under Th2 for 8 hours':3,
            'Stimulated for 16 hours and than Actinomysin added':4,
            'Stimulated under Th2 for 20 hours':3,
            'Purified Unstimulated CD4 T cells':0,
            'Stimulated under Th0 for 4 hours':1,
            'Stimulated under Th1 for 4 hours':2,
            'Stimulated under Th1 for 12 hours':2,
            'Stimulated under Th0 for 16 hours':1,
            'Stimulated under Th2 for 12 hours':3,
            'Stimulated under Th0 for 8 hours':1,
            'Stimulated under Th0 for 24 hours':1,
            'Stimulated under Th2 for 24 hours':3,
            'Stimulated under Th1 for 16 hours':2,
            'Stimulated under Th1 for 20 hours':2,
            'Stimulated under Th2 for 4 hours':3,
            'Stimulated under Th1 for 8 hours':2,
            'Stimulated under Th2 for 16 hours':3,
            'Stimulated under Th1 for 24 hours':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYanez2019 = getYanez2019


def getSpurlock2015(self):
    self.prepareData("MAC52")
    atype = self.h.getSurvName("c polarizing conditions")
    atypes = ['TH1', 'TH2', 'TH17']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSpurlock2015 = getSpurlock2015


def getMicosse2019(self):
    self.prepareData("MAC53")
    atype = self.h.getSurvName("c Characteristics[cell type]")
    atypes = ['TH1', 'TH2', 'TH17', 'TH9']
    ahash = {'T-helper 1 cell':0,
            'T-helper 17 cell':2,
            'T-helper 2 cell':1,
            'T-helper 9 cell':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMicosse2019 = getMicosse2019


def getGonzalezLeal2014(self, tn=1):
    self.prepareData("MAC54")
    atype = self.h.getSurvName("c treatment")
    atypes = ['CLIK', 'DMSO', 'CA074', 'pTH1', 'aTH1', 'pTH2', 'aTH2']
    ahash = {'pro TH2':5,
            'anti TH2':6,
            'CLIK control':0,
            'DMSO control':1,
            'CA074 control':2,
            'pro TH1':3,
            'anti TH1':4}
    if (tn == 2):
        atypes = ['C', 'TNFa']
        ahash = {'pro TH2':1,
                'CLIK control':0,
                'DMSO control':0,
                'CA074 control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGonzalezLeal2014 = getGonzalezLeal2014


def getBelmont(self, tn=1):
    self.prepareData("PLP41")
    atype = self.h.getSurvName("c genotype")
    atypes = ['Normal', 'Adenoma', 'Cancer']
    atypes = ['Normal', 'Adenoma']
    ahash = {'wild type':0, 'APC mutant': 1, 'APC KRAS mutant': 2}
    ahash = {'wild type':0, 'APC mutant': 1}
    if (tn == 2):
        atype = self.h.getSurvName("c src1")
        atypes = ['Normal', 'Cancer']
        ahash = {'primary colon cancer':1, 'normal colon control':0}
    if (tn == 3):
        atype = self.h.getSurvName("c genotype")
        atypes = ['Normal', 'Cancer']
        ahash = {'wild type':0, 'APC BRAF mutant':1,
                'APC BRAF P53 mutant':1}
    if (tn == 4):
        atype = self.h.getSurvName("c genotype")
        batch = self.h.getSurvName("c batch")
        atype = [str(atype[i])+ " " + str(batch[i]) for i in
                range(len(atype))]
        atypes = ['Normal', 'Cancer']
        ahash = {'wild type 3':0, 'APC KRAS mutant 3': 1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    expg = [ i for i in self.h.aRange() if aval[i] is not None]
    self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
    self.ad = [ i for i in self.h.aRange() if aval[i] == 1]
    self.aval = aval
    self.atype = atype
    self.atypes = atypes
    self.order = expg
    self.ibd = self.normal + self.ad
    self.printInfo()

def getGerling2016(self, tn=1):
    self.prepareData("PLP42")
    atype = self.h.getSurvName("c genotype")
    atypes = ['td', 'fl/+', 'fl/fl', 'fl/fl td']
    ahash = {'R26-tdTomato':0,
            'Ptch1fl/+':1,
            'Col1a2CreER; Ptch1fl/fl; R26-LSL-tdTomato':3,
            'Col1a2CreER; Ptch1fl/fl':2}
    if (tn == 2):
        atype = self.h.getSurvName("c Title")
        atype = [re.split("[ _]", str(i))[0] for i in atype]
        atypes = ['wt', 'Hh']
        ahash = {'wt': 0, 'Hh':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    expg = [ i for i in self.h.aRange() if aval[i] is not None]
    self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
    self.ad = [ i for i in self.h.aRange() if aval[i] == 1]
    self.aval = aval
    self.atype = atype
    self.atypes = atypes
    self.order = expg
    self.ibd = self.normal + self.ad
    self.printInfo()

def getMcNeil2016(self, tn=1):
    self.prepareData("PLP43")
    atype = self.h.getSurvName("c Stage")
    atypes = ['U', 'P', 'I', 'LT']
    ahash = {'p52':3, 'p5':1, 'p8':2, 'p51':3, 'p0':0, 'p3':1,
            'p66':3, 'p6':2, 'p54':3, 'p7':2, 'p55':3, 'p58':3}
    if (tn == 2):
        atypes = ['U', 'LT']
        ahash = {'p52':1, 'p51':1, 'p0':0,
                'p66':1, 'p54':1, 'p55':1, 'p58':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMcNeil2016 = getMcNeil2016


def getNeufert2013(self, tn=1):
    self.prepareData("PLP44")
    atype = self.h.getSurvName("c src1")
    atypes = ['CS', 'CC', 'TS', 'TC']
    ahash = {'colorectal control epithelium_colitis-associated':1,
            'colorectal tumor_colitis-associated':2,
            'colorectal control epithelium_sporadic':0,
            'colorectal tumor_sporadic':3}
    if (tn == 2):
        atypes = ['N', 'T']
        ahash = {'colorectal control epithelium_colitis-associated':0,
                'colorectal tumor_colitis-associated':1,
                'colorectal control epithelium_sporadic':0,
                'colorectal tumor_sporadic':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNeufert2013 = getNeufert2013


def getLeclerc2004(self):
    self.prepareData("PLP45")
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'T']
    ahash = {'Normal intestine':0, 'Tumor':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLeclerc2004 = getLeclerc2004


def getZhu2016(self):
    self.prepareData("PLP46")
    atype = self.h.getSurvName("c sample type")
    atypes = ['N', 'T']
    ahash = {'Mouse tumor':1, 'Mouse primary cells':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhu2016 = getZhu2016


def getPaoni(self, tn=1):
    self.prepareData("PLP51")
    atype = self.h.getSurvName("c Title")
    atypes = ['Normal', 'Adenoma', 'Cancer']
    atypes = ['Normal', 'Adenoma']
    if (tn == 2):
        atypes = ['Normal', 'Adenoma', 'Cancer']
    ahash = {}
    aval = [ahash[i] if i in ahash else None for i in atype]
    normal = [ i for i in self.h.aRange() if atype[i].find("-WT") > 0]
    for i in normal:
        aval[i] = 0
    ad = [ i for i in self.h.aRange() if atype[i].find("-adenoma") > 0]
    for i in ad:
        aval[i] = 1
    ca = [ i for i in self.h.aRange() if atype[i].find("-carc") > 0]
    for i in ca:
        aval[i] = 2
    self.aval = aval
    self.atype = atype
    self.atypes = atypes
    self.order = normal + ad
    if (tn == 2):
        self.order = normal + ad + ca
    self.normal = normal
    self.ad = ad
    self.ca = ca
    self.ibd = self.normal + self.ad
    self.printInfo()

def getKaiser2007(self):
    self.prepareData("PLP54")
    atype = self.h.getSurvName("c src1")
    atypes = ['b', 'C']
    ahash = {' b':0, '3 b':0, '2 b':0, 'Colon tissue':1, '1 b':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKaiser2007 = getKaiser2007


def getNoble2010(self, tn=1):
    self.prepareData("PLP58")
    atype = self.h.getSurvName("c src1")
    atypes = ['A', 'D', 'S', 'T', 'A_CD', 'D_CD', 'S_CD', 'T_CD']
    ahash = {'ascending colon biopsy from healthy subject':0,
            'descending colon biopsy from healthy subject':1,
            'sigmoid colon biopsy from healthy subject':2,
            'terminal ileum biopsy from healthy subject':3,
            'ascending colon biopsy from crohns disease subject':4,
            'descending colon biopsy from crohns disease subject':5,
            'sigmoid colon biopsy from crohns disease subject':6,
            'terminal ileum biopsy from crohns disease subject':7}
    if (tn == 2):
        atypes = ['A', 'D']
        ahash = {'ascending colon biopsy from healthy subject':0,
                'descending colon biopsy from healthy subject':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNoble2010 = getNoble2010


def getColonGEO(self):
    self.prepareData("CRC115")
    self.normal = normal
    self.ad = ad
    self.ca = ca
    self.ibd = self.normal + self.ad
    self.printInfo()

def getColonGEO(self):
    self.prepareData("CRC115")
    atype = self.h.getSurvName("c Histology")
    atypes = ['Normal', 'Adenoma', 'Carcinoma']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getColonGEO = getColonGEO


def getSchridde2017(self):
    self.prepareData("MAC55")
    atype = self.h.getSurvName("c desc")
    atype = [str(i).split(" ")[4] if len(str(i).split(" ")) > 5 \
            else None for i in atype]
    atypes = ['P1', 'P2', 'P3', 'P4']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSchridde2017 = getSchridde2017


def getQu2016(self, tn=1):
    self.prepareData("PLP50")
    atype = self.h.getSurvName("c tissue type")
    atypes = ['NS', 'NC', 'A', 'C', 'M']
    ahash = {'Metastasis':4,
            'Carcinoma':3,
            'Normal crypt epithelium':1,
            'Adenoma':2,
            'Normal surface epithelium':0}
    if (tn == 2):
        atypes = ['NC', 'NS']
        ahash = {'Normal crypt epithelium':0,
                'Normal surface epithelium':1}
    if (tn == 3):
        atypes = ['NC', 'A']
        ahash = {'Normal crypt epithelium':0, 'Adenoma':1}
    if (tn == 4):
        atypes = ['NC', 'C']
        ahash = {'Normal crypt epithelium':0, 'Carcinoma':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQu2016 = getQu2016


def getGalamb2008(self):
    self.prepareData("PLP83")
    atype = self.h.getSurvName("c desc")
    atype = [ " ".join(str(i).split(" ")[-2:]) for i in atype]
    atypes = ['N', 'IBD', 'A', 'C']
    ahash = {'bowel disease':1,
            'colorectal cancer':3,
            'colon adenoma':2,
            'healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGalamb2008 = getGalamb2008


def getArendt2015(self):
    self.prepareData("LIV4")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['HC', 'FL', 'SH']
    ahash = {'SS':1, 'NASH':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArendt2015 = getArendt2015


def getTung2011(self):
    self.prepareData("LIV1")
    atype = self.h.getSurvName("c src1")
    atypes = ['H', 'NT', 'C', 'HCC']
    ahash = {'non_tumor':1, 'tumor':3, 'cirrhotic':2, 'healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTung2011 = getTung2011


def getBaker2010(self):
    self.prepareData("LIV7")
    atype = self.h.getSurvName("c disease state")
    atypes = ['H', 'SH']
    ahash = {'non-alcoholic steatohepatitis (NASH)':1,
            'normal (control)':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBaker2010 = getBaker2010


def getAlao2016(self):
    self.prepareData("LIV16")
    atype = self.h.getSurvName("c treatment")
    atypes = ['Base', 'Wk2', 'Wk4']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAlao2016 = getAlao2016


def getMeissner2016HCV(self, tn=1):
    self.prepareData("LIV17")
    atype = self.h.getSurvName("c time point")
    ahash = {'pre-treatment':0, 'post-treatment':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c tissue")
    ahash = {'blood':0, 'liver':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['none', 'simtuzumab']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeissner2016HCV = getMeissner2016HCV


def getMeissner2016Cir(self, tn=1):
    self.prepareData("LIV18")
    atype = self.h.getSurvName("c treatment time")
    atypes = ['Pre', 'Post-DAA']
    ahash = {'EOT':1, 'PRE':0}
    if (tn == 2):
        paired = self.h.getSurvName('c paired mate')
        paired = [re.sub("sample ", "", str(k)) for k in paired]
        phash = {'1', '5', '6', '3', '11', '15', '16', '13'}
        atype = [atype[i] if paired[i] in phash else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeissner2016Cir = getMeissner2016Cir


def getPeters2017(self, dtype=0):
    self.prepareData("PLP18")
    atype = self.h.getSurvName("c src1")
    ahash = {'Transverse colon':0, 'Ascending colon':0, 'Rectum':0,
            'Descending colon':0, 'Sigmoid colon':0, 'Normal':0,
            'Terminal Ileum':1, 'Blood':2,
            'DELETED':3, 'Not Available':3, 'Not Collected':3}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c inflammation")
    atypes = ['Normal', 'NI', 'I', 'NA']
    ahash = {'Non-involved area':1, 'NA':3, 'na':3, 'Inflamed area':2,
            'Normal':0, 'DELETED':3, 'N/A (Sample not received)':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeters2017 = getPeters2017


def getKF(self):
    self.prepareData("CRC11")
    kras = self.h.getSurvName('c KRAS Mutation')
    khash = {'WT':1, 'c.35G>T':2, 'c.35G>A':2, 'c.35G>C':2,
            'c.38G>A':2, 'c.34G>A':2}
    kras_m = [khash[i] if i in khash else 0 for i in kras]
    atype = self.h.getSurvName('c Best Clinical Response Assessment')
    atypes = ['DCG', 'PD']
    ahash = {"CR": 0, "PR":0, "SD":0, "CR + PR":0, "PD":1}
    self.rval = kras_m
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKF = getKF


def getDelRoy(self):
    self.prepareData("CRC101")
    atype = self.h.getSurvName('c Status')
    atypes = ['R', 'NR']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelRoy = getDelRoy


def getLissner(self, tn=1):
    self.prepareData("MAC56")
    atype = self.h.getSurvName('c Agent')
    ahash = {'Lm':0, 'LPS':1}
    bval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Timepoint')
    ahash = {'6hr':6, '2hr':2, '1hr':1, '0hr':0}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Type')
    atypes = ['N', 'A', 'O']
    ahash = {'Neonate':0, 'Adult':1, 'OlderAdult':2}
    if (tn == 2):
        atype = [atype[i] if bval[i] == 0 and rval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if bval[i] == 1 and rval[i] <= 1 else None
                for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if bval[i] == 0 and rval[i] == 0 else None
                for i in range(len(atype))]
        atypes = ['A', 'N']
        ahash = {'Neonate':1, 'Adult':0}
    if (tn == 5):
        atype = [atype[i] if bval[i] == 0 and rval[i] == 0 else None
                for i in range(len(atype))]
        atypes = ['A', 'O']
        ahash = {'OlderAdult':1, 'Adult':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLissner = getLissner


def getHennenberg(self):
    self.prepareData("PLP100")
    atype = self.h.getSurvName('c src1')
    atypes = ['WT', 'KO']
    ahash = {'MDR1A KO, Tumor, Colon, AOM/DSS':1,
            'WT, Tumor, Colon, AOM/DSS':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHennenberg = getHennenberg


def getTakahara2015(self, ri=0):
    self.prepareData("PLP47")
    atype = self.h.getSurvName('c genotype')
    ahash = {'WT':0, 'ATF6bKO':1, 'ATF6aKO':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['Un', 'DSS']
    ahash = {'DSS for 3 days':1, 'untreated':0}
    atype = [atype[i] if rval[i] == ri else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTakahara2015 = getTakahara2015


def getDevHs(self):
    self.prepareData("PLP103")
    atype = self.h.getSurvName('c Type')
    atypes = ['YS', 'Pl', 'ES', 'ST', 'SI', 'CO']
    ahash = {'yolk sac':0, 'placenta':1,
            'esophagus':2, 'stomach':3,
            'small intestine':4, 'colon':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDevHs = getDevHs


def getDevMm(self):
    self.prepareData("PLP104")
    atype = self.h.getSurvName('c Type')
    atypes = ['YS', 'E13.5 Ca', 'E13.5 I', 'E13.5 Co']
    ahash = {'yolk sac':0, 'Caecum':1, 'Ileum':2, 'Colon':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDevMm = getDevMm


def getSmith2009(self, ri=0):
    self.prepareData("MAC57")
    atype = self.h.getSurvName('c Sample Characteristic[stimulus]')
    ahash = {'none':0,
            'heat-killed Escherichia coli':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Sample Characteristic[disease]')
    atypes = ['N', 'iCD', 'cCD']
    ahash = {'normal':0,
            "ileal Crohn's disease":1,
            "colonic Crohn's disease":2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2009 = getSmith2009


def getDeSchepper2018(self):
    self.prepareData("MAC58")
    atype = self.h.getSurvName('c Factor Value[organism part]')
    ahash = {'lamina propria of small intestine':0,
            'muscularis externa layer of small intestine':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Characteristics[phenotype]')
    atypes = ['LPMR', 'LPSM', 'MEMR', 'MESM']
    ahash = {'YFP negative':0, 'YFP positive':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeSchepper2018 = getDeSchepper2018


def getGrover2018(self):
    self.prepareData("PLP105")
    atype = self.h.getSurvName('c sample type')
    atypes = ['NDNG', 'DNG', 'DG', 'IG']
    ahash = {'diabetic gastroparetics':2,
            'non-diabetic non-gastroparetic':0,
            'diabetic non-gastroparetics controls':1,
            'idiopathic gastroparetics':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGrover2018 = getGrover2018


def getXu2014(self):
    self.prepareData("GS2.7")
    atype = self.h.getSurvName('c src1')
    atypes = ['LGD', 'HGD', 'inf', 'EGC']
    ahash = {'inflammation':2, 'HGD':1, 'EGC':3, 'LGD':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXu2014 = getXu2014


def getBadawi2019(self, tn=1):
    self.prepareData("MHP5")
    atype = self.h.getSurvName('c time post surgery')
    atypes = ['0', '45m', '24h']
    ahash = {'45 minutes':1, '24 hours':2, '0 minutes':0}
    if (tn == 2):
        atypes = ['0', '24h']
        ahash = {'45 minutes':0, '24 hours':1, '0 minutes':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBadawi2019 = getBadawi2019


def getRock2005(self):
    self.prepareData("MAC73")
    atype = self.h.getSurvName('c Group')
    atypes = ['Control', 'IFNG']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRock2005 = getRock2005


def getPolak2014(self, tn=1):
    self.prepareData("MAC69")
    atype = self.h.getSurvName('c treatment')
    ahash = {'TNF-alpha':0, 'TSLP':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c time')
    atypes = ['0', '2', '8', '24']
    ahash = {'8h':2, '2h':1, '24h':3, '0h':0}
    if (tn == 2):
        atypes = ['C', 'TNF-a']
        ahash = {'8h':1, '24h':1, '0h':0}
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPolak2014 = getPolak2014


def getSchirmer2009(self, tn=1):
    self.prepareData("KR5")
    atype = self.h.getSurvName('c Disease Status')
    ahash = {'control':0, 'patient':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c src1')
    atypes = ['RM', 'T', 'SM', 'SC', 'M']
    ahash = {'resting monocyes':0,
            'T cells':1,
            'stimulated monocyes':2,
            'stem cells':3,
            'macrophages':4}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSchirmer2009 = getSchirmer2009


def getPritchard2011(self):
    self.prepareData("KR7")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'Ath']
    ahash = {'Male Subject with Carotid Disease':1, 'Male Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPritchard2011 = getPritchard2011


def getHagg2008(self):
    self.prepareData("KR1")
    atype = self.h.getSurvName('c Sample Type')
    atypes = ['MC', 'MD', 'FC', 'FD']
    ahash = {'Baseline macrophages without atherosclerosis':0,
            'Baseline macrophages with atherosclerosis':1,
            'Foam cells with atherosclerosis':3,
            'Foam cells without atherosclerosis':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHagg2008 = getHagg2008


def getBondar2017(self, tn=1):
    self.prepareData("MAC76")
    atype = self.h.getSurvName('c gender')
    ahash = {'Female':0, 'Male':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c underlying disease')
    atypes = ['NICM', 'ICM', 'PPCM', 'NCIM', 'ChemoCM']
    if (tn == 2):
        atypes = ['NICM', 'ICM']
    if (tn == 3):
        atypes = ['NICM', 'ICM']
        atype = [atype[i] if tval[i] == 1 else None
            for i in range(len(atype))]
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBondar2017 = getBondar2017


def getMaciejak2015(self, tn=1):
    self.prepareData("MAC75")
    atype = self.h.getSurvName('c samples collection')
    ahash = {'on the 1st day of MI (admission)':1,
            '1 month after MI':3, 'after 4-6 days of MI (discharge)':2,
            '6 months after MI':4, 'N/A':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c hf progression')
    atypes = ['non-HF', 'HF', 'stable CAD']
    if (tn == 2):
        atypes = ['non-HF', 'HF']
    if (tn == 3):
        atypes = ['non-HF', 'HF']
        atype = [atype[i] if tval[i] == 3 else None
            for i in range(len(atype))]
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaciejak2015 = getMaciejak2015


def getNorona2019(self, tn=1):
    self.prepareData("MAC78")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'TGF', 'MTX']
    ahash = {'TGF':1, 'MTX':2, 'VEH':0}
    if (tn == 2):
        atypes = ['C', 'TGF']
        ahash = {'TGF':1, 'VEH':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNorona2019 = getNorona2019


def getEijgelaar2010(self, tn=1):
    self.prepareData("MAC79")
    atype = self.h.getSurvName('c Sample')
    atypes = ['LI', 'LU', 'SP', 'CA']
    if (tn == 2):
        atypes = ['CA', 'LI']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEijgelaar2010 = getEijgelaar2010


def getKeller2009(self, tn=1):
    self.prepareData("MAC80")
    atype = self.h.getSurvName('c time')
    atype = [str(k).split(" ")[-1] for k in atype]
    atypes = ['CT0', 'CT4', 'CT8', 'CT12', 'CT16', 'CT20']
    if (tn == 2):
        atype = self.h.getSurvName('c time')
        atype = [str(k).split(" ")[0] for k in atype]
        atypes = ['first', 'second']
    ahash = {}
    if (tn == 3):
        atype = self.h.getSurvName('c time')
        atypes = ['f0', 'f4', 'f8', 'f12', 'f16', 'f20',
                's0', 's4', 's8', 's12', 's16', 's20']
        ahash = {'first day CT0':0,
                'first day CT4':1,
                'first day CT8':2,
                'first day CT12':3,
                'first day CT16':4,
                'first day CT20':5,
                'second day CT0':6,
                'second day CT4':7,
                'second day CT8':8,
                'second day CT12':9,
                'second day CT16':10,
                'second day CT20':11}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKeller2009 = getKeller2009


def getImmGenULI(self, tn=1):
    self.prepareData("MAC81")
    atype = self.h.getSurvName('c Title')
    atypes = ['M0', 'M1']
    ahash = {}
    if (tn == 1):
        ahash = {'MF.11cpSigFp.BAL.1':0,
                'MF.11cpSigFp.BAL.2':0,
                'MF.SSChipSigFn.LPS.d3.BAL.1':1,
                'MF.SSChipSigFn.LPS.d3.BAL.2':1,
                'MF.11cpSigFp.LPS.d6.BAL.1':1,
                'MF.11cpSigFp.LPS.d6.BAL.2':1,
                'MF.SSChipSigFn.LPS.d6.BAL.1':1,
                'MF.SSChipSigFn.LPS.d6.BAL.2':1}
        ahash = {'MF.11cpSigFp.BAL.1':0,
                'MF.11cpSigFp.BAL.2':0,
                'MF.SSChipSigFn.LPS.d3.BAL.1':1,
                'MF.SSChipSigFn.LPS.d3.BAL.2':1}
    if (tn == 2):
        atypes = ['M0', 'M1']
        ahash = {'MF.64p6Cn206nIIp.LPS.d3.Lu.1':1,
                'MF.64p6Cn206nIIp.LPS.d3.Lu.2':1,
                'MF.64p6Cn206nIIp.LPS.d6.Lu.1':1,
                'MF.64p6Cn206nIIp.LPS.d6.Lu.2':1,
                'MF.64p6Cn206nIIp.Lu.1':0,
                'MF.64p6Cn206nIIp.Lu.2':0}
    if (tn == 3):
        atypes = ['M0', 'M1']
        ahash = {'MF.64p6Cn206pIIn.LPS.d3.Lu.1':1,
                'MF.64p6Cn206pIIn.LPS.d3.Lu.2':1,
                'MF.64p6Cn206pIIn.LPS.d6.Lu.1':1,
                'MF.64p6Cn206pIIn.LPS.d6.Lu.2':1,
                'MF.64p6Cn206pIIn.Lu.1':0,
                'MF.64p6Cn206pIIn.Lu.2':0}
    if (tn == 4):
        atypes = ['M0', 'M1']
        ahash = {'MF.F.10kIFN.PC#1':1,
                'MF.F.10kIFN.PC#2':1,
                'MF.F.10kIFN.PC#3':1,
                'MF.F.PC#1':0,
                'MF.F.PC#1_':0,
                'MF.F.PC#2':0,
                'MF.F.PC#2_':0,
                'MF.F.PC#3':0,
                'MF.F.PC.1':0,
                'MF.F.PC.2':0,
                'MF.F.PC.3':0,
                'MF.Fem.PC#1_RNA-seq':0,
                'MF.Fem.PC#2_RNA-seq':0}
    if (tn == 5):
        atypes = ['M0', 'M1']
        ahash = {'MF.M.10kIFN.PC#1':1,
                'MF.M.10kIFN.PC#2':1,
                'MF.M.10kIFN.PC#3':1,
                'MF.M.PC#1':0,
                'MF.M.PC#2':0,
                'MF.M.PC#3':0}
        for k in ['MF.PC#1.1', 'MF.PC#1.10', 'MF.PC#1.11', 'MF.PC#1.12',
                'MF.PC#1.13', 'MF.PC#1.14', 'MF.PC#1.15', 'MF.PC#1.16',
                'MF.PC#1.2', 'MF.PC#1.3', 'MF.PC#1.4', 'MF.PC#1.5', 'MF.PC#1.6',
                'MF.PC#1.7', 'MF.PC#1.8', 'MF.PC#1.9', 'MF.PC#1_RNA-seq',
                'MF.PC#2.1', 'MF.PC#2.2', 'MF.PC#2.3', 'MF.PC#2.4', 'MF.PC#2.5',
                'MF.PC#2.6', 'MF.PC#2.7', 'MF.PC#2.8', 'MF.PC#2_RNA-seq',
                'MF.PC#3', 'MF.PC#3_RNA-seq', 'MF.PC#4', 'MF.PC#4_RNA-seq',
                'MF.PC.01', 'MF.PC.02', 'MF.PC.03', 'MF.PC.04', 'MF.PC.05',
                'MF.PC.06', 'MF.PC.07', 'MF.PC.08', 'MF.PC.09', 'MF.PC.10',
                'MF.PC.11', 'MF.PC.12', 'MF.PC.13', 'MF.PC.14', 'MF.PC.15',
                'MF.PC.17', 'MF.PC.18', 'MF.PC.19', 'MF.PC.20', 'MF.PC.21',
                'MF.PC.23', 'MF.PC.24', 'MF.PC.25', 'MF.PC.26', 'MF.PC.37'
                'MF.PC.38', 'MF.PC.39', 'MF.PC.40']:
            ahash[k] = 0
    if (tn == 6):
        atypes = ['M0', 'M1']
        ahash = {'Mo.6Chi11bp.APAP.36h.Lv.1':1,
                'Mo.6Chi11bp.APAP.36h.Lv.2':1,
                'Mo.6Chi11bp.APAP.36h.Lv.3':1,
                'Mo.6Chi11bp.APAP.36h.Lv.4':1,
                'Mo.6Chi11bp.PBS.Lv.1':0,
                'Mo.6Chi11bp.PBS.Lv.2':0,
                'Mo.6Chi11bp.PBS.Lv.3':0,
                'Mo.6Chi11bp.PBS.Lv.4':0}
    if (tn == 7):
        atypes = ['M0', 'M1']
        ahash = {
                'NKT.F.Sp#1':0,
                'NKT.F.Sp#2':0,
                'NKT.F.Sp#3':0,
                'NKT.M.Sp#1':0,
                'NKT.M.Sp#2':0,
                'NKT.M.Sp#3':0,
                'NKT.Sp#3_RNA-seq':0,
                'NKT.Sp.LPS.18hr#1_RNA-seq':1,
                'NKT.Sp.LPS.18hr#2_RNA-seq':1,
                'NKT.Sp.LPS.3hr#1_RNA-seq':1,
                'NKT.Sp.LPS.3hr#2_RNA-seq':1}
    if (tn == 8):
        atypes = ['M0', 'M1']
        ahash = {
                'T4.F.10kIFN.Sp#1':1,
                'T4.F.10kIFN.Sp#2':1,
                'T4.F.10kIFN.Sp#3':1,
                'T4.F.Sp#1':0,
                'T4.F.Sp#2':0,
                'T4.F.Sp#3':0,
                'T4.M.10kIFN.Sp#1':1,
                'T4.M.10kIFN.Sp#2':1,
                'T4.M.10kIFN.Sp#3':1,
                'T4.M.Sp#1':0,
                'T4.M.Sp#2':0,
                'T4.M.Sp#3':0}
    if (tn == 9):
        atypes = ['M0', 'M1']
        ahash = {
                'B.17m.F.Sp#1':0,
                'B.20m.Sp#1':0,
                'B.2m.F.Sp#1':0,
                'B.2m.Sp#1':0,
                'B.6m.F.Sp#1':0,
                'B.6m.F.Sp#2':0,
                'B.6m.Sp#1':0,
                'B.6m.Sp#2':0,
                'B.F.10kIFN.Sp#1':1,
                'B.F.10kIFN.Sp#2':1,
                'B.F.10kIFN.Sp#3':1,
                'B.F.1kIFN.Sp#1':1,
                'B.F.1kIFN.Sp#2':1,
                'B.F.1kIFN.Sp#3':1,
                'B.F.Sp#1':0,
                'B.F.Sp#1_':0,
                'B.F.Sp#2':0,
                'B.F.Sp#2_':0,
                'B.F.Sp#3':0,
                'B.Fem.Sp#1_RNA-seq':0,
                'B.Fo.Sp#1_RNA-seq':0,
                'B.Fo.Sp#2_RNA-seq':0,
                'B.Fo.Sp#3_RNA-seq':0,
                'B.Fo.Sp#4_RNA-seq':0}
    if (tn == 10):
        atypes = ['M0', 'M1']
        ahash = {
                'GN.17m.F.Sp#1':0,
                'GN.20m.Sp#1':0,
                'GN.F.10kIFN.Sp#1':1,
                'GN.F.10kIFN.Sp#2':1,
                'GN.F.10kIFN.Sp#3':1,
                'GN.F.Sp#1':0,
                'GN.F.Sp#2':0,
                'GN.F.Sp#3':0,
                'GN.M.10kIFN.Sp#1':1,
                'GN.M.10kIFN.Sp#2':1,
                'GN.M.10kIFN.Sp#3':1,
                'GN.M.Sp#1':0,
                'GN.M.Sp#2':0,
                'GN.M.Sp#3':0,
                'GN.Sp#3_RNA-seq':0,
                'GN.Sp#4_RNA-seq':0}
    if (tn == 11):
        atypes = ['M0', 'M1', 'M2']
        ahash = {
                'MF.KC.Clec4FpTim4p64p.APAP.12h.Lv.1':2,
                'MF.KC.Clec4FpTim4p64p.APAP.12h.Lv.2':2,
                'MF.KC.Clec4FpTim4p64p.APAP.12h.Lv.4':2,
                'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.1':2,
                'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.2':2,
                'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.3':2,
                'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.4':2,
                'MF.KC.Clec4FpTim4p64p.Lv.2':0,
                'MF.KC.Clec4FpTim4p64p.Lv.3':0,
                'MF.KC.Clec4FpTim4p64p.Lv.4':0,
                'MF.KC.Clec4FpTim4p64p.PBS.Lv.1':0,
                'MF.KC.Clec4FpTim4p64p.PBS.Lv.2':0,
                'MF.KC.Clec4FpTim4p64p.PBS.Lv.3':0,
                'MF.KC.Clec4FpTim4p64p.PBS.Lv.4':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getImmGenULI = getImmGenULI


def getZigmond2014(self, tn=1):
    self.prepareData("MAC82")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['C', 'Il10', 'Il10ra']
    ahash = {'wild type; Cx3cr1gfp/+':0,
            'Interleukin-10 deficient; IL10-/- CX3CR1gfp/+':1,
            'macrophage-restricted interleukin-10 receptor deficient; CX3CR1cre:IL10Raflox/flox':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZigmond2014 = getZigmond2014


def getRamanan2016(self, tn=1):
    self.prepareData("PLP108")
    atype = self.h.getSurvName('c Title')
    gtype = [str(k).split(" ")[0] if len(str(k).split(" ")) > 0 else None
            for k in atype]
    ttype = [str(k).split(" ")[1] if len(str(k).split(" ")) > 1 else None
            for k in atype]
    atype = [ str(gtype[i]) + " " + str(ttype[i]) for i in
            range(len(atype))]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'Nod2-/-, IL-13':2,
            'Nod2-/-, untreated,':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRamanan2016 = getRamanan2016


def getMishra2019(self, tn=1):
    self.prepareData("MAC83")
    atype = self.h.getSurvName('c treatment')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'exposed to IL4/IL13':2,
            'exposed to LPS/IL4/IL13':0,
            'exposed to IFN gamma/LPS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMishra2019 = getMishra2019


def getGuler2015(self, tn=1):
    self.prepareData("MAC84.3")
    atype = self.h.getSurvName('c Title')
    atype = ["_".join(str(k).split("_")[0:-2]) 
            if len(str(k).split("_")) > 2 else None for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'IFNg':1,
            'IL4IL13':2,
            'M.tb_IL4IL13':0,
            'M.tb':1,
            'Ust':0,
            'M.tb_IFNg':1,
            'M.tb_IL41L13':0}
    if (tn == 2):
        atype = self.h.getSurvName('c Title')
        atype = ["_".join(str(k).split("_")[0:-1]) 
                if len(str(k).split("_")) > 2 else None for k in atype]
        ahash = {
                'Ust_28h':0,
                'IL4IL13_28h':2,
                'M.tb_IFNg_28h':1,
                'IFNg_28h':1,
                'M.tb_28h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuler2015 = getGuler2015


def getZhang2010(self, tn=1):
    self.prepareData("MAC85")
    atype = self.h.getSurvName('c Title')
    atype = [str(k).split("_")[0] if len(str(k).split("_")) > 0 else None
            for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'IFNG':1,
            'TNF':1,
            'L.':1,
            'Control':0,
            'IFNB':1,
            'IL4':2,
            'LPS':1,
            'T.':1,
            'IL10':2,
            'IL17':2}
    if (tn == 2):
        atype = self.h.getSurvName('c Title')
        atype = [re.sub("_rep.*", "", str(k)) for k in atype]
        ahash = {'IFNG_12h':1,
                'TNF_12h':1,
                'Control_6h':0,
                'IFNB_24h':1,
                'Control_2h':0,
                'Control_12h':0,
                'LPS_6h':1,
                'Control_24h':0,
                'IL4_12h':2,
                'LPS_2h':1,
                'IFNB_2h':1,
                'Control_0h_for_IL10':0,
                'Control_0h_for_Control':0,
                'TNF_6h':1,
                'IFNG_24h':1,
                'LPS_12h':1,
                'IL10_24h':2,
                'Control_0h_for_IL4':0,
                'Control_0h_for_IL17':0,
                'IL17_6h':2,
                'IFNB_12h':1,
                'Control_0h_for_TNF':0,
                'IL10_12h':2,
                'IL4_6h':2,
                'IFNB_6h':1,
                'IL4_2h':2,
                'LPS_24h':1,
                'IL17_12h':2,
                'IL17_2h':2,
                'IL10_2h_2nd_scan':2,
                'Control_0h_for_LPS':0,
                'IFNG_6h':1,
                'IL10_6h':2,
                'Control_0h_for_IFNG':0,
                'IFNG_2h':1,
                'LPS_12h_2nd_scan':1,
                'TNF_2h':1,
                'Control_0h_for_IFNB':0,
                'IL17_24h':2,
                'IL10_2h':2,
                'TNF_24h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2010 = getZhang2010


def getHealy2016(self, tn=1):
    self.prepareData("MAC86")
    atype = self.h.getSurvName('c src1')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'Human adult brain-derived microglia, M2a':2,
            'Human adult brain-derived microglia, M1':1,
            'Human adult brain-derived microglia, M2c':2,
            'Human adult brain-derived microglia, M0':0,
            'Human adult brain-derived microglia, Mtgf':2}
    if (tn == 2):
        ahash = {'Human adult peripheral blood-derived macrophages, M1':1,
                'Human adult peripheral blood-derived macrophages, M2a':2,
                'Human adult peripheral blood-derived macrophages, M0':0,
                'Human adult peripheral blood-derived macrophages, M2c':2,
                'Human adult peripheral blood-derived macrophages, Mtgf':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHealy2016 = getHealy2016


def getSvensson2011(self, tn=1):
    self.prepareData("MAC87")
    atype = self.h.getSurvName('c factors')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'M-CSF + E/P/IL-10/4/13':2, 'M-CSF + IL-10':2,
            'GM-CSF':0, 'M-CSF':0, 'M-CSF + IL4/13':2,
            'GM-CSF + M-CSF':0, 'GM-CSF/LPS/IFN':1, 'GM-CSF + IL4/13':2,
            'GM-CSF + IL-10':2, 'Decidual macrophages':0, 'Blood monocytes':0,
            'GM-CSF + M-CSF +IL-10':2, 'M-CSF + M-CSF':0, 'GM-CSF/LPS/IFN D6':1,
            'GM-CSF + E/P/IL-10/4/13':2}
    if (tn == 2):
        atypes = ['GMCSF', 'MCSF']
        ahash = {'GM-CSF':0, 'GM-CSF + M-CSF':0, 'M-CSF':1, 'M-CSF + M-CSF':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSvensson2011 = getSvensson2011


def getChandriani2014(self, tn=1):
    self.prepareData("MAC88")
    source = self.h.getSurvName('c src1')
    atype = self.h.getSurvName('c treatment')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'unstimulated':0, 'IL13':2, 'TGFb':1, 'IL10':2, 'IL4':2,
            'Dex':1}
    if (tn == 2):
        atype = source
        ahash = {'Monocytes, unstimulated, 24h':0, 'Monocytes, IL13, 24h':2,
                'Monocytes, unstimulated, 6h':0, 'Monocytes, IL4, 24h':2,
                'Monocytes, IL4, 6h':2, 'Monocytes, IL13, 6h':2}
    if (tn == 3):
        atype = source
        ahash = {'Macrophages, unstimulated, 24h':0, 'Macrophages, IL13, 24h':2,
                'Macrophages, IL10, 24h':2, 'Macrophages, TGFb, 24h':1,
                'Macrophages, IL4, 24h':2, 'Macrophages, Dex, 24h':1}
    if (tn == 4):
        atype = source
        ahash = {'Normal lung fibroblasts, TGFb, 24h':1,
                'Normal lung fibroblasts, IL13, 24h':2,
                'Normal lung fibroblasts, IL4, 24h':2,
                'Normal lung fibroblasts, unstimulated, 24h':0}
    if (tn == 5):
        atype = source
        atypes = ['Mono', 'Mac']
        ahash = {'Monocytes, unstimulated, 24h':0,
                'Monocytes, unstimulated, 6h':0,
                'Macrophages, unstimulated, 24h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChandriani2014 = getChandriani2014


def getMartinez2015(self, tn=1):
    self.prepareData("MAC89")
    atype = self.h.getSurvName('c src1')
    atypes = ['M0', 'M1', 'M2']
    ahash = { 'Monocyte-derived macrophages polarized with IL-4 for 5 days':2,
            'Monocyte-derived macrophages polarized with IL-10 for 5 days':2,
            'Monocyte-derived macrophages polarized with IFNgTNFa for 5 days':1,
            'Monocyte-derived macrophages':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinez2015 = getMartinez2015


def getDas2018(self, tn=1):
    self.prepareData("MAC90")
    atype = self.h.getSurvName('c src1')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'cultured for 4 hrs':0,
            'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 12 hrs':1,
            'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 4 hrs':1,
            'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 24 hrs':1,
            'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 1 hr':1,
            'treated with LPS (100 ng/ml) for 4 hrs':1,
            'treated with IL-13 (10\xe2\x80\x89ng/ml) for 12 hrs':2,
            'treated with IL-4 (10\xe2\x80\x89ng/ml) for 12 hrs':2}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDas2018 = getDas2018


def getDaniel2018(self, tn=1):
    self.prepareData("MAC91")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(".rep.*", "", str(k)) for k in atype]
    atype = [re.sub("mm_BMDM_", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'Wt_2nd_stim_ctrl_RNA':0,
            'Wt_1st_stim_3hIL4_RNA':2,
            'PpargKO_2nd_stim_3hIL4_RNA':2,
            'ctrl_24hVeh_RNA':0,
            'PpargKO_1st_stim_3hIL4_RNA':2,
            'PpargKO_1st_stim_ctrl_RNA':0,
            'Wt_1st_stim_ctrl_RNA':0,
            'Wt_2nd_stim_3hIL4_RNA':2,
            'PpargKO_2nd_stim_ctrl_RNA':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDaniel2018 = getDaniel2018


def getPiccolo2017(self, tn=1):
    self.prepareData("MAC92")
    atype = self.h.getSurvName('c Name')
    atype = [re.sub("_R.*", "", str(k)) for k in atype]
    atype = [re.sub("^_", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'IL4_4h':2,
            'IFNy_2h':1,
            'UT':0,
            'IFNy_IL4_2h':1,
            'IL4_2h':2,
            'IFNy_4h':1,
            'IFNy_IL4_4h':1}
    if (tn == 2):
        atype = self.h.getSurvName('c Name')
        atype = [re.sub("_R.*", "", str(k)) for k in atype]
        ahash = {'_shMYC_UT':0,
                '_scramble_IL-4_4h':2,
                '_scramble_IL-4_2h':2,
                '_scramble_UT':0,
                '_shMYC_IL-4_4h':2,
                '_shMYC_IL-4_2h':2}
        if (tn == 3):
            atype = self.h.getSurvName('c Name')
        atype = [re.sub("_R.*", "", str(k)) for k in atype]
        ahash = {'shCEBP-beta_UT':0,
                'shJunB_UT':0,
                'shJunB_IFNy_4h':1,
                'scramble_UT':0,
                'scramble_IFNy_4h':1,
                'shCEBP-beta_IFNy_4h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPiccolo2017 = getPiccolo2017


def getOstuni2013(self, tn=1):
    self.prepareData("MAC93")
    atype = self.h.getSurvName('c treatment')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'TGFb1 (1 ng/ml) for 4hrs':1,
            'IL4 (10 ng/ml) for 4hrs':2,
            'No treatment':0,
            'IFNg (100 ng/ml) for 4hrs':1,
            'TNFa (10 ng/ml) for 4hrs':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOstuni2013 = getOstuni2013


def getRochaResende(self, tn=1):
    self.prepareData("MAC94.2")
    atype = self.h.getSurvName('c treatment')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'none':0, 'LPS':1, 'IL-4':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRochaResende = getRochaResende


def getHill2018(self, tn=1):
    self.prepareData("MAC95")
    atype = self.h.getSurvName('c condition')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'Cd9+ macrophage transfer':1,
            'Ly6c+ macrophage transfer':1,
            'PBS transfer':0,
            'High fat diet':2,
            'IL4':2,
            'LPS':1,
            'Veh':0}
    if (tn == 2):
        ahash = {'PBS transfer':0,
                'IL4':2,
                'LPS':1,
                'Veh':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHill2018 = getHill2018


def getFreemerman2019(self, tn=1):
    self.prepareData("MAC96")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("\..*", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'GLUT1 WT M1':1,
            'GLUT1 KO M2':2,
            'GLUT1 WT M0':0,
            'GLUT1 WT M2':2,
            'GLUT1 KO M1':1,
            'GLUT1 KO M0':0}
    if (tn == 2):
        ahash = {'GLUT1 WT M1':1,
                'GLUT1 WT M0':0,
                'GLUT1 WT M2':2}
        if (tn == 3):
            ahash = {'GLUT1 KO M2':2,
                    'GLUT1 KO M1':1,
                    'GLUT1 KO M0':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFreemerman2019 = getFreemerman2019


def getEl2010(self, tn=1):
    self.prepareData("MAC97")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" [A-D]$", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'Irf4+/- IL4 4h':2,
            'Irf4-/- IL4 18h':2,
            'Irf4-/- Mock':0,
            'Irf4+/- IL4 18h':2,
            'Irf4+/- Mock':0,
            'Irf4-/- IL4 4h':2}
    if (tn == 2):
        ahash = {'Irf4+/- IL4 4h':2,
                'Irf4+/- IL4 18h':2,
                'Irf4+/- Mock':0}
        if (tn == 3):
            ahash = {'Irf4-/- IL4 18h':2,
                    'Irf4-/- Mock':0,
                    'Irf4-/- IL4 4h':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEl2010 = getEl2010


def getLi2015(self, tn=1):
    self.prepareData("MAC98")
    atype = self.h.getSurvName('c src1')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'Mouse macrophage at M2 treated with nutlin-3a':2,
            'Mouse macrophage at M2':2,
            'Mouse macrophage at M2 treated with 10058F4':2,
            'Mouse macrophage at M1':1,
            'Mouse macrophage at M0':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2015 = getLi2015


def getRamsey2014(self, tn=1):
    self.prepareData("MAC99")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("\..$", "", str(k)) for k in atype]
    atypes = ['Rp', 'Rs', 'Lp']
    ahash = {'LdlrKO.male.polyic':2,
            'Reversa.female.polyic':0,
            'Reversa.female.saline':1,
            'Reversa.male.saline':1,
            'Reversa.male.polyic':0}
    if (tn == 2):
        ahash = {'LdlrKO.male.polyic':2,
                'Reversa.male.saline':1,
                'Reversa.male.polyic':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRamsey2014 = getRamsey2014


def getKuo2011(self, tn=1):
    self.prepareData("MAC100")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("_m[0-9]*_", "_", str(k)) for k in atype]
    atypes = ['ML', 'MML', 'aML', 'aL']
    ahash = {'macrophage_C57BL/6-Ldlr-/-':0,
            'macrophage_C57BL/6.MOLFc4(51Mb)-Ldlr-/-':1,
            'aorta_C57BL/6.MOLFc4(51Mb)-Ldlr-/-':2,
            'aorta_C57BL/6-Ldlr-/-':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKuo2011 = getKuo2011


def getPrice2017(self, tn=1):
    self.prepareData("MAC101")
    atype = self.h.getSurvName('c genotype')
    atypes = ['DKO', 'LDLR', 'miR-33', 'WT']
    ahash = {'LDLR-/-/miR-33-/-':0, 'LDLR-/-':1, 'Wildtype':3,
            'miR-33-/-':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPrice2017 = getPrice2017


def getNicolaou2017(self, tn=1):
    self.prepareData("MAC102")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("_rep.*", "", str(k)) for k in atype]
    atype = [re.sub("mouse_", "", str(k)) for k in atype]
    atype = [re.sub("primary macrophages", "mac", str(k)) for k in atype]
    atype = [re.sub("_B6.Ldlr....Adam17", "_", str(k)) for k in atype]
    atypes = ['mwm', 'awm', 'mem', 'aem', 'mwf', 'awf', 'mef', 'aef']
    ahash = {'mac_wt/wt_male':0,
            'aorta_wt/wt_male':1,
            'mac_wt/wt_female':4,
            'aorta_wt/wt_female':5,
            'mac_ex/ex_female':6,
            'aorta_ex/ex_male':3,
            'mac_ex/ex_male':2,
            'aorta_ex/ex_female':7}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNicolaou2017 = getNicolaou2017


def getPG2019lpsRep(self, tn=1):
    #self.prepareData("MACPH1", cfile="/Users/mahdi/public_html/Hegemon/explore.conf")
    self.prepareData("MAC125")
    ttype = self.h.getSurvName("c type")
    mtype = self.h.getSurvName("c times")
    atype = [ str(ttype[i]) + " " + str(mtype[i]) for i in
            range(len(ttype))]
    atypes = ['KO 0', 'WT 0', 'WT 6hr', 'KO 6hr']
    ahash = {}
    if (tn == 2):
        atypes = ['WT 6hr', 'KO 6hr']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2019lpsRep = getPG2019lpsRep


def getPG2019lpsII(self, tn=1):
    self.prepareData("MAC125.2")
    atype = self.h.getSurvName("c Group")
    atypes = ['KO-0', 'WT-0', 'WT-6h', 'KO-6h']
    ahash = {}
    if (tn == 2):
        atypes = ['WT-6h', 'KO-6h']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2019lpsII = getPG2019lpsII


def getZhang2014(self, tn=1):
    self.prepareData("GL19")
    atype = self.h.getSurvName('c Collection time')
    atypes = ['CT22', 'CT28', 'CT34', 'CT40', 'CT46', 'CT52', 'CT58',
            'CT64']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2014 = getZhang2014


def getMure2018(self, tn=1):
    self.prepareData("GL21.2")
    tissue = self.h.getSurvName('c tissue')
    time = self.h.getSurvName('c time')
    atype = time
    atypes = ['ZT00', 'ZT02', 'ZT04', 'ZT06', 'ZT08', 'ZT10',
            'ZT12', 'ZT14', 'ZT16', 'ZT18', 'ZT20', 'ZT22']
    ahash = {}
    if tn == 2:
        t1 = "Descending colon"
        atype = [ ",".join([str(k[i]) for k in [tissue, time]]) 
                for i in range(len(atype))]
        for k in range(len(atypes)):
            ahash[",".join([t1, atypes[k]])] = k
    if tn == 3:
        t1 = "Ascending colon"
        t1 = "Retina"
        atype = [ ",".join([str(k[i]) for k in [tissue, time]]) 
                for i in range(len(atype))]
        for k in range(len(atypes)):
            ahash[",".join([t1, atypes[k]])] = k
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMure2018 = getMure2018


def getWu2018(self, tn=1):
    self.prepareData("GL22")
    atype = self.h.getSurvName('c collection time')
    atypes = ['0', '6', '12', '18']
    ahash = {'0':0, '1200':2, '600':1, '1800':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWu2018 = getWu2018


def getBraun2018(self, tn=1):
    self.prepareData("GL25")
    atype = self.h.getSurvName('c hour')
    atypes = ['1', '3', '5', '7', '9', '11', '13', '15', '17',
            '19', '21', '23', '25', '27', '29']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBraun2018 = getBraun2018


def getKervezee2019(self, tn=1):
    self.prepareData("GL26")
    atype = self.h.getSurvName('c relclocktime')
    atypes = [ '0', '4-8', '12-16', '17-18', '19-22', '24-26',
            '30-34', '36-38', '39-42']
    ahash = {'0.2':0, '0.3':0, '0.4':0, '0.5':0, '0.8':0, 
            '4.2':1, 
            '4.3':1, '6.2':1, '8.2':1, '8.3':1, '8.4':1, '8.5':1, '8.8':1,
            '12.2':2, '12.3':2, '13.2':2, '14.4':2, '16.2':2, '16.3':2,
            '17.2':3,
            '18.2':3, '18.3':3, '18.4':3, '18.9':3, '19':4, '20.2':4, '20.3':4,
            '21.8':4, '22.2':4, '22.3':4, '24.2':5, '24.3':5, '24.4':5,
            '26.2':5, '26.3':5, '26.5':5,
            '30.2':6, '30.3':6, '32.1':6, '34.2':6, '34.3':6, '34.4':6,
            '36.3':7, '38.2':7, '38.3':7, '38.8':7, '39.2':8, 
            '42.1':8, '42.2':8, '42.3':8}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKervezee2019 = getKervezee2019


def getMargerie2019(self, tn=1):
    self.prepareData("MAC103")
    atype = self.h.getSurvName('c statin')
    atypes = [ 'NS', 'R', 'A', 'S', 'P', 'F']
    ahash = {'No Statin':0, 'Rosuvastatin':1, 'Atorvastatin':2, 'Simvastatin':3,
            'Pravastatin':4, 'Fluvastatin':5}
    if (tn == 2):
        atypes = ['NS', 'R']
        ahash = {'No Statin':0, 'Rosuvastatin':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMargerie2019 = getMargerie2019


def getScicluna2015(self, tn=1):
    self.prepareData("MAC104")
    atype = self.h.getSurvName('c diabetes_mellitus')
    atypes = ['No_DM', 'NA', 'DM']
    ahash = {}
    if (tn == 2):
        atypes = ['No_DM', 'DM']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getScicluna2015 = getScicluna2015


def getKallionpaa2014I(self, tn=1):
    self.prepareData("MAC105.1")
    atype = self.h.getSurvName('c time from t1d diagnosis')
    atypes = ['N', 'E', 'L']
    ahash = {'no T1D diagnosis':0}
    for i in self.h.aRange():
        if atype[i] != 'no T1D diagnosis':
            v = float(atype[i])
            if (v >= -24 and v <= 0):
                ahash[atype[i]] = 1
            else:
                ahash[atype[i]] = 2
    if (tn == 2):
        atypes = ['N', 'T1D']
        ahash = {'no T1D diagnosis':0}
        for i in self.h.aRange():
            if atype[i] != 'no T1D diagnosis':
                ahash[atype[i]] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKallionpaa2014I = getKallionpaa2014I


def getKallionpaa2014II(self, tn=1):
    self.prepareData("MAC105.2")
    atype = self.h.getSurvName('c time from t1d diagnosis')
    atypes = ['N', 'E', 'L']
    ahash = {'no diagnosis':0}
    for i in self.h.aRange():
        if atype[i] != 'no diagnosis':
            v = float(atype[i])
            if (v >= -24 and v <= 0):
                ahash[atype[i]] = 1
            else:
                ahash[atype[i]] = 2
    if (tn == 2):
        atypes = ['N', 'T1D']
        ahash = {'no diagnosis':0}
        for i in self.h.aRange():
            if atype[i] != 'no diagnosis':
                ahash[atype[i]] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKallionpaa2014II = getKallionpaa2014II


def getKallionpaa2014III(self, tn=1):
    self.prepareData("MAC105.3")
    atype = self.h.getSurvName('c time from t1d diagnosis')
    atypes = ['N', 'E', 'L']
    ahash = {'no T1D diagnosis':0}
    for i in self.h.aRange():
        if atype[i] != 'no T1D diagnosis':
            v = float(atype[i])
            if (v >= -24 and v <= 0):
                ahash[atype[i]] = 1
            else:
                ahash[atype[i]] = 2
    if (tn == 2):
        atypes = ['N', 'T1D']
        ahash = {'no T1D diagnosis':0}
        for i in self.h.aRange():
            if atype[i] != 'no T1D diagnosis':
                ahash[atype[i]] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKallionpaa2014III = getKallionpaa2014III


def getRam2016(self, tn=1):
    self.prepareData("MAC106")
    btype = self.h.getSurvName('c type1_diabetes')
    atype = self.h.getSurvName('c condition')
    atypes = ['PMA 6h', 'Basal', 'CD8+', 'CD4+']
    ahash = {}
    if (tn == 2):
        atypes = ['1', '2']
        atype = btype
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRam2016 = getRam2016


def getAlmon2009(self, tn=1):
    self.prepareData("MAC107")
    atype = self.h.getSurvName("c desc")
    ttype = [str(k).split(" ")[5] if len(str(k).split(" ")) > 5 else None
            for k in atype]
    atype = [re.sub(".*feeded with ", "", str(k)) for k in atype]
    mtype = [re.sub("and .*", "", str(k)) for k in atype]
    atype = [str(ttype[i]) + " " + mtype[i] for i in range(len(atype))]
    atypes = ['aN', 'aH', 'lN', 'lH', 'mN', 'mH']
    ahash = {'adipose ND ':0,
            'adipose HFD ':1,
            'liver ND ':2,
            'liver HFD ':3,
            'muscle ND ':4,
            'muscle HFD ':5,
            'livadipose ND ':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAlmon2009 = getAlmon2009


def getChen2014(self, tn=1):
    self.prepareData("MAC109")
    atype = self.h.getSurvName("c stimulated with")
    atype = [re.sub("auto.*\(", "(", str(k)) for k in atype]
    atype = [re.sub(" HLA risk sibling", "", str(k)) for k in atype]
    atype = [re.sub(" plasma", "", str(k)) for k in atype]
    atype = [re.sub(" series", "", str(k)) for k in atype]
    atype = [re.sub("Longitudinal ", "", str(k)) for k in atype]
    atypes = ['nl', 'nh', 'phNp', 'o', 'nhNp', 'nlNp', 'P']
    ahash = {'(AA-) low':0,
            '(AA-) high':1,
            '(AA+) high non-progressor':2,
            'recent onset cultured with IL1RA':3,
            '(AA-) high non-progressor':4,
            '(AA-) low non-progressor':5,
            'progressor':6}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChen2014 = getChen2014


def getLefebvre2017(self, tn=1):
    self.prepareData("LIV3")
    atype = self.h.getSurvName("c Title")
    atype = [str(k).replace("liver biopsy ", "") for k in atype]
    self.patient = [str(k).split(" ")[0] for k in atype]
    atype = [str(k).split(" ")[1] if len(str(k).split(" ")) > 1 else "" for
            k in atype]
    self.paired = [re.sub("\((.*)\)", "\\1", str(k)) for k in atype]
    self.time = self.h.getSurvName("c time")
    self.treatment = self.h.getSurvName("c type of intervention")
    atype = self.h.getSurvName("c src1")
    atypes = ['b', 'f', 'Nb', 'Nf', 'ub', 'uf']
    ahash = {'no NASH liver baseline':0,
            'no NASH liver follow-up':1,
            'NASH liver follow-up':3,
            'NASH liver baseline':2,
            'undefined liver baseline':4,
            'undefined liver follow-up':5}
    self.aval = [ahash[i] if i in ahash else None for i in atype]
    phash = {}
    for i in self.h.aRange():
        phash[self.patient[i]] = i
    self.rtype = [None, None] + [str(self.aval[phash[self.paired[i]]])+" "+\
            str(self.aval[i])+" " + str(self.treatment[i]) \
            if self.paired[i] in phash and self.paired[i] != "" \
            else str(self.aval[i])+" " + str(self.treatment[i]) \
            for i in self.h.aRange()]
    fhash = {}
    for i in self.h.aRange():
        if self.paired[i] in phash and self.paired[i] != "":
            fhash[self.paired[i]] = i
    self.ftype = [None, None] + [str(self.aval[i])+" "+str(self.treatment[i])\
            + " " + str(self.aval[fhash[self.patient[i]]])\
            if self.patient[i] in fhash\
            else str(self.aval[i])+" " + str(self.treatment[i]) \
            for i in self.h.aRange()]
    if (tn == 2):
        atypes = ['b', 'Nb', 'ub']
        ahash = {'no NASH liver baseline':0,
                'NASH liver baseline':1,
                'undefined liver baseline':2}
    if (tn == 3):
        atypes = ['b', 'Nb']
        ahash = {'no NASH liver baseline':0,
                'NASH liver baseline':1}
    if (tn == 4):
        atypes = ['f', 'Nf']
        ahash = {'no NASH liver follow-up':0,
                'NASH liver follow-up':1}
    if (tn == 5):
        atypes = ['R', 'NR']
        atype = self.ftype
        ahash = {'2 Diet 1': 0, '2 Diet 3': 1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLefebvre2017 = getLefebvre2017


def getGallagher2010(self, tn=1):
    self.prepareData("MAC110")
    atype = self.h.getSurvName("c Type")
    atypes = ['N', 'GI', 'D']
    ahash = {'diabetic':2, 'glucoseIntolerant':1, 'normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGallagher2010 = getGallagher2010


def getDuPlessis2015(self, tn=1):
    self.prepareData("MAC111")
    atype = self.h.getSurvName("c src1")
    atypes = ['S1', 'V1', 'S2', 'V2', 'S3', 'V3', 'S4', 'V4']
    ahash = {'Subc Fat, Histology class 1':0,
            'Visceral Fat, Histology class 1':1,
            'Subc Fat, Histology class 2':2,
            'Visceral Fat, Histology class 2':3,
            'Subc Fat, Histology class 3':4,
            'Visceral Fat, Histology class 3':5,
            'Subc Fat, Histology class 4':6,
            'Visceral Fat, Histology class 4':7}
    if (tn == 2):
        atypes = ['N', 'F']
        ahash = {'Subc Fat, Histology class 1':0,
                'Visceral Fat, Histology class 1':0,
                'Subc Fat, Histology class 2':0,
                'Visceral Fat, Histology class 2':0,
                'Subc Fat, Histology class 3':0,
                'Visceral Fat, Histology class 3':0,
                'Subc Fat, Histology class 4':1,
                'Visceral Fat, Histology class 4':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuPlessis2015 = getDuPlessis2015


def getWu2007IS(self, tn=1):
    self.prepareData("MAC112")
    atype = self.h.getSurvName("c status")
    atypes = ['IS', 'IR', 'D']
    ahash = {'insulin sensitive':0, 'diabetic':2, 'insulin resistant':1}
    if (tn == 2):
        atypes = ['IS', 'IR']
        ahash = {'insulin sensitive':0, 'insulin resistant':1}
    if (tn == 3):
        atypes = ['IS', 'D']
        ahash = {'insulin sensitive':0, 'diabetic':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWu2007IS = getWu2007IS


def getStenvers2019(self, tn=1):
    self.prepareData("MAC113")
    atype = self.h.getSurvName("c subject status")
    atypes = ['H', 'T2D']
    ahash = {'Type 2 diabetes':1, 'Healthy':0}
    if (tn >= 2):
        stype = self.h.getSurvName("c subject status")
        ttype = self.h.getSurvName("c timepoint")
        atype = [ str(stype[i]) + " " + str(ttype[i]) for i in
                range(len(stype))]
        atypes = ['H1', 'H2', 'H3', 'H4', 'D1', 'D2', 'D3', 'D4']
        ahash = {'Type 2 diabetes D2_ZT_15:30':4,
                'Type 2 diabetes D3_ZT_0:15':5,
                'Type 2 diabetes D3_ZT_5:45':6,
                'Type 2 diabetes D3_ZT_11:15':7,
                'Healthy D2_ZT_15:30':0,
                'Healthy D3_ZT_0:15':1,
                'Healthy D3_ZT_5:45':2,
                'Healthy D3_ZT_11:15':3}
    if (tn == 3):
        atypes = ['H1', 'H2', 'H3', 'H4']
        ahash = {'Healthy D2_ZT_15:30':0,
                'Healthy D3_ZT_0:15':1,
                'Healthy D3_ZT_5:45':2,
                'Healthy D3_ZT_11:15':3}
    if (tn == 4):
        atypes = ['D1', 'D2', 'D3', 'D4']
        ahash = {'Type 2 diabetes D2_ZT_15:30':0,
                'Type 2 diabetes D3_ZT_0:15':1,
                'Type 2 diabetes D3_ZT_5:45':2,
                'Type 2 diabetes D3_ZT_11:15':3}
    if (tn == 5):
        atypes = ['H1', 'H2', 'D1', 'D2']
        ahash = {'Type 2 diabetes D2_ZT_15:30':2,
                'Type 2 diabetes D3_ZT_0:15':3,
                'Healthy D2_ZT_15:30':0,
                'Healthy D3_ZT_0:15':1}
    if (tn == 6):
        atypes = ['H', 'D']
        ahash = {'Type 2 diabetes D2_ZT_15:30':1,
                'Type 2 diabetes D3_ZT_0:15':1,
                'Healthy D2_ZT_15:30':0,
                'Healthy D3_ZT_0:15':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getStenvers2019 = getStenvers2019


def getDAmore2018(self, tn=1):
    self.prepareData("MAC114")
    atype = self.h.getSurvName("c disease state")
    atypes = ['C', 'MetS']
    ahash = {'control':0, 'Metabolic Syndrome':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDAmore2018 = getDAmore2018


def getArnardottir2014(self, tn=1):
    self.prepareData("GL27")
    rtype = self.h.getSurvName("c responder")
    htype = self.h.getSurvName("c hour")
    btype = self.h.getSurvName("c biological group")
    atype = [ str(btype[i]) + " " + str(rtype[i]) + " " + str(htype[i])
                     for i in range(len(rtype))]
    atypes = ['bl0', 'bl4', 'bl8', 'bl12', 'bl16', 'bl20',
            'bh0', 'bh4', 'bh8', 'bh12', 'bh16', 'bh20',
            'sl0', 'sl4', 'sl8', 'sl12', 'sl16', 'sl20',
            'rl0', 'rl4', 'rl8',
            'sh0', 'sh4', 'sh8', 'sh12', 'sh16', 'sh20',
            'rh0', 'rh4', 'rh8']
    ahash = {
            'baseline low 0':0,
            'baseline low 4':1,
            'baseline low 8':2,
            'baseline low 12':3,
            'baseline low 16':4,
            'baseline low 20':5,
            'baseline high 0':6,
            'baseline high 4':7,
            'baseline high 8':8,
            'baseline high 12':9,
            'baseline high 16':10,
            'baseline high 20':11,
            'sleep deprivation low 0':12,
            'sleep deprivation low 4':13,
            'sleep deprivation low 8':14,
            'sleep deprivation low 12':15,
            'sleep deprivation low 16':16,
            'sleep deprivation low 20':17,
            'recovery low 0':18,
            'recovery low 4':19,
            'recovery low 8':20,
            'sleep deprivation high 0':21,
            'sleep deprivation high 4':22,
            'sleep deprivation high 8':23,
            'sleep deprivation high 12':24,
            'sleep deprivation high 16':25,
            'sleep deprivation high 20':26,
            'recovery high 0':27,
            'recovery high 4':28,
            'recovery high 8':29}
    if (tn == 2):
        atype = [ str(btype[i]) + " " + str(rtype[i])
                for i in range(len(rtype))]
        atypes = ['bl', 'bh', 'sl', 'rl', 'sh', 'rh']
        ahash = {'baseline low':0,
                'baseline high':1,
                'sleep deprivation low':2,
                'recovery low':3,
                'sleep deprivation high':4,
                'recovery high':5}
    if (tn == 3):
        atype = rtype
        atypes = ['high', 'low']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArnardottir2014 = getArnardottir2014


def getADPooledDyn(self, tn=1):
    self.prepareData("AD7")
    atype = self.h.getSurvName('c AD specific');
    ahash = {'0':0, '1':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'AD']
    ahash = {'Normal':0, "Alzheimer's Disease":1, 'normal':0,
            'definite AD':0, 'Control':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getADPooledDyn = getADPooledDyn


def getLiang2007(self):
    self.prepareData("AD2")
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'AD']
    ahash = {'normal\xa0':0, "Alzheimer's Disease\xa0":1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiang2007 = getLiang2007


def getFriedman2017(self):
    self.prepareData("AD3")
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['N', 'AD']
    ahash = {'control':0, "Alzheimer's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFriedman2017 = getFriedman2017


def getBerchtold2018(self, tn = 1):
    self.prepareData("AZ3", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c physical activity tier');
    atypes = ['H', 'M', 'L']
    ahash = {'high activity':0,
            'low activity':2,
            'moderate activity':1}
    if tn == 2:
        atypes = ['H', 'L']
        ahash = {'high activity':0,
                'low activity':1,
                'moderate activity':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2018 = getBerchtold2018


def getSimpson2015(self, tn = 1):
    self.prepareData("AZ5", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c neuronal ddr');
    atypes = ['H', 'L']
    ahash = {'High DNA damage':0, 'Low DNA damage':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSimpson2015 = getSimpson2015


def getPiras2019(self, tn = 1):
    self.prepareData("AZ1", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    age = self.h.getSurvName('c expired_age');
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['N', 'AD']
    ahash = {'ND':0, 'AD':1}
    if tn == 2:
        atype = [atype[i] if int(age[i]) > 90
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPiras2019 = getPiras2019


def getPatel2019(self, tn = 1):
    self.prepareData("AD8")
    atype = self.h.getSurvName('c tissue');
    ahash = {'Temporal_Cortex':3,
            'Cerebellum':4,
            'Frontal_Cortex':5,
            'Entorhinal_Cortex':6}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state');
    atypes = ['N', 'A', 'AD']
    ahash = {'AsymAD':1, 'AD':2, 'control':0}
    if (tn >= 2):
        atypes = ['N', 'AD']
        ahash = {'AD':1, 'control':0}
    if (tn > 2):
        atype = [atype[i] if rval[i] == tn else None for i in range(len(atype))]
    self.rval = rval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPatel2019 = getPatel2019


def getNarayanan2014(self, tn = 1):
    self.prepareData("AZ8", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c disease status');
    atypes = ['N', 'AD', 'HD']
    ahash = {'non-demented':0,
            "Alzheimer's disease":1,
            "Huntington's disease":2}
    if tn == 2:
        atypes = ['N', 'AD']
        ahash = {'non-demented':0,
                "Alzheimer's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNarayanan2014 = getNarayanan2014


def getZhang2013m(self, dbid = 'AZ12', tn = 1):
    self.db = hu.Database("/Users/mgosztyl/public_html/Hegemon/explore.conf")
    self.dbid = dbid
    atype = self.h.getSurvName('c disease');
    atypes = ['N', 'A']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2013m = getZhang2013m


def getBerchtold2014(self, tn = 1):
    self.prepareData("AZ11", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c src1")
    atype = [str(i) for i in atype]
    res = []
    for k in atype:
        l1 = k.split(",")
        if (len(l1) != 4):
            res.extend([k])
        else:
            res.extend([l1[1].strip() + "_" + l1[2].strip()])
    atype = res
    atypes = ['N', 'AD']
    ahash = {'entorhinal cortex_male':0,
            'entorhinal cortex_male_AD':1,
            'entorhinal cortex_female':0,
            'entorhinal cortex_female_AD':1,
            'superior frontal gyrus_male':0,
            'superior frontal gyrus_male_AD':1,
            'superior frontal gyrus_female':0,
            'superior frontal gyrus_female_AD':1,
            'postcentral gyrus_male':0,
            'post-central gyrus_male_AD':1,
            'postcentral gyrus_female':0,
            'post-central gyrus_female_AD':1,
            'hippocampus_male':0,
            'hippocampus_male_AD':1,
            'hippocampus_female':0,
            'hippocampus_female_AD':1}
    if (tn == 2):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'hippocampus_male':0,
                'hippocampus_male_AD':1}
    if (tn == 3):
        ahash = {'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    if (tn == 4):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1}
    if (tn == 5):
        ahash = {'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1}
    if (tn == 6):
        ahash = {'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1}
    if (tn == 7):
        ahash = {'hippocampus_male':0,
                'hippocampus_male_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2014 = getBerchtold2014


def getWang2016(self, dbid="AD9", tn = 1):
    self.dbid = dbid
    atype = self.h.getSurvName('c brain region')
    ahash = {'Inferior Temporal Gyrus':3,
            'Parahippocampal Gyrus':4,
            'Middle Temporal Gyrus':5,
            'Occipital Visual Cortex':6,
            'Prefrontal Cortex':7,
            'Hippocampus':8,
            'Caudate Nucleus':9,
            'Frontal Pole':10,
            'Precentral Gyrus':11,
            'Posterior Cingulate Cortex':12,
            'Superior Temporal Gyrus':13,
            'Superior Parietal Lobule':14,
            'Temporal Pole':15,
            'Anterior Cingulate':16,
            'Inferior Frontal Gyrus':17,
            'Dorsolateral Prefrontal Cortex':18,
            'Putamen':19}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c neuropathological category')
    atypes = ['N', 'AD']
    ahash = {'definite AD':1, 'Possible AD':1,
            'Normal':0, 'Probable AD':1}
    if (tn >= 2):
        atypes = ['N', 'AD']
        ahash = {'definite AD':1, 'Normal':0}
    if (tn >= 3):
        atype = [atype[i] if rval[i] == tn else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2016 = getWang2016


def getBerchtold2014RMA(self, tn=1):
    self.prepareData("AD11")
    atype = self.h.getSurvName('c AD specific');
    ahash = {'0':0, '1':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c src1")
    atype = [str(i) for i in atype]
    res = []
    for k in atype:
        l1 = k.split(",")
        if (len(l1) != 4):
            res.extend([k])
        else:
            res.extend([l1[1].strip() + "_" + l1[2].strip()])
    atype = res
    atypes = ['N', 'AD']
    ahash = {'entorhinal cortex_male':0,
            'entorhinal cortex_male_AD':1,
            'entorhinal cortex_female':0,
            'entorhinal cortex_female_AD':1,
            'superior frontal gyrus_male':0,
            'superior frontal gyrus_male_AD':1,
            'superior frontal gyrus_female':0,
            'superior frontal gyrus_female_AD':1,
            'postcentral gyrus_male':0,
            'post-central gyrus_male_AD':1,
            'postcentral gyrus_female':0,
            'post-central gyrus_female_AD':1,
            'hippocampus_male':0,
            'hippocampus_male_AD':1,
            'hippocampus_female':0,
            'hippocampus_female_AD':1}
    if (tn == 4):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'hippocampus_male':0,
                'hippocampus_male_AD':1}
    if (tn == 5):
        ahash = {'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    if (tn == 6):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1}
    if (tn == 7):
        ahash = {'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1}
    if (tn == 8):
        ahash = {'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1}
    if (tn == 9):
        ahash = {'hippocampus_male':0,
                'hippocampus_male_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2014RMA = getBerchtold2014RMA


def getGelman2012(self, tn = 1):
    self.prepareData("DE39", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c phenotype');
    atypes = ['N', 'HAND']
    ahash = {'HIV Infected':0,
            'HIV Infected with neurocognitive impairment (HAD: HIV-associated dementia)':1,
            'HIV Infected with HAD and HIV encephalitis (HIVE)':1,
            'normal (control)':0,
            'HIV Infected with HAD':1,
            'HIV Infected with HAD and encephalitis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGelman2012 = getGelman2012

    
def getChenPlotkin2008(self, tn = 1):
    self.prepareData("DE1", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c src1")
    atype = [str(k).split(" ")[1] if len(str(k).split(" ")) > 1 else None
            for k in atype]
    atype = [str(k).split("-")[0] for k in atype]
    atypes = ['N', 'P', 'S']
    ahash = {'Normal':0,
            'Progranulin':1,
            'Sporadic':2}
    if tn == 2:
        atypes = ['N', 'FTD']
        ahash = {'Normal':0,
                'Progranulin':1}
    if tn == 3:
        atype = self.h.getSurvName('c disease and tissue')
        atypes = ['N', 'FTD']
        ahash = {'Normal hippocampus':0,
            'Progranulin hippocampus':1,
            'Sporadic hippocampus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChenPlotkin2008 = getChenPlotkin2008

    
def getOlmosSerrano2016(self, tn = 1):
    self.prepareData("DE49", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c disease status');
    atypes = ['N', 'DS']
    ahash = {'CTL':0,
            'DS':1}
    if tn == 2:
        atype = self.h.getSurvName('c disease and tissue');
        atypes = ['N', 'DS']
        ahash = {'CTL ITC':0,
            'CTL STC':0,
            'CTL HIP':0,
            'DS ITC':1, 'DS STC':1, 'DS HIP':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOlmosSerrano2016 = getOlmosSerrano2016

    
def getWilliams2009(self, tn = 1):
    self.prepareData("DE51", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'MCI']
    ahash = {'Control':0, 'MCI':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWilliams2009 = getWilliams2009

    
def getBartolettiStella2019 (self, tn = 1):
    self.prepareData("DE52", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c condition');
    atypes = ['N', 'CJD']
    ahash = {'Control':0, 'sCJD affected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBartolettiStella2019  = getBartolettiStella2019 

    
def getWes2014a (self, tn = 1):
    self.prepareData("AZ78", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', 'rTg4510']
    ahash = {'Dbl Neg':0, 'Tg4510':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWes2014a  = getWes2014a 


def getWes2014b (self, tn = 1):
    self.prepareData("AZ77", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', 'rTg4510']
    ahash = {'Dbl Neg':0, 'Tg4510':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWes2014b  = getWes2014b 

    
def getWes2014c (self, tn = 1):
    self.prepareData("AZ78", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', 'rTg4510']
    ahash = {'Dbl Neg':0, 'Tg4510':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWes2014c  = getWes2014c 

    
def getHokama2013 (self, tn = 1):
    self.prepareData("AZ49", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', '3xTg']
    ahash = {'non-Tg':0, '3xTg-Homo':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHokama2013  = getHokama2013 

    
def getWakutani2014 (self, tn = 1):
    self.prepareData("AZ60", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genetic background');
    atypes = ['N', 'TgCRND8']
    ahash = {'non-transgenic littermate mouse':0,
            'TgCRND8 transgenic mouse':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWakutani2014  = getWakutani2014 

    
def getMeyer2019 (self, dbid = "AZ29", tn = 1):
    self.db = hu.Database("/Users/mgosztyl/public_html/Hegemon/explore.conf")
    self.dbid = dbid
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['N', 'AD']
    ahash = {'normal':0,
            'sporadic Alzheimer\'s disease':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeyer2019  = getMeyer2019 

    
def getScheckel2019 (self, tn = 1):
    self.prepareData("AZ95", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c source_name (ch1)');
    atypes = ['N', 'AD']
    ahash = {'wiltype iPSC-derived neurons':0,
            'APP/PSEN1 double mutant iPSC-derived neurons':1}
    if tn == 2:
        ahash = {'wiltype iPSC-derived neurons':0,
            'APP mutant iPSC-derived neurons':1}
    if tn == 3:
        ahash = {'wiltype iPSC-derived neurons':0,
            'PSEN1 mutant iPSC-derived neurons':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getScheckel2019  = getScheckel2019 


def getTsalik2015(self, tn=1):
    self.prepareData("ms32.0", "/Users/sataheri/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c sirs outcomes (ch1)")
    atypes = ['SHK', 'SS', 'SIRS', 'US', 'SD']
    ahash = {'Septic shock':0,
            'severe sepsis':1,
            'SIRS':2,
            'Uncomplicated sepsis':3,
            'sepsis death':4}
    if (tn == 2):
        atype = self.h.getSurvName("c sirs vs sepsis (ch1)")
        atypes = ['Sepsis', 'SIRS']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTsalik2015 = getTsalik2015


def getBarcella2018I(self, tn=1):
    self.prepareData("ms33.1", "/Users/sataheri/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c clinical classification (ch1)")
    atypes = ['R', 'NR']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBarcella2018I = getBarcella2018I


def getBarcella2018II(self, tn=1):
    self.prepareData("ms33.2", "/Users/sataheri/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c clinical classification (ch1)")
    atypes = ['R', 'NR']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBarcella2018II = getBarcella2018II


def getBarcella2018(self, tn=1):
    self.prepareData("MAC115.1")
    ctype = self.h.getSurvName("c clinical classification")
    ttype = self.h.getSurvName("c timepoint")
    atype = [ str(ctype[i]) + " " + str(ttype[i]) for i in
                            range(len(ctype))]
    atypes = ['R T1', 'R T2', 'NR T1', 'NR T2']
    ahash = {}
    if (tn == 2):
        atypes = ['R', 'NR']
        ahash = {'R T2':0, 'NR T1':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBarcella2018 = getBarcella2018


def getDolinay2012(self, tn=1):
    self.prepareData("MAC116")
    atype = self.h.getSurvName("c src1")
    atypes = ['U', 'RS0', 'S0', 'S7', 'SA0', 'SA7']
    ahash = {'SIRS Day 0':1, 'Sepsis Day 7':3, 'se/ARDS Day 7':5,
            'se/ARDS Day 0':4, 'untreated':0, 'Sepsis Day 0':2}
    if (tn == 2):
        atypes = ['U', 'S']
        ahash = {'SIRS Day 0':1, 'Sepsis Day 7':1, 'se/ARDS Day 7':1,
                'se/ARDS Day 0':1, 'untreated':0, 'Sepsis Day 0':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDolinay2012 = getDolinay2012


def getJuss2016(self, tn=1):
    self.prepareData("MAC117")
    atype = self.h.getSurvName("c src1")
    atypes = ['HVT', 'ARDS', 'C0', 'C6', 'pan6', 'E6', 'd6', 'g6', 'p6']
    ahash = {'PMNs from HVT':0,
            'PMNs from ARDS patient':1,
            'Cultured PMNs from HVT, DMSO, t=0':2,
            'Cultured PMNs from HVT, DMSO, t=6':3,
            'Cultured PMNs from HVT + panPI3K Inhibitor, t=6h':4,
            'Cultured PMNs from HVT + GM-CSF, t=6h':5,
            'Cultured PMNs from HVT + GM-CSF + PI3Kd Inhibitor, t=6h':6,
            'Cultured PMNs from HVT + GM-CSF + PI3Kg Inhibitor, t=6h':7,
            'Cultured PMNs from HVT + GM-CSF + PanPI3K Inhibitor, t=6h':8}
    if (tn == 2):
        atypes = ['HVT', 'ARDS']
        ahash = {'PMNs from HVT':0,
                'PMNs from ARDS patient':1}
    if (tn == 3):
        atypes = ['C', 'T']
        ahash = {'Cultured PMNs from HVT, DMSO, t=0':0,
                'Cultured PMNs from HVT + GM-CSF, t=6h':1}
                
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJuss2016 = getJuss2016


def getKangelaris2015(self, tn=1):
    self.prepareData("MAC118")
    atype = self.h.getSurvName("c disease state")
    atypes = ['S', 'SA']
    ahash = {'sepsis alone':0, 'sepsis with ARDS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKangelaris2015 = getKangelaris2015


def getLu2014(self, tn = 1):
    self.prepareData("AZ21", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c age")
    atype = [str(k).split(" ")[0] for k in atype]
    atypes = ['Y', 'O']
    ahash = {}
    for i in self.h.aRange():
        if float(atype[i]) > 60:
            ahash[atype[i]] = 1
        else:
            ahash[atype[i]] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLu2014 = getLu2014


def getMarttinen2019(self, tn=1):
    self.prepareData("AD12")
    atype = self.h.getSurvName("c braak stage")
    atypes = ['0', '1', '2', '3', '4', '5', '6']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c age")
        atypes = ['Y', 'O']
        for i in self.h.aRange():
            if float(atype[i]) > 70:
                ahash[atype[i]] = 1
            else:
                ahash[atype[i]] = 0
    if (tn == 3):
        atypes = ['L', 'H']
        ahash = {'0':0, '2':0, '3':1, '4':1, '5':1, '6':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMarttinen2019 = getMarttinen2019


def getDonega2019(self, tn=1):
    self.prepareData("AD13.2")
    atype = self.h.getSurvName("c Type")
    ahash = {'SVZ':2, 'CD271':3, 'CD11b':4}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Disease")
    atypes = ['C', 'PD']
    ahash = {'Cntr':0, 'PD':1}
    if (tn > 1):
        atype = [atype[i] if rval[i] == tn
                else None for i in range(len(atype))]
    self.rval = rval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDonega2019 = getDonega2019


def getWatson2017(self, tn=1):
    self.prepareData("MAC120")
    atype = self.h.getSurvName("c sleep duration")
    atypes = ['long', 'short']
    if (tn == 2):
        atypes = ['short', 'long']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWatson2017 = getWatson2017


def getUyhelji2018(self, tn=1):
    self.prepareData("MAC119")
    atype = self.h.getSurvName("c subject group")
    atypes = ['Sleep Deprived', 'Control']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getUyhelji2018 = getUyhelji2018


def getMaret2007(self, tn=1):
    self.prepareData("MAC121")
    atype = self.h.getSurvName("c src1")
    strain = [re.split(", *", str(i))[0] for i in atype]
    expt = [re.split(", *", str(i))[1] if len(re.split(", *", str(i))) > 1
            else None for i in atype]
    time = [re.split(", *", str(i))[2] if len(re.split(", *", str(i))) > 2
            else None for i in atype]
    rep = [re.split(", *", str(i))[3] if len(re.split(", *", str(i))) > 3
            else None for i in atype]
    tissue = [re.split(", *", str(i))[4] if len(re.split(", *", str(i))) >
            4 else None for i in atype]
    atype = expt
    atypes = ['C', 'SD']
    ahash = {'sleep deprived':1,
            'control':0,
            '6 hrs sleep deprivation':1,
            '6hrs sleep deprivation':1}
    if (tn == 2):
        atype = [ ",".join([str(k[i]) for k in [strain, expt, time,
            tissue]]) for i in range(len(atype))]
        atypes = ['C', 'SD']
        ahash = {'C57BL/6J,sleep deprived,time of sacrifice ZT 0,None':1,
                'C57BL/6J,control,time of sacrifice ZT 0,None':0}
    if (tn == 3):
        atype = [ ",".join([str(k[i]) for k in [strain, expt, time,
            tissue]]) for i in range(len(atype))]
        atypes = ['C', 'SD']
        ahash = {'AKR/J,control,time of sacrifice ZT 0,None':0,
                'AKR/J,sleep deprived,time of sacrifice ZT 0,None':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaret2007 = getMaret2007


def getLaing2017(self, tn=1):
    self.prepareData("MAC122")
    time = self.h.getSurvName("c time sample taken")
    atype = self.h.getSurvName("c sleep protocol")
    atype = [ ",".join([str(k[i]) for k in [atype, time]])
                     for i in range(len(atype))]
    atypes = ['R1', 'R2', 'R3', 'R3', 'E1', 'E2', 'E3']
    ahash = {'Sleep Restriction,1':0,
            'Sleep Restriction,6':1,
            'Sleep Restriction,7':2,
            'Sleep Restriction,10':3,
            'Sleep Extension,1':4,
            'Sleep Extension,6':5,
            'Sleep Extension,10':6}
    if (tn == 2):
        atypes = ['R1', 'R2', 'R3', 'E1', 'E2', 'E3']
        ahash = {'Sleep Restriction,1':0,
                'Sleep Restriction,6':1,
                'Sleep Restriction,10':2,
                'Sleep Extension,1':3,
                'Sleep Extension,6':4,
                'Sleep Extension,10':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLaing2017 = getLaing2017


def getResuehr2019(self, tn=1):
    self.prepareData("MAC123")
    atype = self.h.getSurvName("c work shift")
    atypes = ['D', 'N']
    ahash = {'Day-Shift':0, 'Night-Shift':1}
    if (tn == 2):
        shift = self.h.getSurvName("c work shift")
        time = self.h.getSurvName("c time of sample")
        atype = [time[i] if shift[i] == 'Day-Shift' else None
                for i in range(len(time))]
        atypes = ['9', '12', '15', '18', '21', '24', '27', '30']
        ahash = {}
    if (tn == 3):
        shift = self.h.getSurvName("c work shift")
        time = self.h.getSurvName("c time of sample")
        atype = [time[i] if shift[i] == 'Night-Shift' else None
                for i in range(len(time))]
        atypes = ['9', '12', '15', '18', '21', '24', '27', '30']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getResuehr2019 = getResuehr2019


def getKervezee2018(self, tn=1):
    self.prepareData("GL28")
    atype = self.h.getSurvName('c condition')
    atypes = ['D', 'N']
    ahash = {'NightOrientedSchedule':1, 'DayOrientedSchedule':0}
    if (tn == 2):
        atype = self.h.getSurvName("c relclocktime")
        atypes = ['D1', 'D2', 'N1', 'N2', 'N3']
        v1 = [6, 12, 18, 24, 48]
        ahash = {}
        for i in self.h.aRange():
            t = atype[i]
            if t is None:
                continue
            ahash[t] = np.searchsorted(v1, float(t), 'left')
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKervezee2018 = getKervezee2018


def getChristou2019(self, tn=1):
    self.prepareData("GL29")
    atype = self.h.getSurvName('c circadian phase')
    atypes = ['D1', 'D2', 'D3', 'N1', 'N2', 'N3']
    v1 = [-4, 0, 6, 12, 18]
    ahash = {}
    for i in self.h.aRange():
        t = atype[i]
        if t is None:
            continue
        ahash[t] = np.searchsorted(v1, float(t), 'left')
    if (tn == 2):
        atypes = ['N', 'D']
        atype = [ahash[i] if i in ahash else None for i in atype]
        ahash = { 1: 0, 2: 0, 3:1 }
    if (tn == 3):
        atypes = ['D', 'N']
        atype = [ahash[i] if i in ahash else None for i in atype]
        ahash = { 1: 1, 2: 1, 3:0 }
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChristou2019 = getChristou2019


def getZieker2010(self, tn=1):
    self.prepareData("GL30")
    atype = self.h.getSurvName('c time point')
    atypes = ['0h', '6h', '12h', '18h']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZieker2010 = getZieker2010


def getSteidl2010(self, tn=1):
    self.prepareData("LM7")
    atype = self.h.getSurvName('c Title')
    atype = [ re.split(",", str(i))[0] for i in atype]
    atypes = ['class S', 'class F']
    if (tn == 2):
        atype = self.h.getSurvName('c disease stage')
        atypes = ['1', '2', '3', '4']
    if (tn == 3):
        stage = self.h.getSurvName('c disease stage')
        atype = self.h.getSurvName('c Title')
        atype = [ re.split(",", str(atype[i]))[0] if stage[i] == '4' \
                else None for i in range(len(atype))]
        atypes = ['class S', 'class F']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSteidl2010 = getSteidl2010


def getVicente2012(self, tn=1):
    self.prepareData("LM10")
    atype = self.h.getSurvName('c disease status')
    atypes = ['HC', 'FCL', 'SMZL', 'ML', 'DLBCL']
    ahash = {'MALT lymphoma':3,
            'DLBCL':4,
            'SMZL':2,
            'FCL':1,
            'healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVicente2012 = getVicente2012


def getBrech2020(self, tn=1):
    self.prepareData("MACV1")
    atype = self.h.getSurvName('c cell type')
    atypes = ['rccM', 'M1', 'M2', 'Mo', 'rcceD', 'Dc']
    ahash = {'macrophages':0,
            'M1 macrophages':1,
            'M2 macrophages':2,
            'Monocytes':3,
            'ercDCs':4,
            'CD1c+ DCs':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBrech2020 = getBrech2020


def getMontoya2009(self, tn=1):
    self.prepareData("MACV2")
    atype = self.h.getSurvName('c src1')
    atypes = ['BT', 'LL', 'RR']
    ahash = {'lepromatous leprosy':1,
            'borderline tuberculoid leprosy':0,
            'reversal reaction leprosy':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMontoya2009 = getMontoya2009


def getBurckstummer2009(self, tn=1):
    self.prepareData("MACV3")
    atype = self.h.getSurvName('c Title')
    atypes = ['U', 'S']
    ahash = {'L929 unstimulated':0,
            'NIH3T3 unstimulated':0,
            'RAW 264.7 unstimulated':0,
            'L929 + IFN-beta (4hrs)':1,
            'NIH3T3 + IFN-beta (4hr)':1,
            'RAW264.7 IFN-beta (4hrs)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurckstummer2009 = getBurckstummer2009


def getMan2015(self, tn=1):
    self.prepareData("MACV4")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Irf1', 'Aim2', 'Ifnar1']
    ahash = {'wild-type':0, 'Irf1-/-':1, 'Aim2-/-':2, 'Ifnar1-/-':3}
    if (tn == 2):
        atypes = ['Wt', '', 'Irf1']
        ahash = {'wild-type':0, 'Irf1-/-':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMan2015 = getMan2015


def getGray2016Hs(self, tn=1):
    self.prepareData("MACV5")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', '4', '12']
    ahash = {'AP1 dimerizer drug (12h)':2,
            'AP1 dimerizer drug (4h)':1,
            'Mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGray2016Hs = getGray2016Hs


def getGray2016Mm(self, tn=1):
    self.prepareData("MACV6")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'T']
    ahash = {'CT-DNA':1, 'mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGray2016Mm = getGray2016Mm


def getKarki2018(self, tn=1):
    self.prepareData("MACV7")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Irf8']
    ahash = {'Wild-type':0, 'IRF8-/-':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKarki2018 = getKarki2018


def getPan2017(self, tn=1):
    self.prepareData("MACV8")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Tet2']
    ahash = {'wildtype':0, 'Tet2 knockout':1}
    if (tn == 2):
        ctype = self.h.getSurvName('c cell type')
        ahash = {'bone marrow derived macrophages (BMDMs)':0,
                'tumor associated macrophages (TAMs)':1}
        rval = [ahash[i] if i in ahash else None for i in ctype]
        atype = self.h.getSurvName('c genotype/variation')
        atype = [atype[i] if rval[i] == 1 else None \
                for i in range(len(atype))]
        atypes = ['Wt', 'Tet2']
        ahash = {'wildtype':0, 'Tet2 knockout':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPan2017 = getPan2017


def getJura2008(self, tn=1):
    self.prepareData("MACV9")
    atype = self.h.getSurvName('c Title')
    atype = [ re.split(",", str(i))[0] for i in atype]
    atypes = ['Control', 'IL-1', 'IL-6']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJura2008 = getJura2008


def getGustafsson2008(self, tn=1):
    self.prepareData("MACV10")
    atype = self.h.getSurvName('c tissue')
    atypes = ['Decidual', 'Blood']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGustafsson2008 = getGustafsson2008


def getAndreu2017(self, tn=1):
    self.prepareData("MACV11")
    atype = self.h.getSurvName('c infection')
    atypes = ['U', 'D', 'L']
    ahash = {'uninfected':0, 'Dead (Irradiated) Mtb':1, 'Live Mtb':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAndreu2017 = getAndreu2017


def getBurke2020(self, tn=1):
    self.prepareData("MACV12")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Irf3/7', 'Irf1', 'Irf5']
    ahash = {'WT BMDMs':0, 'Irf3/7-/- BMDMs':1,
            'Irf1-/- BMDMs':2, 'Irf5-/- BMDMs':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurke2020 = getBurke2020


def getAaronson2017(self, tn=1):
    self.prepareData("MACV13")
    atype = self.h.getSurvName('c src1')
    atypes = ['Bv', 'Bi', 'Xv', 'Xi', 'Sv', 'Si', 'Lv', 'Li']
    ahash = {'Cerebellum, Vehicle':0,
            'Cerebellum, CHDI-00340246, 10mg/kg':1,
            'Cortex, Vehicle':2,
            'Cortex, CHDI-00340246, 10mg/kg':3,
            'Striatum, Vehicle':4,
            'Striatum, CHDI-00340246, 10mg/kg':5,
            'Liver, Vehicle':6,
            'Liver, CHDI-00340246, 10mg/kg':7}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAaronson2017 = getAaronson2017


def getFensterl2012(self, tn=1):
    self.prepareData("MACV14")
    atype = self.h.getSurvName('c Title')
    atypes = ['WV6', 'IV6', 'WV2', 'IV2']
    ahash = {'wtbrain-VSV-6d-rep1':0,
            'Ifit2KObrain-VSV-6d-rep1':1,
            'wtbrain-VSV-2d-rep1':2,
            'Ifit2KObrain-VSV-2d-rep1':3}
    if (tn == 2):
        atypes = ['W', '', 'Ifit2']
        ahash = {'wtbrain-VSV-2d-rep1':0,
                'Ifit2KObrain-VSV-2d-rep1':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFensterl2012 = getFensterl2012


def getHorsch2015(self, tn=1):
    self.prepareData("MACV15")
    atype = self.h.getSurvName('c genotype/treatment')
    atypes = ['Wt', 'Mut']
    ahash = {'homozygote':1, 'wild type':0}
    if (tn == 2):
        atype = self.h.getSurvName('c treatment protocol')
        atypes = ['C', 'T']
        ahash = {'control':0, 'OVA challenge':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHorsch2015 = getHorsch2015


def getHorsch2015II(self, tn=1):
    self.prepareData("MACV15.2")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Mut']
    ahash = {'Cox4-2 mutant line':1, 'wild type':0}
    if (tn == 2):
        atype = self.h.getSurvName('c treatment')
        atype = [ re.split(" ", str(i))[0] for i in atype]
        atypes = ['C', 'T']
        ahash = {'challenged':1, 'none':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHorsch2015II = getHorsch2015II


def getHorsch2015III(self, tn=1):
    self.prepareData("MACV15.3")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Mut']
    ahash = {'Prdm11 knockout mice':1, 'wild type':0}
    if (tn == 2):
        atype = self.h.getSurvName('c treatment')
        atype = [ re.split(" ", str(i))[0] for i in atype]
        atypes = ['C', 'T']
        ahash = {'challenged':1, 'none':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHorsch2015III = getHorsch2015III


def getLee2012(self, tn=1):
    self.prepareData("MACV16")
    atype = self.h.getSurvName('c aim2 expression')
    atypes = ['Wt', 'Mut']
    ahash = {'persistent expression':1, 'absent expression':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLee2012 = getLee2012


def getVasudevan2019(self, tn=1):
    self.prepareData("MACV17")
    atype = self.h.getSurvName('c genotype')
    atypes = ['Wt', 'Mut']
    ahash = {'control':0, 'Sp110 knockdown':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVasudevan2019 = getVasudevan2019


def getIrey2019Hs(self, tn=1):
    self.prepareData("MACV18")
    atype = self.h.getSurvName('c Title')
    atype = [ re.split(",", str(i))[0] for i in atype]
    atypes = ['V', '7', '231']
    ahash = {'MCF7-conditioned-medium':1,
            'vehicle-conditioned-medium':0,
            'MDA-MB-231-conditioned-medium':2}
    media = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName('c Title')
        atype = [ re.split(", ", str(atype[i]))[1] \
                if len(re.split(", ", str(atype[i]))) > 1 and \
                media[i] != 2 \
                else "" for i in range(len(atype))]
        atypes = ['V', 'I', 'I']
        ahash = {'vehicle':0, 'ruxolitinib':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getIrey2019Hs = getIrey2019Hs


def getIrey2019Mm(self, tn=1):
    self.prepareData("MACV18.2")
    atype = self.h.getSurvName('c phenotype')
    atypes = ['Wt', '', 'STAT3']
    ahash = {'STAT3 wild type':0, 'STAT3 knockout':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getIrey2019Mm = getIrey2019Mm


def getOakes2017Hs(self, tn=1):
    self.prepareData("MACV19")
    atype = self.h.getSurvName('c src1')
    atypes = ['W-D', 'W+D', 'M-D', 'M+D']
    ahash = {'T47D cell line, WT mouse Oas2, -DOX':0,
            'T47D cell line, MUT mouse Oas2, +DOX':3,
             'T47D cell line, WT mouse Oas2, +DOX':1,
             'T47D cell line, MUT mouse Oas2, -DOX':2}
    if (tn == 2):
        atypes = ['W-D', 'W+D']
        ahash = {'T47D cell line, WT mouse Oas2, -DOX':0,
                 'T47D cell line, WT mouse Oas2, +DOX':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOakes2017Hs = getOakes2017Hs


def getOakes2017Mm(self, tn=1):
    self.prepareData("MACV19.2")
    atype = self.h.getSurvName('c src1')
    atypes = ['W2', 'W18', 'M2', 'M18']
    ahash = {'Mammary gland, MUT Oas2, 2dpp':2,
            'Mammary gland, WT Oas2, 2dpp':0,
            'Mammary gland, WT Oas2, 18dpc':1,
            'Mammary gland, MUT Oas2, 18dpc':3}
    if (tn == 2):
        atypes = ['W2', '', 'M2']
        ahash = {'Mammary gland, MUT Oas2, 2dpp':2,
                'Mammary gland, WT Oas2, 2dpp':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOakes2017Mm = getOakes2017Mm


def getLinke2017(self, tn=1):
    self.prepareData("MACV20")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Mut']
    ahash = {'Tsc2fl/fl LysM+/+':0, 'Tsc2fl/fl LysM+/cre':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLinke2017 = getLinke2017


def getQi2017(self, tn=1):
    self.prepareData("MACV21")
    atype = self.h.getSurvName('c src1')
    atype = [ re.sub(".*\((.*)\).*", "\\1", str(i)) for i in atype]
    atypes = ['NC', 'KD']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQi2017 = getQi2017


def getLi2019(self, tn=1):
    self.prepareData("MACV22")
    atype = self.h.getSurvName('c mouse genotype')
    atypes = ['Wt', 'Rnf5']
    ahash = {'WT':0, 'Rnf5 KO':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2019 = getLi2019


def getScortegagna2020(self, tn=1):
    self.prepareData("MACV23")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Siah2']
    ahash = {'Siah2 WT':0, 'Siah2 KO':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getScortegagna2020 = getScortegagna2020


def getGoldmann2015(self, tn=1):
    self.prepareData("MACV24")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Ifnar', 'Usp18', 'C61A', 'KD', 'NC', 'DKO']
    ahash = {'WT':0,
            'IFNARko':1,
            'Usp18_C61A':3,
            'Usp18 ko':2,
            'si RNA Usp18':4,
            'si RNA control':5,
            'USP18ko:IFNARko':6}
    if (tn == 2):
        atypes = ['Wt', '', 'DKO']
        ahash = {'WT':0,
                'USP18ko:IFNARko':2}
    if (tn == 3):
        atypes = ['Ifnar', '', 'DKO']
        ahash = {'IFNARko':0,
                'USP18ko:IFNARko':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGoldmann2015 = getGoldmann2015


def getKurata2011(self, tn=1):
    self.prepareData("MACV25")
    atype = self.h.getSurvName('c Title')
    atypes = ['Wt', 'S', 'U']
    ahash = {'32Dcl3 Xbp1S':1, '32Dcl3 pMIG':0, '32Dcl3 Xbp1U':2}
    if (tn == 2):
        atypes = ['Wt', 'Xbp1']
        ahash = {'32Dcl3 Xbp1S':1, '32Dcl3 pMIG':0, '32Dcl3 Xbp1U':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKurata2011 = getKurata2011


def getWang2019Mac(self, tn=1):
    self.prepareData("MACV26")
    atype = self.h.getSurvName('c treatment')
    atypes = ['None', 'IL-15', 'IL-4', 'Media']
    ahash = {}
    if (tn == 2):
        atypes = ['None', 'IL-15']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2019Mac = getWang2019Mac


def getUckelmann2020(self, tn=1):
    self.prepareData("MACV27")
    rtype = self.h.getSurvName('c cell line')
    ahash = {'IMS-M2':0, 'OCI-AML-3':1}
    line = [ahash[i] if i in ahash else None for i in rtype]
    ttype = self.h.getSurvName('c timepoint')
    ahash = {'day 3':3, 'day 5':5, 'day 7':7}
    time = [ahash[i] if i in ahash else None for i in ttype]
    atype = self.h.getSurvName('c treatment')
    ahash = {'':0, '330nM VTP50469':1}
    treat = [ahash[i] if i in ahash else None for i in atype]
    atype = ["-".join([str(line[i]), str(treat[i]), str(time[i])])
                                for i in range(len(atype))]
    atypes = ['N', 'T']
    ahash = {'0-0-3':0, '0-0-5':0, '0-0-7':0,
            '1-0-3':0, '1-0-5':0, '1-0-7':0,
            '0-1-3':1, '0-1-5':1, '0-1-7':1,
            '1-1-3':1, '1-1-5':1, '1-1-7':1}
    if (tn == 2):
        ahash = {'1-0-3':0, '1-0-5':0, '1-0-7':0, '1-1-5':1}
    if (tn == 3):
        ahash = {'0-0-3':0, '0-0-5':0, '0-0-7':0,
                '0-1-3':1, '0-1-5':1, '0-1-7':1}
    if (tn == 4):
        ahash = {'1-0-3':0, '1-0-5':0, '1-0-7':0,
                '1-1-3':1, '1-1-5':1, '1-1-7':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getUckelmann2020 = getUckelmann2020


def getUckelmann2020Mm(self, tn=1):
    self.prepareData("MACV28")
    rtype = self.h.getSurvName('c cell type')
    rtype = [str(i).replace('mouse ', '').replace(' cells', '') \
            for i in rtype]
    ahash = {'':0, 'LSK':1, 'GMP':2, 'LSK-derived GMP':3,
            'GMP-derived GMP':4, 'LSK-derived LSK':5,
            'GMP-derived LSK':6, 'long-term GMP-derived GMP':7,
            'LSK-derived GMP-like':8, 'GMP-derived GMP-like':9}
    line = [ahash[i] if i in ahash else None for i in rtype]
    ttype = self.h.getSurvName('c timepoint')
    ahash = {'4 weeks':28, '':0, '9 months post transplant':270,
            'day 5':5, 'day 3':3}
    time = [ahash[i] if i in ahash else None for i in ttype]
    atype = self.h.getSurvName('c treatment')
    ahash = {'pIpC induction':1, '':0, '1% VTP50469 chow':2, 
            '30nM VTP50469':3}
    treat = [ahash[i] if i in ahash else None for i in atype]
    ctype = ["-".join([str(line[i]), str(treat[i]), str(time[i])])
                                for i in range(len(atype))]
    atype = treat
    atypes = [0, 1, 2, 3]
    ahash = {}
    if (tn == 2):
        atype = [treat[i] if (line[i] == 0) else None
                for i in range(len(atype))]
        atypes = [0, 3]
    if (tn == 3):
        atype = [treat[i] if (line[i] == 2) else None
                for i in range(len(atype))]
        atypes = [0, 1, 2, 3]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getUckelmann2020Mm = getUckelmann2020Mm


def getDong2019(self, tn=1):
    self.prepareData("MACV29")
    atype = self.h.getSurvName('c genotype')
    atypes = ['W', 'M']
    ahash = {'IRE1 KO':1, 'C57BL/6':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDong2019 = getDong2019


def getGlmozzi2019 (self, tn=1):
    self.prepareData("MACV30")
    atype = self.h.getSurvName('c genotype')
    atypes = ['W', 'M']
    ahash = {'wild type':0, 'PGRMC2 adipose tissue knockout':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlmozzi2019  = getGlmozzi2019 


def getElDahr2019 (self, tn=1):
    self.prepareData("MACV31")
    atype = self.h.getSurvName('c src1')
    atype = [re.sub("\s*[0-9]$", "", str(i)) for i in atype]
    atypes = ['W', '2', '21het', '21']
    ahash = {'WT':0,
            'EZH2-KO & EZH1-WT':1,
            'EZH2-KO & EZH1-hetKO':2,
            'EZH2-KO & EZH1-homoKO':3}
    if (tn == 2):
        atypes = ['W', 'M']
        ahash = {'WT':0,
                'EZH2-KO & EZH1-WT':0,
                'EZH2-KO & EZH1-hetKO':1,
                'EZH2-KO & EZH1-homoKO':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getElDahr2019  = getElDahr2019 


def getEncode2020 (self, tn=1):
    self.prepareData("MACV32")
    atype = self.h.getSurvName('c Type')
    ahash = {'K562':1, 'HepG2':2, '':0}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Target')
    atypes = ['W', 'M']
    ahash = {'':0, 'EEF2':1, 'NFX1':1, 'HMBOX1':1, 
            'HNRNPA1':1, 'NFYB':1, 'PCBP2':1}
    for k in atype:
        if k != '':
            ahash[k] = 1
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None
                for i in range(len(atype))]
        ahash = {'':0, 'EEF2':1, 'NFX1':1, 'HMBOX1':1, 
                'HNRNPA1':1, 'NFYB':1, 'PCBP2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEncode2020  = getEncode2020 


def getBeheshti2015(self, tn=1):
    self.prepareData("MACV33")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'T']
    ahash = {'Spleen from LLC tumor bearing mice':1,
            'Spleen from non-tumor bearing control mice':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBeheshti2015 = getBeheshti2015


def getPlatten2020(self, tn=1):
    self.prepareData("MACV34")
    atype = self.h.getSurvName('c phenotype')
    atypes = ['R', 'NR']
    ahash = {'PD-1 and CTLA-4 non-responder':1,
            'PD-1 and CTLA-4 responder':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPlatten2020 = getPlatten2020


def getVolkmer2020(self, tn=1):
    self.prepareData("MACV35")
    atype = self.h.getSurvName('c src1')
    atype = [ re.split(",", str(i))[0] for i in atype]
    ahash = {'MC38 whole tumor':0,
            'B16-OVA whole tumor':1,
            'PDA30364 cell line pellet':2,
            'PDA30364 whole tumor':3,
            'B16-OVA cell line pellet':4,
            'MC38 cell line pellet':5}
    tumor = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment groups')
    atypes = ['C1', 'T1', 'C2', 'T2', 'C3', 'T3', 'C4', 'T4', 'C5', 'T5']
    ahash = {'24h DMSO':0, '24h GDC-0623':1,
            '72h DMSO':2, '72h GDC-0623':3,
            'Vehicle+control IgG':4, 'Vehicle+CD40 mIgG1':5,
            'GEM+control IgG':6, 'GEM+CD40 mIgG1':7,
            'MEKi+control IgG':8, 'MEKi+CD40 mIgG1':9}
    treat = [ahash[i] if i in ahash else None for i in atype]
    if (tn >= 2):
        atype = ["-".join([str(tumor[i]), str(treat[i])])
                                for i in range(len(atype))]
        atypes = ['C', 'T']
        list1 = ['0-4', '0-5', '0-8', '0-9',
                '1-4', '1-5', '1-8', '1-9',
                '2-0', '2-1', '2-2', '2-3',
                '3-4', '3-5', '3-6', '3-7',
                '4-0', '4-1', '4-2', '4-3',
                '5-0', '5-1', '5-2', '5-3']
        index = (tn - 2) * 2;
        ahash = { list1[index] : 0, list1[index + 1]:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVolkmer2020 = getVolkmer2020


def getVolkmer2020II(self, tn=1):
    self.prepareData("MACV36")
    atype = self.h.getSurvName('c src1')
    atype = [ re.split(",", str(i))[0] for i in atype]
    ahash = {'MC38 whole tumor':0,
            'B16-OVA whole tumor':1,
            'PDA30364 cell line pellet':2,
            'PDA30364 whole tumor':3,
            'B16-OVA cell line pellet':4,
            'MC38 cell line pellet':5}
    tumor = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment groups')
    atypes = ['C1', 'T1', 'C2', 'T2', 'C3', 'T3', 'C4', 'T4', 'C5', 'T5']
    ahash = {'24h DMSO':0, '24h GDC-0623':1,
            '72h DMSO':2, '72h GDC-0623':3,
            'Vehicle+control IgG':4, 'Vehicle+CD40 mIgG1':5,
            'GEM+control IgG':6, 'GEM+CD40 mIgG1':7,
            'MEKi+control IgG':8, 'MEKi+CD40 mIgG1':9}
    treat = [ahash[i] if i in ahash else None for i in atype]
    if (tn >= 2):
        atype = ["-".join([str(tumor[i]), str(treat[i])])
                                for i in range(len(atype))]
        atypes = ['C', 'T']
        list1 = ['0-4', '0-5', '0-8', '0-9',
                '1-4', '1-5', '1-8', '1-9',
                '2-0', '2-1', '2-2', '2-3',
                '3-4', '3-5', '3-6', '3-7',
                '4-0', '4-1', '4-2', '4-3',
                '5-0', '5-1', '5-2', '5-3']
        index = (tn - 2) * 2;
        ahash = { list1[index] : 0, list1[index + 1]:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVolkmer2020II = getVolkmer2020II


def getZhou2020(self, tn=1):
    self.prepareData("MACV37")
    atype = self.h.getSurvName('c treatment')
    atypes = ['T1', 'T2']
    ahash = {'anti-gp120':0, 'anti-MerTK':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhou2020 = getZhou2020


def getZhou2020II(self, tn=1):
    self.prepareData("MACV38")
    atype = self.h.getSurvName('c treatment')
    atypes = ['T1', 'T2']
    ahash = {'anti-gp120':0, 'anti-MerTK':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhou2020II = getZhou2020II


def getSteenbrugge2019(self, tn=1):
    self.prepareData("MACV39")
    atype = self.h.getSurvName('c src1')
    atypes = ['N', 'T', 'N1', 'T1']
    ahash = {'C57BL/6-derived mammary gland':0,
            '4T1 tumor':1,
            'BALB/c-derived mammary gland':2,
            'Py230 tumor':3}
    if (tn == 2):
        atypes = ['N', 'T']
        ahash = {'C57BL/6-derived mammary gland':0,
            'Py230 tumor':1}
    if (tn == 3):
        atypes = ['N', 'T']
        ahash = {'BALB/c-derived mammary gland':0,
                '4T1 tumor':1}
    if (tn == 4):
        atypes = ['BALB/c', 'C57BL/6']
        ahash = {'BALB/c-derived mammary gland':0,
                'C57BL/6-derived mammary gland':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSteenbrugge2019 = getSteenbrugge2019


def getHollern2019(self, tn=1):
    self.prepareData("MACV40")
    atype = self.h.getSurvName('c class')
    atypes = ['S', 'R']
    ahash = {'sensitive':0, 'resistant':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHollern2019 = getHollern2019


def getDas2015(self, tn=1):
    self.prepareData("MACV41")
    atype = self.h.getSurvName('c agent')
    atypes = ['aCTLA4', 'Combo', 'Seq', 'aPD1']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName('c time point')
        atypes = ['pre', 'post', 'pre-s', 'post-s']
        ahash = {'post':1, 'pre':0,
                'pre sequential therapy sample':2,
                'post sequential therapy sample':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDas2015 = getDas2015


def getTaube2015(self, tn=1):
    self.prepareData("MACV42")
    atype = self.h.getSurvName('c pd-l1 status')
    atypes = ['pos', 'neg']
    ahash = {'PD-L1 positive':0, 'PD-L1 negative':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTaube2015 = getTaube2015


def getBalachandran2019(self, tn=1):
    self.prepareData("MACV43")
    atype = self.h.getSurvName('c genotype')
    atypes = ['W', 'M']
    ahash = {'wild type':0, 'IL33-/-':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBalachandran2019 = getBalachandran2019


def getHwang2020(self, tn=1):
    self.prepareData("MACV44")
    atype = self.h.getSurvName('c tumor type')
    atypes = ['L', 'A', 'AS', 'S']
    ahash = {'Large cell neuroendocarine carcinoma':0,
            'adenocarcinoma':1,
            'adenocarcinoma-squamous cell carcinoma':2,
            'squamous cell carcinoma':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHwang2020 = getHwang2020


def getSimpson2019(self, tn=1):
    self.prepareData("MACV45")
    atype = self.h.getSurvName('c patient outcome')
    atypes = ['C', 'I', 'F']
    ahash = {'control':0,
            'acute liver failure':2,
            'acute liver injury':1}
    if (tn == 2):
        atype = self.h.getSurvName('c survival')
        atypes = ['S', 'D']
        ahash = {'spontaneously survived':0, 'dead or transplanted':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSimpson2019 = getSimpson2019


def getHuang2019(self, tn=1):
    self.prepareData("MACV46")
    atype = self.h.getSurvName('c fvc_progressor_10%')
    atypes = ['0', '1']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName('c dlco_progressor_15%')
        atypes = ['0', '1']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHuang2019 = getHuang2019


def getYu2019(self, tn=1):
    self.prepareData("MACV47")
    course = self.h.getSurvName('c course')
    ahash = {'Acute myocarditis':0,
            'Acute myocardial infarction':1,
            'Aortic dissection':2,
            'Congestive heart failure':3,
            'Dilated cardiomyopathy':4,
            'Dilated cardiomyopathy, DCMP':5,
            'Arrhythmia':6}
    rval = [ahash[i] if i in ahash else None for i in course]
    atype = self.h.getSurvName('c outcome')
    atypes = ['S', 'F']
    ahash = {'Success':0, 'Failure':1, 'failure':1}
    if (tn == 2):
        atype = course
        atypes = ['Am', 'Ami', 'Ad', 'Chf', 'Dc', 'Dcmp', 'Ar']
        ahash = {'Acute myocarditis':0, 'Acute myocardial infarction':1,
                'Aortic dissection':2, 'Congestive heart failure':3,
                'Dilated cardiomyopathy':4, 'Dilated cardiomyopathy, DCMP':5,
                'Arrhythmia':6}
    if (tn == 3):
        atype = [atype[i] if rval[i] != 4 and rval[i] != 5 else None
                for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if rval[i] == 1 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYu2019 = getYu2019


def getVanhaverbeke2019(self, tn=1):
    self.prepareData("MACV48")
    time = self.h.getSurvName('c time point')
    status = self.h.getSurvName('c patient diagnosis')
    atype = [" ".join([str(status[i]), str(time[i])])
            for i in range(len(time))]
    atypes = ['MI D0', 'MI D30', 'MI Y1', 'CAD D0']
    ahash = {}
    if (tn == 2):
        atypes = ['CAD D0', 'MI D0']
    if (tn == 3):
        atypes = ['MI D30', 'MI D0']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanhaverbeke2019 = getVanhaverbeke2019


def getSuresh2014(self, tn=1):
    self.prepareData("MACV49")
    atype = self.h.getSurvName('c disease status')
    atypes = ['C', 'Nr', 'r']
    ahash = {'normal control':0,
            'patient without recurrent events':1,
            'patient with recurrent events':2}
    if (tn == 2):
        atypes = ['C', 'r']
        ahash = {'normal control':0,
                'patient with recurrent events':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSuresh2014 = getSuresh2014


def getPalmer2019(self, tn=1, tb=0):
    self.prepareData("MACV50")
    atype = self.h.getSurvName('c src1')
    ahash = {'colon':0, 'blood':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease type')
    atypes = ['C', 'IBD', 'UC', 'CD']
    ahash = {'Control':0,
            'Control - infect':0,
            'Control - celiac':0,
            'Control (E.coli) --> CD':0,
            'Control--> CD':0,
            'IBD w/ oral':1,
            'IBDU':1,
            'UC - pancolitis':2,
            'UC - L colitis':2,
            'UC - proctitis':2,
            "Crohn's Disease":3}
    if (tn == 2):
        atypes = ['C', 'IBD']
        ahash = {'Control':0,
                'Control - infect':0,
                'Control - celiac':0,
                'Control (E.coli) --> CD':0,
                'Control--> CD':0,
                'UC - pancolitis':1,
                'UC - L colitis':1,
                'UC - proctitis':1,
                "Crohn's Disease":1}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPalmer2019 = getPalmer2019


def getOstrowski2019(self, tn=1):
    self.prepareData("MACV51")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'UC', 'CD', 'PSC', 'PBC', 'UCc', 'CDc']
    ahash = {'Control, adult':0,
            'Ulcerative colitis, adult':1,
            'Crohn\xe2\x80\x99s disease, adult':2,
            'Primary sclerosing cholangitis, adult':3,
            'Primary biliary cholangitis, adult':4,
            'Ulcerative colitis, child':5,
            'Crohn\xe2\x80\x99s disease, child':6}
    if (tn == 2):
        atypes = ['C', 'UC', 'CD']
        ahash = {'Control, adult':0,
                'Ulcerative colitis, adult':1,
                'Crohn\xe2\x80\x99s disease, adult':2}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOstrowski2019 = getOstrowski2019


def getJeffrey2006(self, tn=1):
    self.prepareData("MACV52")
    src = self.h.getSurvName('c src1')
    ahash = {'Cord blood':0, 'Peripheral blood':1}
    rval = [ahash[i] if i in ahash else None for i in src]
    title = self.h.getSurvName('c Title')
    atype = [1 if str(k).find("unstimulated") >= 0 or \
            str(k).find("control") >= 0 or \
            str(k).find("Immature") >= 0 else 0 for k in title]
    atypes = ['C', 'S']
    ahash = {0:1, 1:0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if tn == 3:
        atype = self.h.getSurvName("c Title")
        atype = [bone.re.sub("_.*", "", str(k)) for k in atype]
        atypes = ['Other', 'Neu', 'Mac', 'Eos', 'Bas', 'NK']
        ahash = {'Neutrophils':1, 'Cord blood-derived mast cells':0,
                 'Th2 cells':0, 'Th1 cells':0, 'Dendritic cells':0,
                 'Immature dendritic cells':0, 'B cells':0, 'NK cells':5,
                 'Macrophages':2, 'Central memory T cells':0, 'Basophils':4,
                 'Eosinophils':3, 'Effector memory T cells':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJeffrey2006 = getJeffrey2006


def getSippel2018(self, tn=1):
    self.prepareData("MACV53")
    atype = self.h.getSurvName('c treatment')
    atypes = ['V', 'IL5', 'PGD2']
    ahash = {'IL5':1, 'Vehicle':0, 'dkPGD2':2}
    if (tn == 2):
        atypes = ['V', 'T']
        ahash = {'IL5':1, 'Vehicle':0, 'dkPGD2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSippel2018 = getSippel2018


def getPuan2017(self, tn=1):
    self.prepareData("MACV54")
    state = self.h.getSurvName('c donor_state')
    ahash = {'reactive':0, 'anergic':1}
    rval = [ahash[i] if i in ahash else None for i in state]
    atype = self.h.getSurvName('c stimulation')
    atypes = ['C', 'S']
    ahash = {'unstimulated':0, 'Fc-epsilon receptor-crosslinking':1}
    self.pair = [ [2, 11], [3, 6], [7, 5], [4, 9], [10, 8]]
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPuan2017 = getPuan2017


def getKlocperk2020(self, tn=1):
    self.prepareData("MACV55")
    src = self.h.getSurvName('c src1')
    src = [re.sub("mo.*of ", "", str(k)) for k in src]
    treat = [re.sub(".* [1-5] *", "", str(k)) for k in src]
    treat = [re.sub("c.* with ", "", str(k)) for k in treat]
    ahash = {'':0,
            'autologous healthy NETs':1,
            'healthy NETs':2,
            'T1D NETs':3,
            'autologous T1D NETs':4}
    rval = [ahash[i] if i in ahash else None for i in treat]
    disease = [k.split(" ")[0] for k in src]
    atype = disease
    atypes = ['healthy', 'T1D']
    ahash = {}
    if (tn == 2):
        atypes = ['C', 'S']
        atype = rval
        ahash = {0:0, 1:1, 2:1, 3:1, 4:1}
        atype = [atype[i] if disease[i] == 'healthy' \
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'S']
        atype = rval
        ahash = {0:0, 1:1, 2:1, 3:1, 4:1}
        atype = [atype[i] if disease[i] == 'T1D' \
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKlocperk2020 = getKlocperk2020


def getLenardo2020(self, tn=1):
    self.prepareData("MACV58")
    atype = self.h.getSurvName('c time')
    atypes = ['D2', ' ', 'D5']
    ahash = {'Day2':0, 'Day5':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLenardo2020 = getLenardo2020


def getNair2015(self, tn=1):
    self.prepareData("MACV59")
    atype = self.h.getSurvName('c protocol')
    atypes = ['C', 'S']
    ahash = {'unstimulated (control)':0, 'stimulated':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNair2015 = getNair2015


def getJohnson2020(self, tn=1):
    self.prepareData("MACV60")
    atype = self.h.getSurvName("c treatment_code1")
    atypes = ['C', 'S']
    ahash = {'NS':0}
    for k in atype:
        if k != 'NS':
            ahash[k] = 1
    if (tn == 2):
        atypes = ['C', 'HIV']
        ahash = {'NS': 0, 'HIV2_WT':1, 'HIV2_P86HA':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJohnson2020 = getJohnson2020


def getAbbas2005(self, tn=1):
    self.prepareData("MACV61")
    atype = self.h.getSurvName("c src1")
    ahash = {'NK cells from PBMC':0,
            'Plasma cells from bone marrow':1,
            'Monocytes from PBMC':2,
            'CD4+ T cells from PBMC':3,
            'B cells from PBMC':4,
            'Neutrophils from PBMC':5,
            'CD14+ cells from PBMC':6,
            'CD4+ CD45RO+ CD45RA- T cells from PBMC':7,
            'CD8+ T cells from PBMC':8,
            'Plasma cells from PBMC':9}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment agent')
    ahash = {'NA':0,
            'macrophage differentiation medium':1,
            'IL2':2,
            'LPS':3,
            'IL15':4,
            'aCD3/aCD28':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = tval
    atypes = ['C', 'S']
    ahash = {0:0, 1:1, 2:1, 3:1, 4:1, 5:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAbbas2005 = getAbbas2005


def getMetcalf2015(self, tn=1):
    self.prepareData("MACV62")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'S']
    ahash = {'Rig I':1, 'PolyIC':1, 'NoTx':0, 'LyoVec_only':1, 'LPS':1, 'CLO_97':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMetcalf2015 = getMetcalf2015


def getBanchereau2014I(self, tn=1):
    self.prepareData("MACV63")
    atype = self.h.getSurvName('c cell population')
    ahash = {'IL4 DC':1, 'IFNa DC':0}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c stimulation')
    atypes = ['C', 'S']
    ahash = {'MDP':1, 'LPS':1, 'Poly I:C-LMW-Lyovec':1,
            'TNFa':1, 'CpG 2216':0, 'Poly I:C':1, 'R837':1,
            'CL097':1, 'IFNa':1, 'IL10':1, 'CpG 2006':0, 'Flagellin':1,
            'PAM3':1, 'IL15':1, 'IL1b':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        atype = rval
        ahash = {0:1, 1:2}
        atype = [atype[i] if aval[i] == 0 \
                else None for i in range(len(atype))]
    if (tn == 3):
        ahash['A'] = 1
        atype = [atype[i] if rval[i] == 1 \
                else 'A' for i in range(len(atype))]
    if (tn == 4):
        atypes = ['C', 'S']
        ahash = {'CpG 2216':0, 'CpG 2006':0, 'IL1b':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBanchereau2014I = getBanchereau2014I


def getBanchereau2014II(self, tn=1):
    self.prepareData("MACV64")
    atype1 = self.h.getSurvName('c culture conditions')
    atype2 = self.h.getSurvName('c culture condition')
    atype3 = self.h.getSurvName('c vaccine abbr.')
    atype = [ " ".join([str(k) for k in [atype1[i], atype2[i], atype3[i]]]) 
                     for i in range(len(atype1))]
    atypes = ['C', 'S']
    ahash = {'  RAB':1, 'LPS  ':1, ' HKSE ':1,
            ' Media ':0, '  Medium':0, 'Medium1  ':0,
            '  FZ':1, '  TDAP':1, ' H1N1 Brisbane ':1, ' HKSA ':1, 'Medium2  ':1,
            'HEPB  ':1, 'HPV  ':1, 'HIB  ':1, '  POL':1, '  VAR':1, 'VAR  ':1,
            'PVX  ':1, 'RAB  ':1, 'MGL  ':1, '  HIB':1, '  HER':1, '  HPV':1,
            'HER  ':1, '  PVX':1, '  JPE':1, '  HEPB':1, 'HEPA  ':1, 'JPE  ':1,
            '  HEPA':1, 'FZ  ':1, 'POL  ':1, '  MGL':1, 'TDAP  ':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName('c cell population')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'IL4 DC':2, 'IFNa DC':1, 'Monocytes':0}
        atype = [atype[i] if aval[i] == 0 \
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBanchereau2014II = getBanchereau2014II


def getHarfuddin2016(self, tn=1):
    self.prepareData("MACV65")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'S']
    ahash = {'None':0, 'CD137-Fc protein':1, 'GM-CSF + IL-4':1,
            'GM-CSF + IL-4; matured with LPS + IFNg':1,
            'Fc protein':1, 'M-CSF':1}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'None':0, 'GM-CSF + IL-4':2,
                'GM-CSF + IL-4; matured with LPS + IFNg':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHarfuddin2016 = getHarfuddin2016


def getIampietro2017(self, tn=1):
    self.prepareData("MACV66")
    atype = self.h.getSurvName('c infection')
    atypes = ['Mock', 'EBOV', 'LPS']
    ahash = {}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'Mock':0, 'EBOV':1, 'LPS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getIampietro2017 = getIampietro2017


def getSinnaeve2009(self, tn=1):
    self.prepareData("MACV68")
    atype = self.h.getSurvName('c CADi')
    atypes = ['C', 'D']
    ahash = {'0':0}
    for k in atype:
        if k != '0':
            ahash[k] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSinnaeve2009 = getSinnaeve2009


def getMartinez2019(self, tn=1):
    self.prepareData("MACV69")
    atype = self.h.getSurvName('c outcome')
    atypes = ['M', 'F']
    ahash = {'Matured':0, 'Failed':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinez2019 = getMartinez2019


def getHollander2013(self, tn=1):
    self.prepareData("MACV70")
    atype = self.h.getSurvName('c src1')
    atypes = ['non-AR', 'AR']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHollander2013 = getHollander2013


def getBaron2007I(self, tn=1):
    self.prepareData("MACV71")
    atype = self.h.getSurvName('c Title')
    cells = [re.split("_", str(i))[1] if len(re.split("_", str(i))) > 1
                            else None for i in atype]
    acute = [re.split("_", str(i))[3] if len(re.split("_", str(i))) > 3
                            else None for i in atype]
    chronic = [re.split("_", str(i))[4] if len(re.split("_", str(i))) > 4
                            else None for i in atype]
    atype = acute
    atypes = ['C', 'D']
    ahash = {'aGVHD+':1, 'aGVHD-':0}
    if (tn == 2):
        atype = chronic
        ahash = {'cGVHD+':1, 'cGVHD-':0}
    if (tn == 3):
        atype = [atype[i] if cells[i] == 'CD4+' \
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if cells[i] == 'CD8+' \
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = chronic
        ahash = {'cGVHD+':1, 'cGVHD-':0}
        atype = [atype[i] if cells[i] == 'CD4+' \
                else None for i in range(len(atype))]
    if (tn == 6):
        atype = chronic
        ahash = {'cGVHD+':1, 'cGVHD-':0}
        atype = [atype[i] if cells[i] == 'CD8+' \
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBaron2007I = getBaron2007I


def getBaron2007II(self, tn=1):
    self.prepareData("MACV72")
    atype = self.h.getSurvName('c Title')
    cells = [re.split("_", str(i))[1] if len(re.split("_", str(i))) > 1
                            else None for i in atype]
    acute = [re.split("_", str(i))[3] if len(re.split("_", str(i))) > 3
                            else None for i in atype]
    chronic = [re.split("_", str(i))[4] if len(re.split("_", str(i))) > 4
                            else None for i in atype]
    atype = acute
    atypes = ['C', 'D']
    ahash = {'aGVHD+':1, 'aGVHD-':0}
    if (tn == 2):
        atype = chronic
        ahash = {'cGVHD+':1, 'cGVHD-':0}
    if (tn == 3):
        atype = [atype[i] if cells[i] == 'CD4+' \
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if cells[i] == 'CD8+' \
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = chronic
        ahash = {'cGVHD+':1, 'cGVHD-':0}
        atype = [atype[i] if cells[i] == 'CD4+' \
                else None for i in range(len(atype))]
    if (tn == 6):
        atype = chronic
        ahash = {'cGVHD+':1, 'cGVHD-':0}
        atype = [atype[i] if cells[i] == 'CD8+' \
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBaron2007II = getBaron2007II


def getSagoo2010(self, tn=1):
    self.prepareData("MACV73")
    atypes = ["Tol-DF", "s-LP", "s-CNI", "s-nCNI", "CR", "HC"]
    atype = self.h.getSurvName("c Tol-DF")
    for k in atypes:
        g1 = self.h.getSurvName("c " + k)
        for i in range(len(g1)):
            if (g1[i] == '1'):
                atype[i] = k
    ahash = {}
    if (tn == 2):
        atypes = ['HC', 'SD', 'CR']
        ahash = {"Tol-DF":1, "s-LP":1, "s-CNI":1, "s-nCNI":1, "CR":2, "HC":0}
    if (tn == 3):
        atypes = ['HC', 'CR', 'SD']
        ahash = {"Tol-DF":2, "CR":1, "HC":0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSagoo2010 = getSagoo2010


def getHerazoMaya2013(self, tn=1):
    self.prepareData("MACV75")
    atype = self.h.getSurvName("c outcome")
    atypes = ['0', '1']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHerazoMaya2013 = getHerazoMaya2013


def getLi2012I(self, tn=1):
    self.prepareData("MACV76")
    atype = self.h.getSurvName("c disease state")
    atype = [str(k)[0:2] for k in atype]
    atypes = ['ST', 'AR']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2012I = getLi2012I


def getLi2012II(self, tn=1):
    self.prepareData("MACV77")
    atype = self.h.getSurvName("c Disease State")
    atype = [str(k)[0:1] for k in atype]
    atypes = ['S', 'A']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2012II = getLi2012II


def getKurian2014(self, tn=1):
    self.prepareData("MACV79")
    atype = self.h.getSurvName("c phenotype")
    atypes = ['S', 'R']
    ahash = {'Acute Kidney Rejection':1,
            'Well-functioning kidney transplant':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKurian2014 = getKurian2014


def getZhang2019(self, tn=1):
    self.prepareData("MACV80")
    atype = self.h.getSurvName("c acr before or at 6m")
    atypes = ['S', 'R']
    ahash = {'ACR or Borderline':1, 'None':0}
    if (tn == 2):
        atype = self.h.getSurvName("c acr after 6m")
        atypes = ['S', 'B', 'R']
        ahash = {'Borderline':1, 'None':0, 'ACR':2}
    if (tn == 3):
        atype = self.h.getSurvName("c acr after 6m")
        atypes = ['S', 'R']
        ahash = {'Borderline':0, 'None':0, 'ACR':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2019 = getZhang2019


def getKhatri2013(self, tn=1):
    self.prepareData("MACV81")
    atype = self.h.getSurvName("c patient group")
    atypes = ['S', 'R']
    ahash = {'stable patient (STA)':0,
            'patient with acute rejection (AR)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKhatri2013 = getKhatri2013


def getEinecke2010(self, tn=1):
    self.prepareData("MACV82")
    atype = self.h.getSurvName("c rejection/non rejection")
    atypes = ['S', 'R']
    ahash = {'rej':1, 'nonrej':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEinecke2010 = getEinecke2010


def getReeve2013(self, tn=1):
    self.prepareData("MACV83")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['S', 'R']
    ahash = {'non-rejecting':0, 'MIXED':1, 'ABMR':1, 'TCMR':1, 'Nephrectomy':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getReeve2013 = getReeve2013


def getRay2007(self, tn=1):
    self.prepareData("MACV84")
    atype = self.h.getSurvName("c Title")
    atype = [re.split(" ", str(i))[2] if len(re.split(" ", str(i))) > 2
                                            else None for i in atype]
    atypes = ['S', 'R']
    ahash = {'not':0, 'PGD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRay2007 = getRay2007


def getVanLoon2019(self, tn=1):
    self.prepareData("MACV88")
    atype = self.h.getSurvName("c tissue")
    ahash = {'BLOOD':0, 'kidney allograft biopsy':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    tcmr = self.h.getSurvName("c tcmr (no:0_borderline:1_TCMR:2)")
    abmr = self.h.getSurvName("c abmr (no:0_Yes:1)")
    atype = [str(tcmr[i]) + " " + str(abmr[i]) for i in range(len(tcmr))]
    atypes = ['S', 'R']
    ahash = {'2 0':1, '0 1':1, '0 0':0, '1 0':1, '2 1':1, '1 1':1}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 \
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 \
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanLoon2019 = getVanLoon2019


def getMorgun2006I(self, tn=1):
    self.prepareData("MACV89")
    atype = self.h.getSurvName("c Final Clinical diagnosis")
    atypes = ['N', 'Ch', 'Pre-Ch', 'R', 'Pre-R', 'Tox']
    ahash = {'toxoplasma myocarditis':5}
    if (tn == 2):
        atypes = ['S', 'R']
        ahash = {'R':1, 'N':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMorgun2006I = getMorgun2006I


def getMorgun2006II(self, tn=1):
    self.prepareData("MACV90")
    atype = self.h.getSurvName("c Final Clinical diagnosis")
    atypes = ['N', 'R', 'Pre-R']
    ahash = {}
    if (tn == 2):
        atypes = ['S', 'R']
        ahash = {'R':1, 'N':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMorgun2006II = getMorgun2006II


def getShannon2012(self, tn=1):
    self.prepareData("MACV91")
    atype = self.h.getSurvName("c rejection status")
    atypes = ['NR', 'AR']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShannon2012 = getShannon2012


def getWells2003(self, tn=1):
    self.prepareData("MACV93")
    atype = self.h.getSurvName("c Title")
    strain = [str(k).split("_")[1] if len(str(k).split("_")) > 1 else None
            for k in atype]
    ahash = {'BalbC':0, 'C57BL6':3, 'C3H/ARC':2, 'C3H/HeJ':1, 'C57/BL6':4}
    rval = [ahash[i] if i in ahash else None for i in strain]
    atype = [re.sub(".*time", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'7h':1, '21h':None, '0.5h':0, '2h':1, '0h':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] >= 3 else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if rval[i] == 2 else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 6):
        atype = strain
        atypes = ['BalbC', 'C3H/HeJ', 'C3H/ARC', 'C57/BL6']
        ahash = {'BalbC':0, 'C57BL6':3, 'C3H/ARC':2, 'C3H/HeJ':1, 'C57/BL6':3}
    if (tn == 7):
        atype = strain
        atypes = ['BalbC', 'C57/BL6']
        ahash = {'BalbC':0, 'C57BL6':1, 'C57/BL6':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWells2003 = getWells2003


def getvanErp2006(self, tn=1):
    self.prepareData("MACV94")
    atype = self.h.getSurvName("c Title")
    strain = [str(k).split("_")[0] if len(str(k).split("_")) > 0 else None
            for k in atype]
    ahash = {'Bc':0, 'Bl6':1}
    rval = [ahash[i] if i in ahash else None for i in strain]
    treat = [str(k).split("_")[1] if len(str(k).split("_")) > 1 else None
            for k in atype]
    ahash = {'WAP':1, 'p60':2, 'MOCK':0, 'Mock':0}
    tval = [ahash[i] if i in ahash else None for i in treat]
    atype = [str(k).split("_")[2] if len(str(k).split("_")) > 2 else None
            for k in atype]
    atypes = ['M0', 'M1']
    ahash = {'IFNy':1, '1':0, '3':0, '2':0, '4':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = strain
        atypes = ['Bc', 'Bl6']
        ahash = {}
    if (tn == 5):
        atype = strain
        atype = [atype[i] if tval[i] == 0 else None for i in range(len(atype))]
        atypes = ['Bc', 'Bl6']
        ahash = {}
    self.rval = [ahash[i] if i in ahash else None for i in strain]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getvanErp2006 = getvanErp2006


def getTang2020(self, tn=1):
    self.prepareData("MACV95")
    atype = self.h.getSurvName("c macrophage phenotype")
    ahash = {'Tissue resident, F4/80hi CD206-':1,
            'Monocyte-derived, F4/80int CD206+':0}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Title")
    strain = [str(k).split("_")[0] if len(str(k).split("_")) > 0 else None
            for k in atype]
    ahash = {'B6':1, 'Balbc':0}
    rval = [ahash[i] if i in ahash else None for i in strain]
    atype = [str(k).split("_")[1] if len(str(k).split("_")) > 1 else None
            for k in atype]
    atypes = ['IL4c', 'ThioIL4c']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = strain
        atype = [atype[i] if mval[i] == 0 else None for i in range(len(atype))]
        atypes = ['Balbc', 'B6']
        ahash = {}
    if (tn == 5):
        atype = strain
        atype = [atype[i] if mval[i] == 1 else None for i in range(len(atype))]
        atypes = ['Balbc', 'B6']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTang2020 = getTang2020


def getSajti2020(self, tn=1):
    self.prepareData("MACV96.2")
    atype = self.h.getSurvName("c Title")
    group = [re.sub("(.*h[_in]*).*", "\\1", str(k)) for k in atype]
    atype = [str(k).split("_")[0] if len(str(k).split("_")) > 0 else None
            for k in atype]
    ahash = {'AM':0, 'IM':1, 'IMo':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment time")
    ahash = {'20h':20, '6h':6, '':0, '2h':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'LPSi', 'LPSn']
    ahash = {'nasal LPS':2, 'IP LPS':1, '':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if rval[i] == 2 else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if tval[i] <= 2 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSajti2020 = getSajti2020


def getSajti2020II(self, tn=1):
    self.prepareData("MACV96.3")
    strain = self.h.getSurvName("c strain")
    ahash = {'C57BL6/J':0, 'DBA2/J':1}
    rval = [ahash[i] if i in ahash else None for i in strain]
    atype = self.h.getSurvName("c Title")
    atype = [str(k).split(" ")[0] if len(str(k).split(" ")) > 0 else None
                            for k in atype]
    atypes = ['AM', 'IM', 'iMo', 'pMo']
    ahash = {'AM':0, 'IM':1, 'iMo':2, 'pMo':3}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = strain
        atypes = ['BL6', 'DBA']
        ahash = {'C57BL6/J':0, 'DBA2/J':1}
    if (tn >= 5 and tn < 9):
        atype = strain
        atype = [atype[i] if aval[i]==(tn-5) else None for i in range(len(atype))]
        atypes = ['BL6', 'DBA']
        ahash = {'C57BL6/J':0, 'DBA2/J':1}
    if (tn == 9):
        atypes = ['IM', 'AM']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSajti2020II = getSajti2020II


def getShemer2018(self, tn=1):
    self.prepareData("MACV97.2")
    atype = self.h.getSurvName("c Title")
    atype = [str(k).split("_")[4] if len(str(k).split("_")) > 4 else None
            for k in atype]
    atypes = ['M0', 'M1']
    ahash = {'LPS':1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShemer2018 = getShemer2018


def getLink2018(self, tn=1):
    self.prepareData("MACV98")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_rep.*", "", str(k)) for k in atype]
    group = [re.sub("BMDM_RNA_polyA_", "", str(k)) for k in atype]
    strain = self.h.getSurvName("c strain")
    ahash = {'C57BL/6J':1, 'SPRET/EiJ':2, 'BALB/cJ':0, 'NOD/ShiLtJ':3, 'PWK/PhJ':4}
    rval = [ahash[i] if i in ahash else None for i in strain]
    atype = self.h.getSurvName("c ligands in culture")
    atypes = ['M0', 'M1']
    ahash = {'no treatment':0, 'KLA 6h':1, 'KLA 1h':1}
    if (tn == 2):
        atype = group
        ahash = {'BALB_notx':0, 'BALB_KLA_1h':1}
    if (tn == 3):
        atype = group
        ahash = {'C57_notx_6h':0, 'C57_KLA_6h':1}
    if (tn == 4):
        atype = group
        atypes = ['Bl6', 'Bl6t', 'Bc', 'Bct']
        ahash = {'C57_notx_6h':0, 'C57_KLA_6h':1,
                'BALB_notx':2, 'BALB_KLA_1h':3}
    if (tn == 5):
        atype = group
        atypes = ['Bc', 'Bl6']
        ahash = {'C57_notx_6h':1, 'C57_KLA_6h':1,
                'BALB_notx':0, 'BALB_KLA_1h':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLink2018 = getLink2018


def getMunck2019(self, tn=1):
    self.prepareData("MACV99")
    strain = self.h.getSurvName("c strain")
    ahash = {'C57BL/6':1, 'BALB/c':0}
    rval = [ahash[i] if i in ahash else None for i in strain]
    atype = self.h.getSurvName("c infection status")
    atypes = ['M0', 'M1']
    ahash = {'infected':1, 'control':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = strain
        atypes = ['Bc', 'B6']
        ahash = {'C57BL/6':1, 'BALB/c':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMunck2019 = getMunck2019


def getElderman2018(self, tn=1):
    self.prepareData("MACV100")
    strain = self.h.getSurvName("c strain")
    ahash = {'Balb/c':0, 'C57bl/6':1}
    rval = [ahash[i] if i in ahash else None for i in strain]
    atype = self.h.getSurvName("c tissue")
    atypes = ['colon', 'ileum']
    ahash = {'distal ileum':1, 'proximal colon':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = strain
        atypes = ['Bc', 'B6']
        ahash = {'Balb/c':0, 'C57bl/6':1}
    if (tn == 5):
        atype = strain
        atype = [atype[i] if aval[i] == 0 else None for i in range(len(atype))]
        atypes = ['Bc', 'B6']
        ahash = {'Balb/c':0, 'C57bl/6':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getElderman2018 = getElderman2018


def getHowes2016(self, tn=1):
    self.prepareData("MACV101")
    strain = self.h.getSurvName("c strain")
    ahash = {'C57BL/6':1, 'BALB/c':0}
    rval = [ahash[i] if i in ahash else None for i in strain]
    gtype = self.h.getSurvName("c genotype/variation")
    ahash = {'WT':0, 'IFNabRKO':1}
    gval = [ahash[i] if i in ahash else None for i in gtype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1']
    ahash = {'HkBps':1, 'media':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 4):
        atype = strain
        atypes = ['Bc', 'B6']
        ahash = {'C57BL/6':1, 'BALB/c':0}
    if (tn == 5):
        atype = strain
        atype = [atype[i] if aval[i] == 0 and gval[i] == 0 \
                else None for i in range(len(atype))]
        atypes = ['Bc', 'B6']
        ahash = {'C57BL/6':1, 'BALB/c':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHowes2016 = getHowes2016


def getJochems2018(self, tn=1):
    self.prepareData("MACV102")
    atype = self.h.getSurvName("c src1")
    atypes = ['TIV', 'LAIV']
    ahash = {'LAIV_nasal cells_D-4':1,
            'TIV_nasal cells_D-4':0,
            'LAIV_nasal cells_D2':1,
            'TIV_nasal cells_D2':0,
            'LAIV_nasal cells_D9':1,
            'TIV_nasal cells_D9':0}
    if (tn == 2):
        atypes = ['TIV', 'LAIV']
        ahash = {'LAIV_nasal cells_D2':1,
                'TIV_nasal cells_D2':0}
    if (tn == 3):
        atype = self.h.getSurvName("c carriage status")
        atypes = ["NEG", "POS"]
        ahash = {}
    if (tn == 4):
        ahash = {'LAIV_nasal cells_D2':1,
                'TIV_nasal cells_D2':0}
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.h.getSurvName("c carriage status")
        atype = [atype[i] if aval[i] is not None else None for i in range(len(atype))]
        atypes = ["NEG", "POS"]
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJochems2018 = getJochems2018


def getZhai2015(self, tn=1):
    self.prepareData("MACV103")
    atype = self.h.getSurvName("c time point")
    atypes = ['C', 'I']
    ahash = {'Baseline':0, 'Spring':0, 'Day0':1, 'Day6':1, 'Day21':1,
            'Day2':1, 'Day4':1}
    if (tn == 2):
        ahash = {'Baseline':0, 'Day0':1}
    if (tn == 3):
        atypes = ['CV', 'AV']
        ahash = {'Day0':1, 'Day21':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhai2015 = getZhai2015


def getMitchell2013(self, tn=1):
    self.prepareData("MACV104")
    time = self.h.getSurvName("c timepoint")
    atype = [re.sub("h.*", "", str(k)) for k in time]
    ahash = {'0':0, '18':3, '12':2, '6':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection code")
    atypes = ['C', 'I']
    ahash = {'BatSRBD':1, 'icSARS':1, 'Mock':0, 'dORF6':1, 'H1N1':1, 'mock':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c infection code")
        ahash = {'BatSRBD':1, 'icSARS':1, 'Mock':0, 'mock':0}
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['C' if aval[i] == 1 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName("c infection code")
        ahash = {'H1N1':1, 'Mock':0, 'mock':0}
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['C' if aval[i] == 1 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMitchell2013 = getMitchell2013


def getWallihan2018(self, tn=1):
    self.prepareData("MACV105")
    atype = self.h.getSurvName('c viral organism')
    ahash = {'Coronavirus':1, 'RSV+Coronavirus':1,
            'Coronavirus+Bocavirus+Adenovirus':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c race')
    ahash = {'Other':3, 'Black or African American':1,
             'White':0, 'Asian':2, 'Unknown':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c condition")
    atypes = ['C', 'I', 'CoV']
    ahash = {'Pneumonia':1, 'Healthy Control':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] != 1 else 'CoV' for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = ['W', 'B', 'A']
        ahash = {0:0, 1:1, 2:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWallihan2018 = getWallihan2018


def getPfaender2020(self, tn=1):
    self.prepareData("MACV106.2")
    atype = self.h.getSurvName("c conditional ly6e knock-out")
    ahash = {'wt':0, 'ko':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c tissue")
    ahash = {'spleen':0, 'liver':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c inoculation")
    atypes = ['PBS', 'MHV']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 and gval[i] == 0 
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 0 and gval[i] == 1 
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if rval[i] == 1 and gval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if rval[i] == 1 and gval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPfaender2020 = getPfaender2020


def getJones2019(self, tn=1):
    self.prepareData("MACV107")
    atype = self.h.getSurvName("c visit")
    ahash = {'AV':0, 'CV':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c src1")
    ahash = {'NMS':0, 'PBMC':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c virus positive at av (1=yes, 0=no, 9=not measured)")
    atypes = ['0', '1']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c human coronavirus at av (1=yes, 0=no, 9=not measured)")
        atype = [atype[i] if rval[i] == 0 and gval[i] == 0
                else None for i in range(len(atype))]
    if (tn >= 3):
        atype = self.h.getSurvName("c visit")
        ahash = {'AV':1, 'CV':0}
        atypes = ['CV', 'AV']
    if (tn == 4):
        atype = [atype[i] if rval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if rval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJones2019 = getJones2019


def getBermejoMartin2010(self, tn=1):
    self.prepareData("MACV108")
    mtype = self.h.getSurvName("c mechanical ventilation")
    atype = self.h.getSurvName("c disease phase")
    mtype = [ " ".join([str(atype[i]), str(mtype[i])])
            for i in range(len(atype))]
    atypes = ['C', 'ENMV', 'LNMV', 'EMV', 'LMV']
    ahash = {'late period no':2, 'early period no':1,
            'early period yes':3, 'late period yes':4, 'control ':0}
    mval = [ahash[i] if i in ahash else None for i in mtype]
    ptype = self.h.getSurvName("c patient")
    ph = {'7':1, '8':1, '11':1, '17':1, '19':1}
    rval = [ph[i] if i in ph else None for i in ptype]
    atype = self.h.getSurvName("c disease phase");
    atypes = ['C', 'E', 'L']
    ahash = {'late period':2, 'early period':1, 'control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c mechanical ventilation");
        atypes = ['NMV', 'MV']
        ahash = {'no':0, 'yes':1}
    if (tn == 3 or tn == 4):
        atype = mtype
        atypes = ['C', 'ENMV', 'LNMV', 'EMV', 'LMV']
        ahash = {'late period no':2, 'early period no':1,
                'early period yes':3, 'late period yes':4, 'control ':0}
    if (tn == 4):
        atype = mtype
        atypes = ['C', 'ENMV', 'LNMV', 'EMV', 'LMV', 'LMVD']
        atype = [mtype[i] if mval[i] != 4 or rval[i] != 1 else "LMVD"
                for i in range(len(atype))]
        ahash = {'late period no':2, 'early period no':1,
                'early period yes':3, 'late period yes':4, 'control ':0}
    if (tn == 5):
        atype = mtype
        atypes = ['LMV', 'LMVD']
        atype = [mtype[i] if mval[i] != 4 or rval[i] != 1 else "LMVD"
                for i in range(len(atype))]
        ahash = {'late period yes':0}
    if (tn == 6):
        atype = mtype
        atypes = ['C', 'A', 'D']
        atype = [mtype[i] if mval[i] != 4 or rval[i] != 1 else "D"
                for i in range(len(atype))]
        ahash = {'late period no':1, 'early period no':1,
                'early period yes':1, 'late period yes':1, 'control ':0}
    if (tn == 7):
        atype = mtype
        atypes = ['CV', 'AV']
        atype = [mtype[i] if mval[i] != 4 or rval[i] != 1 else "D"
                for i in range(len(atype))]
        ahash = {'late period no':0, 'early period no':1,
                'early period yes':1, 'late period yes':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBermejoMartin2010 = getBermejoMartin2010


def getCameron2007(self, tn=1):
    self.prepareData("MACV109")
    atype = self.h.getSurvName("c Status");
    atypes = ['HC', 'C', 'Pre', 'Post']
    ahash = {'pre-pO2 nadir':2, 'post-pO2 nadir':3, 'healthy control':0,
            'convalescent':1}
    if (tn == 2):
        atypes = ['HC', 'I']
        ahash = {'pre-pO2 nadir':1, 'post-pO2 nadir':1, 'convalescent':1,
                'healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCameron2007 = getCameron2007


def getJosset2014(self, tn=1):
    self.prepareData("MACV110")
    atype = self.h.getSurvName("c virus");
    atypes = ['C', 'InfA', 'CoV'];
    ahash = {'MA15':2, 'MOCK':0, 'PR8':1}
    if (tn == 2):
        atypes = ['C', 'InfA']
        ahash = {'MOCK':0, 'PR8':1}
    if (tn == 3):
        atypes = ['C', 'CoV'];
        ahash = {'MA15':1, 'MOCK':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJosset2014 = getJosset2014


def getPrice2020(self, tn=1):
    self.prepareData("MACV111")
    mtype = self.h.getSurvName("c tissue")
    atype = self.h.getSurvName("c infection condtion")
    atype = [ " ".join([str(atype[i]), str(mtype[i])]) for i in
            range(len(atype))]
    atypes = ['S', 'SI', 'L', 'LI'];
    ahash = {'Mock Spleen':0, 'Mock Liver':2,
            'MA-EBOV infected Spleen':1, 'MA-EBOV infected Liver':3}
    if (tn == 2):
        atypes = ['C', 'I'];
        ahash = {'Mock Liver':0, 'MA-EBOV infected Liver':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPrice2020 = getPrice2020


def getReynard2019(self, tn=1):
    self.prepareData("MACV112")
    atype = self.h.getSurvName("c group")
    atypes = ['HC', 'SR', 'VR', 'F'];
    ahash = {'Fatalities':3, 'Healthy controls':0,
            'Survivors in recovery phase':1, 'Viremic survivors':2}
    if (tn == 2):
        atypes = ['HC', 'R', 'I'];
        ahash = {'Fatalities':2, 'Healthy controls':0,
                'Survivors in recovery phase':1, 'Viremic survivors':2}
    if (tn == 3):
        atypes = ['C', 'I'];
        ahash = {'Fatalities':1, 'Healthy controls':0,
                'Survivors in recovery phase':0, 'Viremic survivors':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getReynard2019 = getReynard2019


def getDunning2018(self, tn=1):
    self.prepareData("MACV113")
    atype = self.h.getSurvName("c t1severity")
    atypes = ['HC', '1', '2', '3'];
    ahash = {'HC':0, '1':1, '2':2, '3':3}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'I']
        ahash = {'1':1, '2':1, '3':1}
    if (tn == 3):
        atype = self.h.getSurvName("c ethnicity")
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'White':0, 'Other':3, 'Black':1, 'Asian':2}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDunning2018 = getDunning2018


def getBerdal2011(self, tn=1):
    self.prepareData("MACV114")
    atype = self.h.getSurvName("c control")
    expt = ["P" if k == "" else "" for k in atype]
    atype = self.h.getSurvName("c patient")
    ctrl = ["C" if k == "" else "" for k in atype]
    atype = [ str(expt[i]) + str(ctrl[i]) for i in range(len(atype)) ]
    atypes = ['C', 'P']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerdal2011 = getBerdal2011


def getXie2018(self, tn=1):
    self.prepareData("MACV115")
    atype = self.h.getSurvName("c disease interval")
    atypes = ['Con', 'Cr']
    ahash = {'convalescent':0, 'crisis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXie2018 = getXie2018


def getFriesenhagen2012(self, tn=1):
    self.prepareData("MACV116")
    atype = self.h.getSurvName("c desc")
    atype = [re.sub("Gene.*from ", "", str(k)) for k in atype]
    atype = [re.sub(" mac.*", "", str(k)) for k in atype]
    atypes = ['C', 'PR8', 'H5N1', 'FPV']
    ahash = {'FPV-infected':3, 'H5N1-infected':2, 'control':0, 'PR8-infected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFriesenhagen2012 = getFriesenhagen2012


def getGuan2018(self, tn=1):
    self.prepareData("MACV117")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(", .*", "", str(k)) for k in atype]
    atype = [ re.split(" ", str(k))[1] if len(re.split(" ", str(k))) > 1
                            else None for k in atype]
    ahash = {'2':2, 'control':0, '3':3, '8':8, '7':7, '4':4,
            '9':9, '6':6, '10':10}
    pval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c number of day post-infection")
    atype = [re.sub("NA", "0", str(k)) for k in atype]
    phash = {}
    for i in range(2, len(atype)):
        if pval[i] not in phash:
            phash[pval[i]] = []
        phash[pval[i]].append([i, int(atype[i])])
    before = []
    after = []
    for i in phash.keys():
        if i == 0:
            before += [k[0] for k in phash[i]]
            after += [k[0] for k in phash[i]]
        else:
            ll = sorted(phash[i], key = lambda k: k[1])
            before.append(ll[0][0])
            before.append(ll[1][0])
            after.append(ll[-1][0])
    beforehash = set(before)
    afterhash = set(after)
    ahash = {'2':0, '3':1, '8':1, '7':0, '4':1,
            '9':0, '6':0, '10':1}
    gender = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c subject category")
    atypes = ['C', 'P']
    ahash = {'Patient':1, 'Control':0}
    if (tn == 2):
        atype = [atype[i] if gender[i] == 0 or pval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        #Categories: 1) Control; 2) Mild - non-invasive ventilation (#4,#5);
        #3) Moderate- MV, but discharged within 2 mo (#1, 2, 7, 8);
        #4) Severe-  MV+prolonged hospitalization (#6 and 9);
        #5) Death after MV + ECMO = (#3 and 10)
        atype = pval
        atypes = ['C', 'Mi', 'Mo', 'S', 'D']
        ahash = {0:0, 4:1, 5:1, 1:2, 2:2, 7:2, 8:2, 6:3, 9:3, 3:4, 10:4}
        #atype = [atype[i] if gender[i] == 1 or pval[i] == 0
        #        else None for i in range(len(atype))]
    if (tn == 4):
        atype = pval
        atypes = ['C', 'Mi', 'MV', 'D']
        ahash = {0:0, 4:1, 5:1, 1:2, 2:2, 7:2, 8:2, 6:2, 9:2, 3:3, 10:3}
        atype = [atype[i] if i in beforehash
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = pval
        atypes = ['C', 'Mi', 'Mo', 'S', 'D']
        ahash = {0:0, 4:1, 5:1, 1:2, 2:2, 7:2, 8:2, 6:3, 9:3, 3:4, 10:4}
        atype = [atype[i] if i in afterhash
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuan2018 = getGuan2018


def getKongchanagul2011(self, tn=1):
    self.prepareData("MACV118")
    atype = self.h.getSurvName("c disease state")
    atypes = ['C', 'H5N1']
    ahash = {'H5N1 influenza':1, 'normal':0}
    if (tn == 2):
        atype = self.h.getSurvName("c time of death")
        atypes = ['C', '6', '17']
        ahash = {'Day 17 of illness':2, 'Day 6 of illness':1, 'N/A':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKongchanagul2011 = getKongchanagul2011


def getHo2013(self, tn=1):
    self.prepareData("MACV119")
    atype = self.h.getSurvName("c gender")
    atypes = ['F', 'M']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c age")
        atypes = ['C', '6', '17']
        ahash = {'Day 17 of illness':2, 'Day 6 of illness':1, 'N/A':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHo2013 = getHo2013


def getKuparinen2013(self, tn=1):
    self.prepareData("MACV120")
    atype = self.h.getSurvName("c frailty index")
    atypes = ['nF', 'pF', 'F']
    ahash = {'non-frail':0, 'pre-frail':1, 'frail':2}
    if (tn == 2):
        atype = self.h.getSurvName("c age")
        aha = {'90':1, 'nonagenarian':1}
        age = ['O' if atype[i] in aha else 'Y' for i in range(len(atype))]
        atype = self.h.getSurvName("c Sex")
        ahash = {'Female':'F', 'female':'F', 'Male':'M', 'male':'M'}
        sex = [ahash[k] if k in ahash else None for k in atype]
        atype = [ str(age[i]) + str(sex[i]) for i in range(len(atype)) ]
        atypes = ['YF', 'YM', 'OF', 'OM']
        ahash = {}
    if (tn == 3):
        atype = self.h.getSurvName("c age")
        aha = {'90':1, 'nonagenarian':1}
        atype = ['O' if atype[i] in aha else 'Y' for i in range(len(atype))]
        atypes = ['Y', 'O']
        ahash = {}
    if (tn == 4):
        atype = self.h.getSurvName("c cmv serostatus")
        atypes = ['neg', 'pos']
        ahash = {'pos.':1, 'neg.':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKuparinen2013 = getKuparinen2013


def getWagstaffe2020(self, tn=1):
    self.prepareData("MACV121")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("^._", "", str(k)) for k in atype]
    atypes = ['U-', 'U+', 'E-', 'E+']
    ahash = {'Med_CD14-':0, 'EBOV_CD14+':3, 'Med_CD14+':1, 'EBOV_CD14-':2}
    if (tn == 2):
        atypes = ['U', 'E']
        ahash = {'Med_CD14-':0, 'EBOV_CD14+':1, 'Med_CD14+':0, 'EBOV_CD14-':1}
    if (tn == 3):
        atypes = ['U', 'E']
        ahash = {'Med_CD14-':0, 'EBOV_CD14-':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWagstaffe2020 = getWagstaffe2020


def getPeng2016I(self, tn=1):
    self.prepareData("MACV122")
    atype = self.h.getSurvName("c disease state")
    ahash = {'Chronic Obstructive Lung Disease':2,
            'Interstitial lung disease':1,
            'Control':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Sex")
    atypes = ['F', 'M']
    ahash = {'1-Male':1, '2-Female':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = ['C', 'ILD', 'COPD']
        ahash = {0:0, 1:1, 2:2}
    if (tn == 4):
        atype = self.h.getSurvName("c gold stage")
        atypes = ['R', 'Mi', 'Mo', 'S', 'VS']
        ahash = {'4-Very Severe COPD':4, '0-At Risk':0, '2-Moderate COPD':2,
                '1-Mild COPD':1, '3-Severe COPD':3}
    if (tn == 5):
        atype = tval
        atypes = ['C', 'COPD']
        ahash = {0:0, 2:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeng2016I = getPeng2016I


def getPeng2016II(self, tn=1):
    self.prepareData("MACV122.2")
    atype = self.h.getSurvName("c disease state")
    ahash = {'Chronic Obstructive Lung Disease':2,
            'Interstitial lung disease':1,
            'Control':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Sex")
    atypes = ['F', 'M']
    ahash = {'1-Male':1, '2-Female':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = ['C', 'ILD', 'COPD']
        ahash = {0:0, 1:1, 2:2}
    if (tn == 4):
        atype = self.h.getSurvName("c gold stage")
        atypes = ['R', 'Mi', 'Mo', 'S', 'VS']
        ahash = {'4-Very Severe COPD':4, '0-At Risk':0, '2-Moderate COPD':2,
                '1-Mild COPD':1, '3-Severe COPD':3}
    if (tn == 5):
        atype = tval
        atypes = ['C', 'COPD']
        ahash = {0:0, 2:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeng2016II = getPeng2016II


def getBosse2012(self, tn=1):
    self.prepareData("MACV123")
    atype = self.h.getSurvName("c src1")
    atypes = ['Lung']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBosse2012 = getBosse2012


def getMenachery2017(self, tn=1):
    self.prepareData("COV1.3")
    time1 = self.h.getSurvName("c time")
    time2 = self.h.getSurvName("c time-post-infection")
    time3 = self.h.getSurvName("c time-hrs post infection")
    time4 = self.h.getSurvName("c time point")
    atype = [ "".join([str(k) for k in
                           [time1[i], time2[i], time3[i], time4[i]]])
                                    for i in range(len(time1))]
    atype = [re.sub("[h ].*", "", k) for k in atype]
    ahash = {'00':'0'}
    tval = [ahash[i] if i in ahash else i for i in atype]
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(",.*", "", str(k)) for k in atype]
    ahash = {'Primary human fibroblasts':0,
            'Primary human airway epithelial cells':1,
            'Primary human microvascular endothelial cells':2,
            'Primary human dendritic cells':3}
    rval = [ahash[i] if i in ahash else None for i in atype]
    v1 = self.h.getSurvName("c virus")
    v2 = self.h.getSurvName("c virus infection")
    v3 = self.h.getSurvName("c infection")
    v4 = self.h.getSurvName("c infected with")
    atype = [ "".join([str(k) for k in [v1[i], v2[i], v3[i], v4[i]]])
            for i in range(len(atype))]
    atypes = ['C', 'I']
    ahash = { 'Mock':0, 'icMERS':1, 'RFP-MERS':1, 'mock':0, 'd4B-MERS':1,
            'MOCK':0, 'Mockulum':0, 'dNSP16-MERS':1, 'd3-5-MERS':1,
            'MERS-coronavirus (icMERS)':1}
    if (tn >= 2 and tn <= 5):
        atype = [atype[i] if rval[i] == (tn - 2) else None
                for i in range(len(atype))]
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['C' if aval[i] == 1 and tval[i] == '0' else atype[i]
                for i in range(len(atype))]
    if (tn == 6):
        atype = rval
        atypes = ['FI', 'AE', 'ME', 'DC']
        ahash = {0:0, 1:1, 2:2, 3:3}
    if (tn == 7):
        atype = self.h.getSurvName("c Series")
        atypes = ['FI', 'AE', 'ME', 'DC']
        ahash = {'GSE100496':0, 'GSE100504':1, 'GSE100509':2, 'GSE79172':3}
    if (tn == 8):
        ctype = self.h.getSurvName("c cell type")
        atypes = ['C', 'I']
        ahash = { 'Mock':0, 'mock':0,
                'MOCK':0, 'Mockulum':0, 'EBOV-WT':1}
        atype = [atype[i] if ctype[i] == 'Immortalized Human Hepatocytes (IHH)'
                else None for i in range(len(atype))]
    self.rval = rval
    self.tval = tval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMenachery2017 = getMenachery2017


def getYoshikawa2010(self, tn=1):
    self.prepareData("COV2")
    atype = self.h.getSurvName("c time")
    ahash = {'48 hours post infection':3,
            '24 hours post infection':2,
             '12 hours post infection':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['C', 'CoV', 'DoV']
    ahash = {'Mock-infected':0,
            'SARS-CoV-infected (MOI=0.1)':1,
            'DOHV-infected (MOI=0.1)':2}
    if (tn == 2):
        atypes = ['C', 'CoV']
        ahash = {'Mock-infected':0,
                'SARS-CoV-infected (MOI=0.1)':1}
        atype = [atype[i] if tval[i] == 3 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYoshikawa2010 = getYoshikawa2010


def getKelvin2014(self, tn=1):
    self.prepareData("COV3")
    atype = self.h.getSurvName("c time point")
    ahash = {'Day 3':3, 'Day 2':2, 'Day 1':1, 'Day 5':5,
            'Day 0':0, 'Day 28':28, 'Day 14':14}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'SARS-CoV (TOR-2)':1, 'Uninfected':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] is not None and 
                (tval[i] == 0 or tval[i] <= 3) else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKelvin2014 = getKelvin2014


def getDeDiego2011(self, tn=1):
    self.prepareData("COV4")
    atype = self.h.getSurvName("c hours post infection")
    ahash = {'15 hpi':15, '24 hpi':24, '65 hpi':65, '7 hpi':7, '0 hpi':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'SARS CoV':1, 'SARS CoV DeltaE':1, 'Mock infected':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 or tval[i] == 65 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeDiego2011 = getDeDiego2011


def getSims2013(self, tn=1):
    self.prepareData("COV5")
    atype = self.h.getSurvName("c time")
    ahash = {'30h':30, '36h':36, '7h':7, '60h':60, '48h':48, '0h':0,
            '12h':12, '72h':72, '24h':24, '54h':54, '3h':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['U', 'CoV']
    ahash = {'mock infected':0,
            'SARS CoV Urbani infected':1,
            'SARS delta ORF6 infected':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] is not None and 
                (tval[i] == 0 or tval[i] >= 64) else None
                for i in range(len(atype))]
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['U' if aval[i] == 1 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] is not None and tval[i] > 40 else None
                for i in range(len(atype))]
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['U' if aval[i] == 1 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    self.tval = tval
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSims2013 = getSims2013


def getKatze2012(self, tn=1):
    self.prepareData("COV6")
    atype = self.h.getSurvName("c sample collection time post virus infection")
    ahash = {'36h':36, '72h':72, '7h':7, '60h':60, '54h':54,
            '30h':30, '12h':12, '24h':24, '0h':0, '48h':48}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'SARS CoV Urbani infected':1,
            'SARS Bat SRBD infected':1,
            'mock infected':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] is not None and 
                (tval[i] == 0 or tval[i] >= 18) else None
                for i in range(len(atype))]
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['U' if aval[i] == 1 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKatze2012 = getKatze2012


def getJosset2013(self, tn=1):
    self.prepareData("COV7")
    atype = self.h.getSurvName("c time point")
    ahash = {'18 hpi':18, '12 hpi':12, '3 hpi':3,
            '24 hpi':24, '7 hpi':7, '0 hpi':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infected with")
    atypes = ['U', 'CoV']
    ahash = {'Human Coronavirus EMC 2012 (HCoV-EMC)':1, 'Mock':0}
    if (tn == 2):
        atype = [atype[i] if (tval[i] == 12 and atype[i] == 'Mock') 
                or tval[i] == 18 else None
                for i in range(len(atype))]
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = ['U' if aval[i] == 1 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJosset2013 = getJosset2013


def getKatze2014(self, tn=1):
    self.prepareData("COV8")
    atype = self.h.getSurvName("c time")
    ahash = {'d7':7, 'd2':2, 'd4':4, 'd1':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'Mock':0, 'SARS MA15':1, 'SARS CoV':1, 'SARS BatSRBD mutant':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] >= 4 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKatze2014 = getKatze2014


def getJimenezGuardeno2014(self, tn=1):
    self.prepareData("COV9")
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'Mock':0, 'SARS-CoV-wt':1, 'SARS-CoV-mutPBM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJimenezGuardeno2014 = getJimenezGuardeno2014


def getSelinger2014(self, tn=1):
    self.prepareData("COV10")
    atype = self.h.getSurvName("c timepoint")
    ahash = {'18h':18, '7h':7, '0h':0, '12h':12, '24h':24, '3h':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'Mock  Infected':0, 'LoCoV':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 18 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSelinger2014 = getSelinger2014


def getTotura2015(self, tn=1):
    self.prepareData("COV11")
    atype = self.h.getSurvName("c time")
    ahash = {'4 dpi':4, '7 dpi':7, '2 dpi':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'CoV']
    ahash = {'MA15 virus':1, 'mockulum':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTotura2015 = getTotura2015


def getFerris2017(self, tn=1):
    self.prepareData("COV12")
    atype = self.h.getSurvName("c src1")
    atypes = ['I']
    ahash = {'lung tissue, 4 days post SARS-CoV infection':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFerris2017 = getFerris2017


def getMoodley2016(self, tn=1):
    self.prepareData("COV13")
    v1 = self.h.getSurvName("c jak inhinitor")
    v2 = self.h.getSurvName("c jak inhibitor")
    atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                    for i in range(len(v1))]
    atypes = ['U', 'T']
    ahash = {'IL2 only':2, 'PAN':1, 'JAK1/2':1, 'JAK1':1, 'Untreated':0, 'JAK3':1}
    if (tn == 2):
        ahash = {'PAN':1, 'JAK1/2':1, 'Untreated':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMoodley2016 = getMoodley2016


def getPardanani2013(self, tn=1):
    self.prepareData("COV14")
    atype = self.h.getSurvName("c pre or post-treatment")
    atypes = ['U', 'T']
    ahash = {'post':1, 'pre':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c anemia response to treatment")
        atypes = ['R', 'NR']
        ahash = {'responder':0, 'Non-responder':1}
        atype = [atype[i] if aval[i] == 1 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPardanani2013 = getPardanani2013


def getKuo2019(self, tn=1):
    self.prepareData("COV15")
    atype = self.h.getSurvName("c treatment")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['U', 'F', 'T']
    ahash = {'none':0, 'TNF':2, 'Fib':1}
    if (tn == 2):
        atype = self.h.getSurvName("c cultured with")
        atypes = ['U', 'F', 'T', 'TF']
        ahash = {'none (alone)':0,
                'fibroblasts':1,
                'tumor necrosis factor (TNF)':2,
                'tumor necrosis factor (TNF) +fibroblasts':3}
    if (tn == 3):
        atype = self.h.getSurvName("c treatment")
        atypes = ['U', 'T', 'T+Q', 'T+T']
        ahash = {'none':0, 'TNF':1, 'TNF + Hydrox':2, 'TNF + tofa':3}
    if (tn == 4):
        atype = self.h.getSurvName("c Title")
        atypes = ['MT', 'MTaIL6', 'MTF', 'MTFaIL6']
        ahash = {"P2_MT_aIL6":1, "P2_MTF_aIL6":3, "P1_MTF_aIL6":3,
                "P1_MT_aIL6":1, 'P1_MT':0, 'P1_MTF':2,
                'P2_MT':0, 'P2_MTF':2, 'P3_MT':0, 'P3_MTF':2, 
                'P4_MT':0, 'P4_MTF':2}
    if (tn == 5):
        atype = self.h.getSurvName("c Title")
        atypes = ['MT', 'MTFaIL6']
        ahash = {"P2_MTF_aIL6":1, "P1_MTF_aIL6":1,
                'P1_MT':0, 'P2_MT':0, 'P3_MT':0, 'P4_MT':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKuo2019 = getKuo2019


def getSmyth2016(self, tn=1):
    self.prepareData("COV16")
    atype = self.h.getSurvName("c treated with")
    atypes = ['U', 'Q', 'S', 'S+Q']
    ahash = {'none (untreated)':0,
            '20 \xc2\xb5M hydroxychloroquine (HCQ)':1,
            'group A streptococcus':2,
            'group A streptococcus and hydroxychloroquine':3}
    if (tn == 2):
        atypes = ['S', 'S+Q']
        ahash = {'group A streptococcus':0,
                'group A streptococcus and hydroxychloroquine':1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmyth2016 = getSmyth2016


def getWhyte2007(self, tn=1):
    self.prepareData("COV17")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub("rol, .*", "rol", str(k)) for k in atype]
    atype = [re.sub(".*h, ", "", str(k)) for k in atype]
    atypes = ['U', 'R']
    ahash = {'25 microM resveratrol':1, 'ethanol control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWhyte2007 = getWhyte2007


def getBaty2006(self, tn=1):
    self.prepareData("COV18")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("HC .. ", "", str(k)) for k in atype]
    atype = [re.sub("^. ", "", str(k)) for k in atype]
    atypes = ['C', 'W', 'A', 'G']
    ahash = {'12 wine':1, '12 alcohol':2, '12 grape.juice':3,
            '12 water':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBaty2006 = getBaty2006


def getAuerbach2014(self, tn=1, tv=0, dg=''):
    self.prepareData("COV19")
    atype = self.h.getSurvName("c tissue")
    ahash = {'Liver':0, 'Kidney':1, 'Bone marrow':2, 'Intestine':3,
            'Brain':4, 'Heart':5, 'Spleen':6, 'Skeletal muscle':7,
            'Primary rat hepatocytes':8}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c time")
    ahash = { '0 d':0, '0.25 d':0.25, '0.67 d':0.67, '1 d':1,
            '3 d':3, '4 d':4, '5 d':5, '6 d':6, '7 d':7,
            '14 d':14, '28 d':28, '30 d':30, '91 d':91}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = ['C' if dval[i] is None or dval[i] <= 3
            else 'T' for i in range(len(atype))]
    atypes = ['C', 'T']
    if (tn == 2):
        atype = self.h.getSurvName("c compound")
        atypes = ['C', 'D']
        ahash = {'':0, 'Captopril':1}
        atype = [atype[i] if tval[i] == 1 else None for i in range(len(atype))]
    if (tn == 3):
        ahash = {'Benazepril':1, 'Captopril':1, 'Enalapril':1,
                'Quinapril':1, 'Ramipril':1, 'Lisinopril':1, 'Flutamide':1}
        atype = self.h.getSurvName("c compound")
        ahash = {'':0, 'Captopril':1}
        pval = [ahash[i] if i in ahash else None for i in atype]
        atype = dval
        atype = [0 if pval[i] == 0 else dval[i]
                for i in range(len(atype))]
        if (dg == ''):
            atype = [atype[i] if pval[i] is not None and tval[i] == 1 else None
                    for i in range(len(atype))]
        else:
            atype = [atype[i] if pval[i] is not None and tval[i] == 1 and
                    (dval[i] == 0 or dval[i] == dg)
                    else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    if (tn == 4):
        atype = self.h.getSurvName("c compound")
        ahash = {'':0, 'Chlorpromazine':1}
        pval = [ahash[i] if i in ahash else None for i in atype]
        atype = dval
        atype = [0 if pval[i] == 0 else dval[i]
                for i in range(len(atype))]
        if (dg == ''):
            atype = [atype[i] if pval[i] is not None and tval[i] == 0 
                    else None for i in range(len(atype))]
        else:
            atype = [atype[i] if pval[i] is not None and tval[i] == 0 and
                    (dval[i] == 0 or dval[i] == dg)
                    else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        #atypes = [0, 0.25]
        ahash = {}
    if (tn == 5):
        atype = self.h.getSurvName("c compound")
        atypes = ['C', 'D']
        ahash = {'':0, dg:1}
        atype = [atype[i] if tval[i] == tv else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAuerbach2014 = getAuerbach2014


def getRuiz2016(self, tn=1):
    self.prepareData("COV20")
    atype = self.h.getSurvName("c treatment")
    atypes = ['U', 'O', 'O+C']
    ahash = {'basal condition, untreated':0,
            'combined treatment, oxaliplatin plus curcumin':2,
            'single treatment, oxaliplatin':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRuiz2016 = getRuiz2016


def getMeja2008(self, tn=1):
    self.prepareData("COV21")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(".*cells, ", "", str(k)) for k in atype]
    atypes = ['U', 'C', 'R', 'R+C']
    ahash = {'untreated, 18h':0,
            'ROS exposed, 18h':2,
            'ROS exposed, 4h':2,
            'ROS exposed, 1uM curcumin treated, 4h':3,
            '1uM curcumin treated, 18h':1,
            'untreated, 4h':0,
            'ROS exposed, 1uM curcumin treated, 18h':3,
            '1uM curcumin treated, 4h':1}
    if (tn == 2):
        atypes = ['U', 'C']
        ahash = {'untreated, 18h':0,
                '1uM curcumin treated, 18h':1,
                'untreated, 4h':0,
                '1uM curcumin treated, 4h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeja2008 = getMeja2008


def getGarland2019(self, tn=1):
    self.prepareData("COV22")
    v1 = self.h.getSurvName("c timepoint")
    v2 = self.h.getSurvName("c dosing regimen")
    atype = [ " ".join([str(k) for k in [v1[i], v2[i]]])
                            for i in range(len(v1))]
    atypes = ['B', '12', '13']
    ahash = {'13 weeks intermittent':2,
            '12 weeks continuous':1,
            'baseline continuous':0,
            'baseline intermittent':0,
            '12 weeks intermittent':1,
            '13 weeks continuous':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGarland2019 = getGarland2019


def getGuo2016(self, tn=1):
    self.prepareData("COV23")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'T']
    ahash = {'Aspirin treated for 48 hours':1, 'DMSO treated for 48 hours':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuo2016 = getGuo2016


def getPavuluri2014(self, tn=1):
    self.prepareData("COV24")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'T']
    ahash = {'Treated with 2.0 mM Aspirin':1, 'Untreated with Aspirin':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPavuluri2014 = getPavuluri2014


def getFallahi2013(self, tn=1):
    self.prepareData("COV25")
    atype = self.h.getSurvName("c group")
    atypes = ['N', 'S', 'R']
    ahash = {'aspirin sensitive':1, 'high normal':0, 'aspirin resistant':2}
    if (tn == 2):
        atypes = ['N', 'A']
        ahash = {'aspirin sensitive':1, 'high normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFallahi2013 = getFallahi2013


def getLewis2020(self, tn=1):
    self.prepareData("PLP113")
    atype = self.h.getSurvName("c Sex")
    atypes = ['F', 'M']
    ahash = {'Male':1, 'Female':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLewis2020 = getLewis2020


def getLv2017(self, tn=1):
    self.prepareData("COV26")
    atype = self.h.getSurvName("c group")
    ahash = {'treatment':1, 'control':0}
    atypes = ['C', 'T']
    if (tn == 2):
        atype = self.h.getSurvName("c perturbagen")
        ahash = {'DMSO':0, 'Resveratrol':1}
        atypes = ['C', 'R']
    if (tn == 3):
        atype = self.h.getSurvName("c perturbagen")
        ahash = {'DMSO':0, 'Hyodeoxycholic acid':1,
                'Ursodeoxycholic acid':1, 'Deoxycholic acid':1,
                'Chenodeoxycholic acid':1, 'Artemisinin':3, 'Resveratrol':2}
        atypes = ['C', 'B', 'R', 'A']
    if (tn == 4):
        atype = self.h.getSurvName("c perturbagen")
        ahash = {'DMSO':0, 'Hyodeoxycholic acid':1,
                'Ursodeoxycholic acid':1, 'Deoxycholic acid':1,
                'Chenodeoxycholic acid':1}
        atypes = ['C', 'B']
    if (tn == 5):
        atype = self.h.getSurvName("c perturbagen")
        ahash = {'DMSO':0, 'Artemisinin':1}
        atypes = ['C', 'A']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLv2017 = getLv2017


def getCMAP(self, tn=1):
    self.prepareData("COV27")
    atype = self.h.getSurvName("c type")
    ahash = {'treatment':1, 'control':0}
    atypes = ['C', 'T']
    if (tn == 2):
        atype = self.h.getSurvName("c name")
        ahash = {'null':0, 'chlorpromazine':1}
        atypes = ['C', 'CPZ']
    if (tn == 3):
        atype = self.h.getSurvName("c name")
        ahash = {'null':0, 'resveratrol':1}
        atypes = ['C', 'R']
    if (tn == 4):
        atype = self.h.getSurvName("c name")
        ahash = {'null':0, 'sirolimus':1}
        atypes = ['C', 'S']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCMAP = getCMAP


def getWoo2015(self, tn=1):
    self.prepareData("COV28")
    atype = self.h.getSurvName("c src1")
    atypes = ['OCI-Ly3', 'OCI-Ly7', 'U2932', 'HeLa']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c compound treated")
        ahash = {'DMSO':0, 'ESTRADIOL':1, 'Docetaxel': 2}
        atypes = ['C', 'E', 'D']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWoo2015 = getWoo2015


def getNiculescu2019(self, tn=1):
    self.prepareData("COV29")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['None', 'BP', 'SZA', 'MDD', 'SZ', 'PSYCH', 'PTSD', 'MOOD']
    ahash = {'':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNiculescu2019 = getNiculescu2019


def getBush2017(self, tn=1):
    self.prepareData("COV30")
    atype = self.h.getSurvName("c src1")
    atypes = ['B', 'L', 'U']
    ahash = {'BT20 cell line':0, 'LnCAP cell line':1, 'U87 cell line':2}
    if (tn == 2):
        atype = self.h.getSurvName("c drug")
        atypes = ['C', 'T']
        ahash = {'DMSO':0, 'untreated':0, 'Temsirolimus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBush2017 = getBush2017


def getPlauth2015(self, tn=1):
    self.prepareData("COV31")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'R']
    ahash = {'vehicle':0, '16 hours 50 \xc2\xb5M resveratrol':1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPlauth2015 = getPlauth2015


def getTomeCarneiro2013(self, tn=1):
    self.prepareData("COV32")
    atype = self.h.getSurvName("c time point")
    ahash = {'12 months':12, '6 months':6, 'day 0 (basal)':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c dietary group")
    atypes = ['C', 'G', 'R']
    ahash = {'placebo (A)':0,
            'grape extract (B)':1,
            'resveratrol-enriched grape extract (C)':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = ['C' if aval[i] != 0 and tval[i] == 0 else atype[i]
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 12 else None
                for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if tval[i] == 6 else None
                for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if tval[i] == 12 else None
                for i in range(len(atype))]
        atypes = ['C', 'R']
        ahash = {'placebo (A)':0,
                'resveratrol-enriched grape extract (C)':1}
    if (tn == 6):
        atype = [atype[i] if tval[i] == 12 else None
                for i in range(len(atype))]
        atypes = ['C', 'G']
        ahash = {'placebo (A)':0,
                'grape extract (B)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTomeCarneiro2013 = getTomeCarneiro2013


def getUribe2011(self, tn=1):
    self.prepareData("COV33")
    atype = self.h.getSurvName("c agent")
    atypes = ['C', 'R1', 'R2']
    ahash = {'0.03% ethanol (control)':0,
            'resveratrol_150mM':1,
            'resveratrol_250mM':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getUribe2011 = getUribe2011


def getEuba2015(self, tn=1):
    self.prepareData("COV34")
    atype = self.h.getSurvName("c agent")
    atypes = ['C', 'I']
    ahash = {'Haemophilus influenzae strain NTHi375':1, 'none':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEuba2015 = getEuba2015


def getDavis2019(self, tn=1):
    self.prepareData("COV35")
    atype = self.h.getSurvName("c tissue")
    atypes = ['L', 'BM', 'S', 'K', 'H']
    ahash = {'Liver':0, 'Bone Marrow':1, 'Skin':2, 'Kidney':3, 'Heart':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDavis2019 = getDavis2019


def getPodtelezhnikov2020(self, tn=1):
    self.prepareData("COV36")
    atype = self.h.getSurvName("c Sex")
    atypes = ['F', 'M']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'X', 'C', 'A', 'Q', 'CP']
        ahash = {'-control-':0, 'amoxicillin':1, 'caffeine':2,
                'acetaminophen':3, 'quinidine':4, 'captopril':5}
    if (tn == 3):
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'CP']
        ahash = {'-control-':0, 'captopril':1}
    if (tn == 4):
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'F']
        ahash = {'-control-':0, 'flutamide':1}
    if (tn == 5):
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'Q']
        ahash = {'-control-':0, 'quinidine':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPodtelezhnikov2020 = getPodtelezhnikov2020


def getMonks2018(self, tn=1):
    self.prepareData("GL15")
    atype = self.h.getSurvName("c Title")
    atype = [ re.split("_", str(k))[2] if len(re.split("_", str(k))) > 2
            else None for k in atype]
    conc = [int(re.sub("nM", "", k)) if k is not None else None for k in atype]
    atype = self.h.getSurvName("c tissue")
    ahash = {'Renal':0, 'CNS':1, 'Melanoma':2, 'Lung':3, 'Breast':4,
            'Ovarian':5, 'Colon':6, 'Leukemia':7, 'Prostate':8}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Time")
    ahash = {'2h':2, '6h':6, '24h':24}
    mval = [ahash[i] if i in ahash else None for i in atype]
    drug = self.h.getSurvName("c Drug")
    atype = self.h.getSurvName("c Time")
    atypes = ['2h', '6h', '24h']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c Drug")
        atypes = ['S', 'C', 'SN']
        ahash = {'sunitinib':2, 'sirolimus':0, 'cisplatin':1}
        atype = [atype[i] if tval[i] == 3 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 3 else None
                for i in range(len(atype))]
    if (tn == 4):
        atype = self.h.getSurvName("c Drug")
        atypes = ['sunitinib', 'dasatinib', 'sirolimus',
                'lapatinib', 'doxorubicin', 'sorafenib',
                'bortezomib', 'cisplatin', 'erlotinib']
        atype = [atype[i] if tval[i] == 3 and conc[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 5):
        atype = conc
        atype = [atype[i] if tval[i] == 6 and drug[i] == 'sirolimus' and
                mval[i] == 24 else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([k for k in atype if k is not None]))
    self.conc = conc
    self.mval = mval
    self.tval = tval
    self.drug = drug
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMonks2018 = getMonks2018


def getKornakiewicz2018(self, tn=1):
    self.prepareData("COV37")
    atype = self.h.getSurvName("c src1")
    atypes = ['SC', 'pSC']
    ahash = {'Human Kidney Cancer Stem Cells treated with everolimus':0,
            'Parental treated with everolimus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKornakiewicz2018 = getKornakiewicz2018


def getYu2018(self, tn=1):
    self.prepareData("COV38")
    atype = self.h.getSurvName("c agent")
    atypes = ['C', 'R']
    ahash = {'control':0, 'rapamycin':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYu2018 = getYu2018


def getSabine2010(self, tn=1):
    self.prepareData("COV39")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*\.P", "P", str(k)) for k in atype]
    atype = [re.sub("_.*", "", str(k)) for k in atype]
    atypes = ['Pre', 'Post']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSabine2010 = getSabine2010


def getDreschers2019(self, tn=1):
    self.prepareData("COV40")
    atype = self.h.getSurvName("c agent")
    atypes = ['1', '2', '3', '4', '5']
    ahash = {'IL10':0, 'IFN-\xce\xb3':1, 'M(IFN-gamma)':2,
            'M (IL-10)':3, 'M(IFN-gamma) + Rapamycin':4,
            'M (IL-10) + Rapamycin':5}
    if (tn == 2):
        atypes = ['C', 'C+R']
        ahash = {'M (IL-10)':0, 'M (IL-10) + Rapamycin':1}
    if (tn == 3):
        atypes = ['C', 'C+R']
        ahash = {'M(IFN-gamma)':0, 'M(IFN-gamma) + Rapamycin':1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDreschers2019 = getDreschers2019


def getCardoso2016(self, tn=1):
    self.prepareData("COV41")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'S', 'S+D']
    ahash = {'sunitinib and docetaxel combined treatment':2,
            'prior to any treatment':0,
            'sunitinib':1}
    if (tn == 2):
        atype = self.h.getSurvName("c neosu_response")
        atypes = ['good', 'bad']
        ahash = {}
    if (tn == 3):
        atype = self.h.getSurvName("c docetaxel_response")
        atypes = ['good', 'bad']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCardoso2016 = getCardoso2016


def getKrishnaSubramanian2012(self, tn=1):
    self.prepareData("COV42")
    atype = self.h.getSurvName("c Title")
    atypes = ['N', 'T']
    ahash = {'N1+2':0, 'T3+4':1, 'N3+4':0, 'N5+6':0, 'T5+6':1, 'T1+2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKrishnaSubramanian2012 = getKrishnaSubramanian2012


def getAuerbach2020(self, tn=1, dg=""):
    self.prepareData("COV43")
    atype = self.h.getSurvName("c cell type")
    atypes = ['HepaRG cells']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c chemical")
        atypes = ['D', 'S', 'A', 'C', 'B']
        ahash = {'DMSO':0,'aspirin':2,'sucrose':1,'caffeine':3, 'CDCA':4}
    if (tn == 3):
        drug = self.h.getSurvName("c chemical")
        dose = self.h.getSurvName("c dose")
        dose = [re.sub(" nM", "", str(k)) for k in dose]
        atype = [dose[i] if drug[i] == dg or drug[i] == 'DMSO' else None
                for i in range(len(drug))]
        atype = ['0' if drug[i] == 'DMSO' else atype[i]
                for i in range(len(drug))]
        atype = [float(atype[i]) if atype[i] is not None else None
                for i in range(len(drug))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAuerbach2020 = getAuerbach2020


def getWigger2019(self, tn=1):
    self.prepareData("COV44")
    atype = self.h.getSurvName("c time point")
    ahash = {'24h':24, '4h':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c ligand treatment")
    atypes = ['Control', 'CDCA', 'FXR-L', 'PPARa', 'LXR-L']
    ahash = {}
    if (tn == 2):
        atypes = ['Control', 'CDCA', 'FXR-L']
        atype = [atype[i] if tval[i] == 24 else None
                for i in range(len(atype))]
    if (tn == 3):
        drug = self.h.getSurvName("c ligand treatment")
        atype = self.h.getSurvName("c time point")
        atypes = ['4h', '24h']
        ahash = {}
        atype = [atype[i] if drug[i] == 'CDCA' else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWigger2019 = getWigger2019


def getIjssennagger2016(self, tn=1):
    self.prepareData("COV45")
    atype = self.h.getSurvName("c treatment")
    atypes = ['V', 'OCA']
    ahash = {'vehicle (0.1% DMSO)':0, '1 uM OCA (INT-747)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getIjssennagger2016 = getIjssennagger2016


def getIjssennagger2016Mm(self, tn=1):
    self.prepareData("COV45.2")
    atype = self.h.getSurvName("c genotype/variation")
    ahash = {'FXR KO':1, 'wild type':0}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['V', 'OCA']
    ahash = {'OCA (INT-747) (10 mg/kg/day, dissolved in 1% methyl cellulose':1,
            'vehicle (1% methyl cellulose)':0}
    if (tn == 2):
        atype = [atype[i] if gval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if gval[i] == 1 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getIjssennagger2016Mm = getIjssennagger2016Mm


def getShen2020(self, tn=1):
    self.prepareData("COV46")
    atype = self.h.getSurvName("c disase state")
    atypes = ['H', 'AIDS']
    ahash = {'AIDS':1, 'healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShen2020 = getShen2020


def getLiao2020(self, tn=1):
    self.prepareData("COV47")
    atype = self.h.getSurvName("c infection")
    atypes = ['C', 'Zika']
    ahash = {'control':0, 'Zika virus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiao2020 = getLiao2020


def getSuthar2019(self, tn=1):
    self.prepareData("COV48")
    atype = self.h.getSurvName("c Title")
    time = [ re.split("_", str(k))[2] if len(re.split("_", str(k))) > 2
            else None for k in atype]
    ahash = {'12H':12, '24H':24}
    tval = [ahash[i] if i in ahash else None for i in time]
    atype = [ re.split("_", str(k))[1] if len(re.split("_", str(k))) > 1
            else None for k in atype]
    atypes = ['Mock', 'WNV', 'RIG-I', 'MDA5', 'IFNb']
    ahash = {}
    if (tn == 2):
        atypes = ['Mock', 'WNV']
        atype = [None if atype[i] == 'WNV' and tval[i] == 12 else atype[i]
                for i in range(len(atype))]
    if (tn == 3):
        atypes = ['Mock', 'IFNb']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSuthar2019 = getSuthar2019


def getScott2019(self, tn=1):
    self.prepareData("COV49")
    atype = self.h.getSurvName("c culture condition")
    atypes = ['None', 'EGF', 'RAFT']
    ahash = {'':0, 'mouse EGF':1, '3D raft culture':2}
    if (tn == 2):
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.h.getSurvName("c infection")
        atype = [atype[i] if aval[i] == 1 else None
                for i in range(len(atype))]
        atypes = ['C', 'HPV']
        ahash = {'ECM only':0,
                '7d ECM only':0,
                '7dpi with HPV16':1,
                'transfected/immortalized with HPV16 genome':1,
                '27-33 d ECM only':0,
                '27-33 dpi with HPV16 WT':1,
                '36-55 dpi with HPV16 WT':1,
                '43-63 dpi with HPV16 WT':1,
                'uninfected control':0,
                'HPV16 WT':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getScott2019 = getScott2019


def getvanTol2020(self, tn=1):
    self.prepareData("COV50")
    atype = self.h.getSurvName("c infection")
    atypes = ['C', 'WNV']
    ahash = {'Mock with PBS':0, 'WNV at MOI = 5':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getvanTol2020 = getvanTol2020


def getDeutschman2019(self, tn=1):
    self.prepareData("COV51")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'Non':0, 'CAP-D3':1, 'CAP-H2':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['DMSO', 'VitK']
    ahash = {'DMSO':0, '50uM Menadione':1}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeutschman2019 = getDeutschman2019


def getDeferme2013(self, tn=1):
    self.prepareData("COV52")
    atype = self.h.getSurvName("c time")
    tval = [float(atype[i]) if i > 1 else None for i in range(len(atype))]
    atype = self.h.getSurvName("c compound, dose")
    atypes = ['C', 'Men', 'H2O2', 'TBH']
    ahash = {'100\xc2\xb5M Men':1,
            '50\xc2\xb5M H2O2':2,
            'Control':0,
            '200\xc2\xb5M TBH':3}
    if (tn == 2):
        atypes = ['C', 'VitK']
        ahash = {'100\xc2\xb5M Men':1,
                'Control':0}
        atype = [atype[i] if tval[i] == 8 else None
                for i in range(len(atype))]
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeferme2013 = getDeferme2013


def getBell2017(self, tn=1):
    self.prepareData("COV53")
    atype = self.h.getSurvName("c treatment time")
    ahash = {'48h':0, '7d':1, '14d':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'CPZ', 'Am', 'Af']
    ahash = {'DMSO':0, 'Chlorpromazine':1, 'Amiodarone':2, 'Aflatoxin B1':3}
    if (tn == 2):
        atypes = ['C', 'CPZ']
        ahash = {'DMSO':0, 'Chlorpromazine':1}
        atype = [atype[i] if tval[i] == 1 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBell2017 = getBell2017


def getDeAbrew2015(self, tn=1):
    self.prepareData("COV54")
    drug = self.h.getSurvName("c agent")
    drugs = hu.uniq(drug[2:])
    time = self.h.getSurvName("c time")
    conc = self.h.getSurvName("c concentration")
    comb = [str(time[i]) + " " + str(conc[i]) for i in range(len(time))]
    atype = self.h.getSurvName("c time")
    ahash = {'24 hours':24, '48 hours':48}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['24', '48']
    ahash = {'24 hours':0, '48 hours':1}
    if (tn == 2):
        atype = [comb[i] if drug[i] == 'Chlorpromazine HCl' or
                drug[i] == 'DMSO' or drug[i] == 'Water'
                else None for i in range(len(atype))]
        atype = ['0 ' + drug[i] + ' ' + comb[i] if 
                drug[i] == 'DMSO' or drug[i] == 'Water'
                else atype[i] for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        atypes = ['C', 'CPZ']
        ahash = {'0 DMSO 24 hours 1%':0, '48 hours 0.8 uM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeAbrew2015 = getDeAbrew2015


def getVandenHof2014(self, tn=1):
    self.prepareData("COV55")
    atype = self.h.getSurvName("c treatment time")
    ahash = {'24h':24}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['24']
    ahash = {'24h':0}
    if (tn == 2):
        atype = self.h.getSurvName("c compound, dose")
        atypes = ['C', 'CPZ']
        ahash = {'DMSO, 0.5 %':0, 'CLP, 20 \xc2\xb5M':1}
        atype = [atype[i] if tval[i] == 24 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName("c compound, dose")
        atypes = ['C', 'CP']
        ahash = {'DMSO, 0.5 %':0, 'CP, 2000 \xc2\xb5M':1}
        atype = [atype[i] if tval[i] == 24 else None
                for i in range(len(atype))]
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVandenHof2014 = getVandenHof2014


def getVitins2017(self, tn=1):
    self.prepareData("COV56")
    atype = self.h.getSurvName("c time point")
    ahash = {'25 days':25}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C1', 'C2', 'CPZ', 'EE2']
    ahash = {'EE2':3, 'CPZ':2, 'CPZ vehicle (PBS)':0,
            'EE2 vehicle (sunflower oil)':1}
    if (tn == 2):
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'CPZ']
        ahash = {'CPZ vehicle (PBS)':0, 'CPZ':1}
        atype = [atype[i] if tval[i] == 25 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVitins2017 = getVitins2017


def getBroin2016(self, tn=1):
    self.prepareData("COV57")
    atype = self.h.getSurvName("c treatment")
    atypes = ['control', 'ACEI', 'ARB']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBroin2016 = getBroin2016


def getMatsuuraHachiya2015mm(self, tn=1):
    self.prepareData("COV58")
    atype = self.h.getSurvName("c treatment")
    atype = [re.sub("l .*", "l", str(k)) for k in atype]
    atypes = ['control', 'ACEI']
    ahash = {'applied 30% ethanol':0, 'applied 1% enalapril':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatsuuraHachiya2015mm = getMatsuuraHachiya2015mm


def getAbdAlla2010mm(self, tn=1):
    self.prepareData("COV59")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(".*from ", "", str(k)) for k in atype]
    atypes = ['C', 'Ath', 'Ath+ACEI']
    ahash = {'non-transgenic C57BL/6J control mice':0,
            'captopril-treated APOE-deficient mice':2,
             'atherosclerotic APOE-deficient mice':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAbdAlla2010mm = getAbdAlla2010mm


def getEun2013rat(self, tn=1):
    self.prepareData("COV60")
    atype = self.h.getSurvName("c treatment")
    time = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'7':7, '2':2, '3':3, '10':10}
    tval = [ahash[i] if i in ahash else None for i in time]
    atype = [re.sub(".*after ", "", str(k)) for k in atype]
    atype = [re.sub(" treatment ", "", str(k)) for k in atype]
    atype = [re.sub(" dose", "", str(k)) for k in atype]
    group = [re.sub(".*from ", "", str(k)) for k in atype]
    drug = [re.sub("(.*)\\((.*)\\)", "\\1", str(k)) for k in group]
    dose = [re.sub(".*\\((.*)\\)", "\\1", str(k)) for k in group]
    ahash = {'high':3, 'middle':2, 'low':1, 'corn oil':0}
    dval = [ahash[i] if i in ahash else None for i in dose]
    atype = [time[i] + " " + dose[i] for i in range(len(atype))]
    atypes = sorted(hu.uniq(atype[2:]), reverse=1)
    atype = group
    atypes = ['C', 'T']
    ahash = {'RAN(high)':1, 'PZA(middle)':1, 'RAN(low)':1, 'CPZ(low)':1,
            'CPZ(high)':1, 'CBZ(high)':1, 'RAN(middle)':1, 'PZA(low)':1,
            'CBZ(middle)':1, 'PZA(high)':1, 'CPZ(middle)':1,
            'vehicle (corn oil)':0, 'CBZ(low)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEun2013rat = getEun2013rat


def getDelVecchio2014(self, tn=1):
    self.prepareData("COV61")
    atype = self.h.getSurvName("c agent")
    atypes = ['C', 'VitK', 'PERKi']
    ahash = {'PERKi':2, 'Menadione':1, 'DMSO':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelVecchio2014 = getDelVecchio2014


def getBriede2010I(self, tn=1):
    self.prepareData("COV62")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*FD .*", "B", str(k)) for k in atype]
    group = [re.sub(".* bi.*", "A", str(k)) for k in atype]
    atype = self.h.getSurvName("c Title")
    drug = [re.sub("T.*", "", str(k)) for k in atype]
    atype = self.h.getSurvName("c time")
    atype = [re.sub(" h", "", str(k)) for k in atype]
    time = [float(atype[i]) if i > 1 else None for i in range(len(atype))]
    atypes = sorted(hu.uniq(time[2:]))
    atype = time
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if drug[i] == 'Menadione'  and group[i] == 'B' 
                else None for i in range(len(atype))]
        atypes = ['C', 'VitK', 'T2']
        ahash = {0.08:0, 16.0:2, 0.5:0, 0.25:0, 4.0:1, 
                8.0:2, 2.0:1, 24.0:2, 1.0:1}
    if (tn == 3):
        atype = [atype[i] if drug[i] == 'H2O2' else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBriede2010I = getBriede2010I


def getBriede2010II(self, tn=1):
    self.prepareData("COV63")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*FD .*", "B", str(k)) for k in atype]
    group = [re.sub(".* bi.*", "A", str(k)) for k in atype]
    atype = self.h.getSurvName("c Title")
    drug = [re.sub("T.*", "", str(k)) for k in atype]
    atype = self.h.getSurvName("c time")
    atype = [re.sub(" h", "", str(k)) for k in atype]
    time = [float(atype[i]) if i > 1 else None for i in range(len(atype))]
    atypes = sorted(hu.uniq(time[2:]))
    atype = time
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if drug[i] == 'Menadione'  and group[i] == 'B' 
                else None for i in range(len(atype))]
        atypes = ['C', 'VitK', 'T2']
        ahash = {0.08:0, 16.0:2, 0.5:0, 0.25:0, 4.0:1, 
                8.0:2, 2.0:1, 24.0:2, 1.0:1}
    if (tn == 3):
        atype = [atype[i] if drug[i] == 'H2O2' else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBriede2010II = getBriede2010II


def getYau2008(self, tn=1):
    self.prepareData("COV64")
    atype = self.h.getSurvName("c src1")
    atypes = ['U', 'M', 'H', 'D', 'E', 'i']
    ahash = {'MCF7, E2 deprivation, 72hr':4,
            'MCF7, anti-ERa siRNA treatment, 72hr':5,
            'MCF7, 0.5mM H2O2, 8hr':2,
            'MCF7, 275mM diamide, 8hr':3,
            'MCF7,untreated':0,
            'MCF7, 10mM menadione, 8hr':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYau2008 = getYau2008


def getGusenleitner2014(self, tn=1):
    self.prepareData("COV65")
    tissue = self.h.getSurvName("c tissue")
    time = self.h.getSurvName("c time")
    conc = self.h.getSurvName("c dose")
    comb = [str(time[i]) + " " + str(conc[i]) for i in range(len(time))]
    drug = self.h.getSurvName("c compound")
    atype = self.h.getSurvName("c tissue")
    atypes = ['K', 'L', 'H', 'PH', 'TM']
    ahash = {'Kidney':0,
            'Liver':1,
            'Heart':2,
            'Primary rat hepatocytes':3,
            'Thigh muscle':4}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [drug[i] if aval[i] == 3 and (drug[i] == 'Chloroquine'
                or drug[i] == '') else None
                for i in range(len(drug))]
        atypes = ['Neg', 'Q']
        ahash = {'':0, 'Chloroquine': 1}
    if (tn == 3):
        atype = [drug[i] if aval[i] == 3 and (drug[i] == 'Chlorpromazine'
                or drug[i] == '') else None
                for i in range(len(drug))]
        atypes = ['Neg', 'CPZ']
        ahash = {'':0, 'Chlorpromazine': 1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGusenleitner2014 = getGusenleitner2014


def getReadhead2018I(self, tn=1):
    self.prepareData("COV66")
    atype = self.h.getSurvName("c treatment")
    atypes = ['DMSO', 'LOX', 'MPB']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getReadhead2018I = getReadhead2018I


def getReadhead2018II(self, tn=1):
    self.prepareData("COV67")
    drug = self.h.getSurvName("c perturbagen")
    atype = self.h.getSurvName("c perturbation type")
    ahash = {'vehicle':0, 'test':1, 'poscon':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c src1")
    atypes = ['Ccontrol', 'Cancer', 'SZ']
    ahash = {'hiPSC_control':0, 'cancerCells':1, 'hiPSC_SZ':2}
    if (tn == 2):
        atype = [drug[i] if drug[i] == 'hydroquinine'
                or tval[i] == 0 or tval[i] == 2 else None
                for i in range(len(drug))]
        for i in range(len(atype)):
            if (tval[i] == 0):
                atype[i] = 'Neg'
            if (tval[i] == 2):
                atype[i] = 'Pos'
        atypes = ['Neg', 'Q']
        ahash = {'hydroquinine': 1}
    if (tn == 3):
        atype = [drug[i] if drug[i] == 'Chlorpromazine'
                or tval[i] == 0 or tval[i] == 2 else None
                for i in range(len(drug))]
        for i in range(len(atype)):
            if (tval[i] == 0):
                atype[i] = 'Neg'
            if (tval[i] == 2):
                atype[i] = 'Pos'
        atypes = ['Neg', 'CPZ']
        ahash = {'Chlorpromazine': 1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getReadhead2018II = getReadhead2018II


def getDeGottardi2016(self, tn=1):
    self.prepareData("COV73")
    atype = self.h.getSurvName("c treatment")
    ahash = {'IgG1 Control Ab':0, 'Anti-IL15 Ab':1}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c timepoint")
    atypes = ['Pre', 'Post']
    ahash = {'Post Treatment':1, 'Pre Treatment':0, 'Prost Treatment':1}
    if (tn == 2):
        atype = [atype[i] if dval[i] == 1 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeGottardi2016 = getDeGottardi2016


def getDeGottardi2016Hs(self, tn=1):
    self.prepareData("COV73.2")
    atype = self.h.getSurvName("c treatment")
    ahash = {'IgG1 Control Ab':0, 'Anti-IL15 Ab':1}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c timepoint")
    atypes = ['Pre', 'Post']
    ahash = {'Post Treatment':1, 'Pre Treatment':0, 'Prost Treatment':1}
    if (tn == 2):
        atype = [atype[i] if dval[i] == 1 else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDeGottardi2016Hs = getDeGottardi2016Hs


def getLaszlo2016(self, tn=1):
    self.prepareData("COV74")
    atype = self.h.getSurvName('c stimulated with')
    atypes = ['WT IL2', 'F42K IL2', 'IL15']
    ahash = {'0.4nM F42K IL-2 for 24 hr':1, '0.4nM IL-15 for 24 hr':2,
            '0.4nM WT IL-2 for 24 hr':0}
    if (tn == 2):
        atypes = ['IL2', 'IL15']
        ahash = {'0.4nM IL-15 for 24 hr':1,
                '0.4nM WT IL-2 for 24 hr':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLaszlo2016 = getLaszlo2016


def getShao2013(self, tn=1):
    self.prepareData("COV75")
    atype = self.h.getSurvName("c passage")
    atypes = ['P19', 'P18', 'P17', 'P16']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c agent")
        atypes = ['C', 'S1PR']
        ahash = {'Untreated cells':0, 'Fingolimod 4uM':1}
    if (tn == 3):
        atype = self.h.getSurvName("c agent")
        atypes = ['U', 'C']
        ahash = {'Untreated cells':0, 'Cyclophosphamide 5mM':1}
    if (tn == 4):
        atype = self.h.getSurvName("c agent")
        atypes = ['U', 'C']
        ahash = {'Untreated cells':0, 'Cyclophosphamide S9 treated 3mM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShao2013 = getShao2013


def getMelo2020(self, tn=1):
    self.prepareData("COV76.3")
    atype = self.h.getSurvName("c src1")
    atypes = ['C1', 'C2', 'RSV', 'IAV', 'CoV1', 'CoV2']
    ahash = {'Mock treated NHBE cells':0,
            'SARS-CoV-2 infected NHBE cells':4,
            'Mock treated A549 cells':1,
            'SARS-CoV-2 infected A549 cells':5,
            'RSV infected A549 cells':2,
            'IAV infected A549 cells':3}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['C', 'CoV']
        ahash = {'SARS-CoV-2 infected A549 cells':1,
                'Mock treated A549 cells':0}
    if (tn == 3):
        atypes = ['C', 'CoV']
        ahash = {'Mock treated NHBE cells':0,
                'SARS-CoV-2 infected NHBE cells':1}
    if (tn == 4):
        atype = aval
        atypes = ['C', 'RSV', 'IAV', 'CoV']
        ahash = {1:0, 2:1, 3:2, 5:3}
    if (tn == 5):
        atype = self.h.getSurvName("c subject status")
        atypes = ['C', 'CoV']
        ahash = {'No treatment - healthy 72 years old, male':0,
                'No treatment - healthy 77 years old, male':0,
                'No treatment; >60 years old male COVID-19 deceased patient':1}
    if (tn == 6):
        atypes = ['C', 'CoV']
        ahash = {'SARS-CoV-2 infected A549 cells':1,
                'Mock treated A549 cells':0,
                'Mock treated NHBE cells':0,
                'SARS-CoV-2 infected NHBE cells':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMelo2020 = getMelo2020


def getMelo2020II(self, tn=1, ta=0, tb=1):
    self.prepareData("COV76.2")
    h = self.h
    atype = h.getSurvName("c tissue/cell type")
    ahash = {'Nasal Wash':0, 'Trachea':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName("c time after treatment")
    ahash = {'1 day':1, '3 days':3, '7 days':7, '14 days':14}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c src1')
    atypes = ['C', 'CoV', 'IAV']
    ahash = {'Mock treated 4 month old Ferret':0,
            'SARS-CoV-2 infected 4 month old Ferret':1,
            'IAV infected 4 month old Ferret':2}
    if (tn == 2):
        atype = [atype[i] if gval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if gval[i] == 1 else None
                for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if gval[i] == 0 and 
                (tval[i] == 3 or tval[i] == 7) else None
                for i in range(len(atype))]
    if (tn == 5):
        atypes = ['C', 'CoV']
        ahash = {'Mock treated 4 month old Ferret':0,
                'SARS-CoV-2 infected 4 month old Ferret':1}
        atype = [atype[i] if gval[i] == ta and tval[i] == tb else None
                for i in range(len(atype))]
    if (tn == 6):
        atypes = ['C', 'CoV']
        ahash = {'Mock treated 4 month old Ferret':0,
                'SARS-CoV-2 infected 4 month old Ferret':1}
        atype = [atype[i] if gval[i] == 0 and 
                (tval[i] == 3 or tval[i] == 7) else None
                for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMelo2020II = getMelo2020II


def getPalmieri2017(self, tn=1):
    self.prepareData("COV78")
    atype = self.h.getSurvName("c treatment")
    atypes = ['U', 'T']
    ahash = {'100 mM of trehalose':1, 'Untreated':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPalmieri2017 = getPalmieri2017


def getOzgyin2018(self, tn=1):
    self.prepareData("COV79")
    atype = self.h.getSurvName("c treatment")
    atypes = ['U', 'T']
    ahash = {'non-lyophilized (control)':0, 'Lyophilized':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOzgyin2018 = getOzgyin2018


def getMeugnier2011(self, tn=1):
    self.prepareData("COV80")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(".*after (.*) treatment", "\\1", str(k)) for k in atype]
    atypes = ['E', 'A']
    ahash = {'Etanercept':0, 'Adalimumab':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeugnier2011 = getMeugnier2011


def getZaba2007(self, tn=1):
    self.prepareData("COV81")
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'T']
    ahash = {'Day 5 DC control':0, 'Day 5 DC etanercept':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZaba2007 = getZaba2007


def getSeelbinder2020(self, tn=1):
    self.prepareData("COV83")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'AF', 'E', 'EAF']
    ahash = {'none':0,
            'A. fumigatus ATCC 46645':1,
            'Etanercept 2 ug/mL':2,
            'Etanercept 2 ug/mL; A. fumigatus ATCC 46645':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSeelbinder2020 = getSeelbinder2020


def getBolenhcv(self, tn=1):
    self.prepareData("COV84")
    atype = self.h.getSurvName("c infection (ch1)")
    atypes = ['healthy', 'infected']
    ahash = {'Healthy':0, 'HCV':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBolenhcv = getBolenhcv

    
def getMoalhev(self, tn=1):
    self.prepareData("COV86")
    atype = self.h.getSurvName("c patient type (ch1)")
    atypes = ['healthy', 'infected']
    ahash = {'Control':0, 'Infected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMoalhev = getMoalhev


def getBrandes2013(self, tn=1):
    self.prepareData("COV87")
    atype = self.h.getSurvName("c cell type")
    ahash = {'':0, 'alveolar macrophage':1, 'lymphocyte (BC, TC, NK)':2,
            'Ly6Chi mononuclear myeloid cell':3, 'neutrophil':4,
            'pulmonary CD45neg':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c perturbation")
    atypes = ['C', 'A', 'T', 'P']
    ahash = {'H1N1 influenza A PR8 (100LD50)':3,
            'H1N1 influenza A PR8 (0.2LD50)':3,
            'H1N1 influenza A PR8 (0.6LD50)':3,
            'H1N1 influenza A PR8 (10LD50)':3,
            'H1N1 influenza A 0.6LD50 PR8':3,
            'ALUM (162g)':1,
            'H1N1 influenza A TX91 (10^6PFU)':2,
            'H1N1 influenza A TX91 (10^6 PFU)':2,
            'sham':0}
    if (tn >= 2 and tn < 6):
        atypes = ['NL', 'SL', 'L']
        ahash = {'H1N1 influenza A PR8 (100LD50)':2,
                'H1N1 influenza A PR8 (10LD50)':2,
                'H1N1 influenza A PR8 (0.2LD50)':1,
                'H1N1 influenza A PR8 (0.6LD50)':1,
                'H1N1 influenza A 0.6LD50 PR8':1,
                'H1N1 influenza A TX91 (10^6PFU)':0,
                'H1N1 influenza A TX91 (10^6 PFU)':0}
        atype = [atype[i] if tval[i] == (tn - 1) else None
                for i in range(len(atype))]
    if (tn == 7):
        atype = self.h.getSurvName("c cell type")
        ahash = {'alveolar macrophage':1,
                'Ly6Chi mononuclear myeloid cell':0}
        atypes = ['Mo', 'AM']
    if (tn == 8):
        atype = self.h.getSurvName("c cell type")
        ahash = {'alveolar macrophage':1,
                'neutrophil':0}
        atypes = ['Neu', 'AM']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBrandes2013 = getBrandes2013


def getJung2016(self, tn=1):
    self.prepareData("COV88")
    atype = self.h.getSurvName("c bmi")
    atypes = ['N', 'Obese', 'Overweight']
    ahash = {'18.5~23 kg/m2':0, '27.5~30 kg/m2':2, '25~27.5 kg/m2':1}
    if (tn == 2):
        atypes = ['N', 'Obese']
        ahash = {'18.5~23 kg/m2':0, '25~27.5 kg/m2':1}
    if (tn == 3):
        atypes = ['N', 'Overweight']
        ahash = {'18.5~23 kg/m2':0, '27.5~30 kg/m2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJung2016 = getJung2016


def getEckle2013(self, tn=1):
    self.prepareData("COV89")
    atype = self.h.getSurvName("c treatment")
    atypes = ['NS', 'S']
    ahash = {'non-strained':0, 'strained':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEckle2013 = getEckle2013


def getOmura2018(self, tn=1):
    self.prepareData("COV90")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("([1a]) .*", "\\1", str(k)) for k in atype]
    atypes = ['CN', 'IN', 'CH', 'IH', 'UV']
    ahash = {'Mock-infected Neuro-2a':0,
            'Mock-infected HL-1':2,
            'HL-1':4,
            'TMEV-infected HL-1':3,
            'TMEV-infected Neuro-2a':1}
    if (tn == 2):
        atypes = ['C', 'I']
        ahash = {'Mock-infected HL-1':0, 'TMEV-infected HL-1':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOmura2018 = getOmura2018


def getOmura2014(self, tn=1):
    self.prepareData("COV91")
    atype = self.h.getSurvName('c time point')
    ahash = {'4 days post infection':4, '7 days post infection':7,
            '60 days post infection':60}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infected with")
    atypes = ['C', 'I']
    ahash = {"Theiler's murine encephalomyelitis virus (TMEV)":1,
            'none (naïve control)':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 4 else None
            for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOmura2014 = getOmura2014


def getCoronado2012(self, tn=1):
    self.prepareData("COV92")
    atype = self.h.getSurvName('c time')
    ahash = {'10 dpi':0, '90 dpi':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c gender')
    ahash = {'Female':0, 'Male':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'I']
    ahash = {'PBS':0, 'CVB3':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 and gval[i] == 1 else None
            for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'FI', 'MI']
        atype = [atype[i] if tval[i] == 0 else None
            for i in range(len(atype))]
        atype = ['MI' if aval[i] == 1 and tval[i] == 0 and
                gval[i] == 1 else atype[i] for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCoronado2012 = getCoronado2012


def getClementeCasares2017(self, tn=1):
    self.prepareData("COV93")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(" rep.*", "", str(k)) for k in atype]
    atypes = ['DC103', 'DC11b', 'MF']
    ahash = {'Cardiac CD103+ DC':0,
            'Cardiac CD11b+ DC':1,
            'Cardiac MHC-IIhi MF':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getClementeCasares2017 = getClementeCasares2017


def getSchinke2004(self, tn=1):
    self.prepareData("COV95")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("PGA-(.*)-.*", "\\1", str(k)) for k in atype]
    atypes = ['CCMP', 'N', 'AS']
    ahash = {}
    if (tn == 2):
        atypes = ['N', 'CCMP']
    if (tn == 3):
        atypes = ['N', 'AS']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSchinke2004 = getSchinke2004


def getSchinke2004II(self, tn=1):
    self.prepareData("COV96")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("PGA[-_](.*)[-_][0-9]+", "\\1", str(k)) for k in atype]
    atype = [re.sub("_", "-", str(k)) for k in atype]
    atypes = ['Hs-V', 'Hs-S', 'PA-D', 'Hs-D', 'PA-N', 'Hs-F',
            'Hs-H', 'PA-S', 'Hs-P']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSchinke2004II = getSchinke2004II


def getClarelli2017(self, tn=1):
    self.prepareData("COV97")
    atype = self.h.getSurvName("c agent")
    atypes = ['C', 'IFNb']
    ahash = {'IFN-beta':1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getClarelli2017 = getClarelli2017


def getHu2012(self, tn=1):
    self.prepareData("COV98")
    atype = self.h.getSurvName("c src1")
    atypes = ['C', 'LPS', 'TNFa', 'S', 'S+L', 'S+Ta']
    ahash = {'BEAS-2B - control':0,
            'BEAS-2B - mechanical stretch plus TNF-\xce\xb1':5,
            'BEAS-2B - mechanical stretch plus LPS':4,
            'BEAS-2B - LPS':1,
            'BEAS-2B - mechanical stretch':3,
            'BEAS-2B - TNF-\xce\xb1':2}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'BEAS-2B - control':0,
                'BEAS-2B - mechanical stretch':1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHu2012 = getHu2012


def getdosSantos2004(self, tn=1):
    self.prepareData("COV99")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*from ", "", str(k)) for k in atype]
    atype = [re.sub("RNA ", "", str(k)) for k in atype]
    atype = [re.sub("hr .*", "", str(k)) for k in atype]
    atype = [re.sub(".ells[ \\-]*", "", str(k)) for k in atype]
    atype = [re.sub(".tatic ", "", str(k)) for k in atype]
    atype = [re.sub(" [JO][ac][nt].*", "", str(k)) for k in atype]
    atype = [re.sub(" [14]", "", str(k)) for k in atype]
    atype = [re.sub("A549 ", "", str(k)) for k in atype]
    atypes = ['C', 'LPS', 'TNFa', 'S', 'S+Ta']
    ahash = {'Control':0, 'LPS':1, 'LPSh':1, 'TNF':2,
            'Stretch':3, 'stretch':3, 'TNF+Stretch':4}
    if (tn == 2):
        atype = self.h.getSurvName("c Title")
        atypes = ['C', 'S']
        ahash = {'A549 TNF+Stretch 1hr Jan 31 2003':1,
                'mRNA from A549 Static Control 1hr Jan 31 2003':0,
                'A549 TNF+Stretch 4hr Jan 31 2003':1,
                'RNA A549 cells Static Control 4hr Replicate 1':0,
                'A549 stretch 4hr Jan 31 2003':1,
                'A549 Stretch 1hr Jan 31 2003':1,
                'RNA from A549 Cells - Static Control 1hr Replicate 1':0,
                'A549 RNA Static Control 4hr Jan 31 2003':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getdosSantos2004 = getdosSantos2004


def getSwanson2012(self, tn=1):
    self.prepareData("COV100")
    atype = self.h.getSurvName("c src1")
    atypes = ['noVAP', 'VAP']
    ahash = {'patients without VAP':0, 'patients with VAP':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSwanson2012 = getSwanson2012


def getGharib2014(self, tn=1):
    self.prepareData("COV101")
    atype = self.h.getSurvName("c tissue type")
    atypes = ['Kidney', 'Lung', 'Liver']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub("_.$", "", str(k)) for k in atype]
        atypes = ['C', 'MV']
        ahash = {'Kidney_Control':0,
                'Lung_MV+SA':1,
                'Kidney_MV+SA':1,
                'Lung_Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGharib2014 = getGharib2014


def getStaudt2018(self, tn=1):
    self.prepareData("COV103")
    atype = self.h.getSurvName("c tissue")
    atypes = ['AM', 'SAE']
    ahash = {'alveolar macrophage':0, 'small airway epithelium brushing':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c smoking status")
        atypes = ['NS', 'ENN', 'EN']
        ahash = {'nonsmoker':0, 'Ecig_no-nicotine':1, 'Ecig+nicotine':2}
        atype = [atype[i] if aval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName("c smoking status")
        atypes = ['NS', 'ENN', 'EN']
        ahash = {'nonsmoker':0, 'Ecig_no-nicotine':1, 'Ecig+nicotine':2}
        atype = [atype[i] if aval[i] == 1 else None
                for i in range(len(atype))]
    if (tn == 4):
        atypes = ['SAE', 'AM']
        ahash = {'alveolar macrophage':1,
                'small airway epithelium brushing':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getStaudt2018 = getStaudt2018


def getQuast2019(self, tn=1):
    self.prepareData("COV104")
    sex = self.h.getSurvName("c Sex")
    atype = self.h.getSurvName("c disease state")
    atypes = ['H', 'COPD']
    ahash = {'healthy':0, 'COPD':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c smoker")
        atypes = ['N', 'Y']
        ahash = {}
    if (tn == 3):
        atype = sex
        atypes = ['F', 'M']
        ahash = {}
    if (tn == 4):
        atype = self.h.getSurvName("c treatment")
        atypes = ['0', '4']
        ahash = {}
        atype = [atype[i] if aval[i] == 0 else None
                for i in range(len(atype))]
    if (tn == 5):
        atype = self.h.getSurvName("c treatment")
        atypes = ['Air', 'CS']
        ahash = {'0':0, '4':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQuast2019 = getQuast2019


def getYu2011(self, tn=1):
    self.prepareData("COV105")
    atype = self.h.getSurvName("c infection")
    atypes = ['U', 'EV71']
    ahash = {'uninfected':0, 'EV71':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYu2011 = getYu2011


def getMolinaNavarro2013(self, tn=1):
    self.prepareData("COV106")
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'DC', 'IC']
    ahash = {'Ischemic cardiomyopathy':2,
            'Dilated cardiomyopathy':1,
            'Normal heart':0}
    if (tn == 2):
        atypes = ['N', 'DC']
        ahash = {'Dilated cardiomyopathy':1,
                'Normal heart':0}
    if (tn == 3):
        atypes = ['N', 'IC']
        ahash = {'Ischemic cardiomyopathy':1,
                'Normal heart':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMolinaNavarro2013 = getMolinaNavarro2013


def getLi2020(self, tn=1):
    self.prepareData("COV107")
    atype = self.h.getSurvName("c src1")
    atypes = ['LV', 'ARV']
    ahash = {"ARVC patients' heart LV":0, 'heart RV tissue':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2020 = getLi2020


def getRen2020(self, tn=1):
    self.prepareData("COV108")
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'HF', 'HCM']
    ahash = {'heart failure':1, 'hypertrophic cardiomyopathy':2,
            'normal heart tissue':0}
    if (tn == 2):
        atypes = ['N', 'HF']
        ahash = {'heart failure':1, 'normal heart tissue':0}
    if (tn == 3):
        atypes = ['N', 'HCM']
        ahash = {'hypertrophic cardiomyopathy':1, 'normal heart tissue':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRen2020 = getRen2020


def getRen2020Mm(self, tn=1):
    self.prepareData("COV108.2")
    atype = self.h.getSurvName("c src1")
    atypes = ['N', 'C', 'CM']
    ahash = {'TAC sham CM':1, 'TAC2W CM':2, 'TAC5W CM':2,
            'TAC8W CM':2, 'TAC11W CM':2, 'normal heart CM':0}
    if (tn == 2):
        atypes = ['C', 'CM']
        ahash = {'TAC sham CM':0, 'TAC2W CM':1, 'TAC5W CM':1,
                'TAC8W CM':1, 'TAC11W CM':1, 'normal heart CM':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRen2020Mm = getRen2020Mm


def getMorley2019(self, tn=1):
    self.prepareData("COV109")
    atype = self.h.getSurvName("c etiology")
    atypes = ['N', 'DCM', 'HCM', 'PPCM']
    ahash = {'Dilated cardiomyopathy (DCM)':1,
            'Non-Failing Donor':0,
            'Hypertrophic cardiomyopathy (HCM)':2,
            'Peripartum cardiomyopathy (PPCM)':3}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['H', 'CM']
        ahash = {'Dilated cardiomyopathy (DCM)':1,
                'Non-Failing Donor':0,
                'Hypertrophic cardiomyopathy (HCM)':1,
                'Peripartum cardiomyopathy (PPCM)':1}
    if (tn == 3):
        atype = self.h.getSurvName("c race")
        atypes = ['W', 'B']
        ahash = {'Caucasian':0, 'African American':1}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMorley2019 = getMorley2019


def getLiu2015(self, tn=1):
    self.prepareData("COV110")
    age = self.h.getSurvName("c age")
    age = [int(age[i]) if i > 1 else None for i in range(len(age))]
    sex = self.h.getSurvName("c gender")
    atype = self.h.getSurvName("c disease status")
    atypes = ['NF', 'ICM', 'CMP']
    ahash = {'ischemic':1, 'non-failing':0, 'idiopathic dilated CMP':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if sex[i] == 'male' and age[i] < 40 
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if sex[i] == 'female' and age[i] < 60
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = age
        ahash = {}
        atypes = ['Y', 'O']
        for i in range(len(atype)):
            if age[i] < 30 and age[i] < 10:
                ahash[atype[i]] = 0
            if age[i] < 30 and age[i] > 10:
                ahash[atype[i]] = 1
        atype = [atype[i] if aval[i] == 0 
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = self.h.getSurvName("c heart failure")
        atypes = ['no', 'yes']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiu2015 = getLiu2015


def getAkat2014(self, tn=1):
    self.prepareData("COV111")
    h = self.h
    atype = h.getSurvName("c disease")
    atypes = ['NF', 'ICM', 'DCM']
    ahash = {'None':0, 'Ischemic Cardiomyopathy':1, 'Dilated Cardiomyopathy':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = h.getSurvName('c treatment')
        atypes = ['N', 'E', 'I']
        ahash = {'None':0, 'LVAD Explantation':1, 'LVAD Implantation':2}
    if (tn == 3):
        atype = h.getSurvName('c treatment')
        atypes = ['N', 'E', 'I']
        ahash = {'None':0, 'LVAD Explantation':1, 'LVAD Implantation':2}
        atype = [atype[i] if aval[i] == 2 or atype[i] == 'None'
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = h.getSurvName('c treatment')
        atypes = ['N', 'E', 'I']
        ahash = {'None':0, 'LVAD Explantation':1, 'LVAD Implantation':2}
        atype = [atype[i] if aval[i] == 1 or atype[i] == 'None'
                else None for i in range(len(atype))]
    if (tn == 5):
        atypes = ['NF', 'CM']
        ahash = {'None':0, 'Ischemic Cardiomyopathy':1, 'Dilated Cardiomyopathy':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAkat2014 = getAkat2014


def getCasey2016(self, tn=1):
    self.prepareData("COV112")
    age = self.h.getSurvName("c age")
    age = [int(age[i]) if i > 1 else None for i in range(len(age))]
    sex = self.h.getSurvName("c Sex")
    atype = sex
    atypes = ['F', 'M']
    ahash = {}
    if (tn == 2):
        atype = age
        ahash = {}
        atypes = ['Y', 'O']
        for i in range(len(atype)):
            if age[i] is None:
                continue
            if age[i] > 0 and age[i] < 55:
                ahash[atype[i]] = 0
            if age[i] > 0 and age[i] > 75:
                ahash[atype[i]] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCasey2016 = getCasey2016


def getHannenhalli2006(self, tn=1):
    self.prepareData("COV113")
    atype = self.h.getSurvName("c src1")
    atypes = ['NF', 'LV']
    ahash = {'explanted heart tissue at time of cardiac transplantation':1,
             'unused donor heart with normal LV function':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHannenhalli2006 = getHannenhalli2006


def getJacobson2016(self, tn=1, ar=None):
    self.prepareData("COV114")
    time = self.h.getSurvName("c timepoint")
    ahash = {'Week12':12, '12wk':12, '0wk':0, 'Week0':0}
    time = [ahash[i] if i in ahash else None for i in time]
    atype = self.h.getSurvName("c src1")
    atypes = ['U', 'T']
    ahash = {'Blood_untreated_chronic_HIV_placebo':0,
            'Blood_antiretroviral_therapy_placebo':1,
            'Blood_antiretroviral_therapy_Chloroquine':1,
            'Blood_untreated_chronic_HIV_Chloroq':0,
            'Blood_untreated_chronic_HIV_Chloroq_placebo':0}
    if (tn == 2):
        ahash = {'Blood_untreated_chronic_HIV_placebo': 'P',
                'Blood_untreated_chronic_HIV_Chloroq': 'Q',
                'Blood_untreated_chronic_HIV_Chloroq_placebo': 'P'}
        treat = [ahash[i] if i in ahash else None for i in atype]
        atype = [ahash[atype[i]] + ' ' + str(time[i]) if atype[i] in ahash
                else None for i in range(len(atype))]
        if (ar is not None):
            atype = [atype[i] if time[i] == ar
                    else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    if (tn == 3):
        ahash = {'Blood_antiretroviral_therapy_placebo': 'P',
                'Blood_antiretroviral_therapy_Chloroquine': 'Q'}
        treat = [ahash[i] if i in ahash else None for i in atype]
        atype = [ahash[atype[i]] + ' ' + str(time[i]) if atype[i] in ahash
                else None for i in range(len(atype))]
        if (ar is not None):
            atype = [atype[i] if time[i] == ar
                    else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    if (tn == 4):
        ahash = {'Blood_untreated_chronic_HIV_placebo': 'P',
                'Blood_untreated_chronic_HIV_Chloroq': 'Q',
                'Blood_untreated_chronic_HIV_Chloroq_placebo': 'P'}
        treat = [ahash[i] if i in ahash else None for i in atype]
        atype = [ahash[atype[i]] + ' ' + str(time[i]) if atype[i] in ahash
                else None for i in range(len(atype))]
        if (ar is not None):
            atype = [atype[i] if treat[i] == ar
                    else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    if (tn == 5):
        ahash = {'Blood_antiretroviral_therapy_placebo': 'P',
                'Blood_antiretroviral_therapy_Chloroquine': 'Q'}
        treat = [ahash[i] if i in ahash else None for i in atype]
        atype = [ahash[atype[i]] + ' ' + str(time[i]) if atype[i] in ahash
                else None for i in range(len(atype))]
        if (ar is not None):
            atype = [atype[i] if treat[i] == ar
                    else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJacobson2016 = getJacobson2016


def getLhakhang2014(self, tn=1):
    self.prepareData("COV115")
    atype = self.h.getSurvName("c treatment")
    atypes = ['U', 'T']
    ahash = {'untreated':0, 'hY3 treated':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLhakhang2014 = getLhakhang2014


def getTakeda2018(self, tn=1):
    self.prepareData("COV116")
    atype = self.h.getSurvName("c treatment")
    atypes = ['V', 'Q', 'M', 'M+O']
    ahash = {'mefloquine with trans-l-diaminocyclohexane oxalatoplatinum':3,
            'vehicle':0,
            'chloroquine':1,
            'mefloquine':2}
    if (tn == 2):
        atypes = ['V', 'M']
        ahash = {'vehicle':0,
                'mefloquine':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTakeda2018 = getTakeda2018


def getBargiela2019(self, tn=1):
    self.prepareData("COV117")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".$", "", str(k)) for k in atype]
    atypes = ['MC', 'MD', 'Q 0.1M', 'Q 10M']
    ahash = {'Myoblast_Control':0,
            'Myoblast_Treated_01_S':2,
            'Myoblast_Treated_10_S':3,
            'Myoblast_Disease':1}
    if (tn == 2):
        atypes = ['MC', 'Q']
        ahash = {'Myoblast_Control':0,
                'Myoblast_Treated_01_S':1,
                'Myoblast_Treated_10_S':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBargiela2019 = getBargiela2019


def getGoldberg2018(self, tn=1):
    self.prepareData("COV118")
    meds = self.h.survhdrs[15:24]
    atype = self.h.getSurvName(meds[0])
    for k in meds[1:]:
        atype1 = self.h.getSurvName(k)
        atype = [atype[i] + atype1[i] for i in range(len(atype))]
    drugs = atype
    atype = self.h.getSurvName('c gender')
    atypes = ['F', 'M']
    ahash = {}
    if (tn == 2):
        atype = drugs
        atypes = ['C', 'M', 'I']
        ahash = {'--Celebrex--------------':0,
                'Methotrexate----------------':1,
                '------Infliximab----------':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGoldberg2018 = getGoldberg2018


def getGrinman2019(self, tn=1):
    self.prepareData("COV119")
    atype = self.h.getSurvName('c stage of mammary development')
    atypes = ['A', 'D']
    ahash = {'Secretory Activation':0,
            'Secretory Differentiation':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGrinman2019 = getGrinman2019


def getThomas2018(self, tn=1):
    self.prepareData("COV120")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'LPS', 'T0+LPS']
    ahash = {'DMSO':0, 'DMSO+LPS':1, 'T0+LPS':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getThomas2018 = getThomas2018


def getSallam2018(self, tn=1):
    self.prepareData("COV121")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'T']
    ahash = {'Macrophage (DMSO)':0, 'Macrophage (GW3965)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSallam2018 = getSallam2018


def getDun2013(self, tn=1):
    self.prepareData("COV122")
    atype = self.h.getSurvName('c src1')
    atype = [re.sub(".*, (.*)h .*", "\\1", str(k)) for k in atype]
    atypes = ['C', 'T']
    ahash = {'72':1, '0':0, '24':1, '12':1, '48':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDun2013 = getDun2013


def getWood2016(self, tn=1):
    self.prepareData("COV123")
    atype = self.h.getSurvName('c disease state')
    atypes = ['C', 'SSc', 'Mor']
    ahash = {'Systemic Sclerosis':1, 'SSc patient':1,
            '':0, 'normal control':0, 'normal':0, 'morphea patient':2}
    if (tn == 2):
        atype = self.h.getSurvName('c treatment')
        atypes = ['U', 'Q']
        ahash = {'':0, 'NA':0, 'Plaquenil':1}
    if (tn == 3):
        atype = self.h.getSurvName('c treatment')
        atypes = ['U', 'MMF']
        ahash = {'':0, 'NA':0, 'mycophenolate mofetil':1, 'MMF':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWood2016 = getWood2016


def getFribourg2020(self, tn=1):
    self.prepareData("COV124")
    atype = self.h.getSurvName('c tacrolimus withdrawal')
    atypes = ['NW', 'WS', 'WR']
    ahash = {'Withdrawal Rejection (WR)':2,
            'Withdrawal Stable (WS)':1,
            'No withdrawal (NW)':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFribourg2020 = getFribourg2020


def getMunz2020(self, tn=1):
    self.prepareData("COV125")
    atype = self.h.getSurvName('c treatment group')
    atypes = ['M', 'T']
    ahash = {'Mock (FK506\xe2\x80\x93)':0, 'FK506 (FK506+)':1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMunz2020 = getMunz2020


def getDorr2015(self, tn=1):
    self.prepareData("COV126")
    atype = self.h.getSurvName('c time of blood draw')
    atypes = ['Pre', '1w', '3m', '6m']
    ahash = {'3 months post-transplant':2,
            'pre-transplant':0,
            '1 week post-transplant':1,
            '6 months post-transplant':3}
    if (tn == 2):
        atype = self.h.getSurvName('c race')
        atypes = ['W', 'AA', 'AI']
        ahash = {'American Indian orAlaskaNative':2,
                'Caucasian / White':0,
                'Black or African American':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDorr2015 = getDorr2015


def getCamus2012(self, tn=1):
    self.prepareData("COV127")
    h = self.h
    atype = h.getSurvName('c Title')
    atype = [re.sub("\xc2\xb5", "u", str(k)) for k in atype]
    atype = [re.sub("\xb5", "u", str(k)) for k in atype]
    time = [re.sub(".*([0-9]+)h.*", "\\1", str(k)) for k in atype]
    treat = [re.sub(".*(CQ).*", "CQd", str(k)) for k in atype]
    treat = [re.sub(".*(uninfe).*", "U", str(k)) for k in treat]
    treat = [re.sub(".*(uninfe).*", "U", str(k)) for k in treat]
    treat = [re.sub("^H.*", "U", str(k)) for k in treat]
    infection = [re.sub(".*(JFH1).*", "I", str(k)) for k in atype]
    infection = [re.sub("^H.*", "U", str(k)) for k in infection]
    atype = [treat[i] + " " + infection[i] for i in range(len(atype))]
    atypes = ['U U', 'CQd U', 'CQd I', 'U I']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCamus2012 = getCamus2012


def getAlao2018(self, tn=1):
    self.prepareData("COV128")
    h = self.h
    atype = h.getSurvName('c treatment')
    atypes = ['U', 'T']
    ahash = {'Base':0, 'Wk2':1, 'Wk4':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAlao2018 = getAlao2018


def getMeissner2014I(self, tn=1):
    self.prepareData("COV129")
    h = self.h
    atype = h.getSurvName('c time point')
    atypes = ['U', 'T']
    ahash = {'day 24':1, 'day 10':1,
            'end of treatment (post)':1,
            'pre treatment':0, 'day 0':0,
            'end of treatment (post-C)':1,
            'pre-C treatment':0}
    if (tn == 2):
        ahash = {'day 24':1, 'day 10':1,
                'end of treatment (post)':1,
                'pre treatment':0, 'day 0':0,
                'end of treatment (post-C)':0,
                'pre-C treatment':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeissner2014I = getMeissner2014I


def getMeissner2014II(self, tn=1):
    self.prepareData("COV130")
    h = self.h
    atype = h.getSurvName('c tissue')
    atypes = ['S', 'R']
    ahash = {'unpaired liver biopsy with sustained virologic response':0,
            'unpaired liver biopsy with relapse':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeissner2014II = getMeissner2014II


def getCostarelli2017(self, tn=1):
    self.prepareData("COV131")
    atype = self.h.getSurvName('c treatment')
    atypes = ['U', 'OMP', 'LSP']
    ahash = {'lansoprazole':2, 'untreated':0, 'omeprazole':1}
    if (tn == 2):
        atype = self.h.getSurvName('c src1')
        atype = [re.sub(".*cells, ", "", str(k)) for k in atype]
        ahash = {'young, control':0,
                'young, omeprazole':1, 'young,omeprazole':1}
        atypes = ['U', 'OMP']
    if (tn == 3):
        atype = self.h.getSurvName('c src1')
        atype = [re.sub(".*cells, ", "", str(k)) for k in atype]
        ahash = {'old, control':0, 'old, omeprazole':1}
        atypes = ['U', 'OMP']
    if (tn == 4):
        atype = self.h.getSurvName('c src1')
        atype = [re.sub(".*cells, ", "", str(k)) for k in atype]
        ahash = {'young, control':0,
                'young, omeprazole':0, 'young,omeprazole':0,
                'old, control':1, 'old, omeprazole':1}
        atypes = ['Y', 'O']
    if (tn == 5):
        atype = self.h.getSurvName('c src1')
        atype = [re.sub(".*cells, ", "", str(k)) for k in atype]
        ahash = {'young, control':0,
                'young, omeprazole':1, 'young,omeprazole':1,
                'old, control':2, 'old, omeprazole':3}
        atypes = ['YU', 'YOMP', 'OU', 'OOMP']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCostarelli2017 = getCostarelli2017


def getAbdelAziz2016(self, tn=1):
    self.prepareData("COV132")
    h = self.h
    atype = h.getSurvName('c Title')
    tissue = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'whole':0, 'esophageal':1}
    tval = [ahash[i] if i in ahash else None for i in tissue]
    treat = [re.sub("whole blood ", "", str(k)) for k in atype]
    treat = [re.sub("esophageal tissue ", "", str(k)) for k in treat]
    treat = [re.sub(" .*", "", str(k)) for k in treat]
    atype = treat
    atypes = ['U', 'Sh', 'ST', 'O']
    ahash = {'STW5':2, 'sham':1, 'untreated':0, 'omeprazole':3}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
        atypes = ['U', 'OMP']
        ahash = {'untreated':0, 'omeprazole':1}
    if (tn == 3):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
        atypes = ['U', 'OMP']
        ahash = {'untreated':0, 'omeprazole':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAbdelAziz2016 = getAbdelAziz2016


def getWickramasinghe2015(self, tn=1):
    self.prepareData("COV133")
    h = self.h
    atype = h.getSurvName('c incubated with')
    atypes = ['GLU', 'HMO', 'LAC']
    ahash = {'glucose-grown B. infantis':0,
            'HMO-grown B. infantis':1,
            'lactose-grown B. infantis':2}
    if (tn == 2):
        atypes = ['C', 'HMO']
        ahash = {'glucose-grown B. infantis':0,
                'HMO-grown B. infantis':1,
                'lactose-grown B. infantis':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWickramasinghe2015 = getWickramasinghe2015


def getHead2011(self, tn=1):
    self.prepareData("COV134")
    h = self.h
    atype = h.getSurvName('c cell type')
    atype = [re.sub("Caco2 cells exposed to ", "", str(k)) for k in atype]
    atype = [re.sub("human milk oligosacharides", "HMO", str(k)) for k in atype]
    atypes = ['C', 'HMO', 'HMO+Bi', 'Bi+Lac']
    ahash = {'B. infantis pre-grown on lactose.':3, 'HMO':1,
            'HMO and to B. infantis pre-grown on HMO.':2,
            'Caco2 cells with no B. infantis and no HMO':0}
    if (tn == 2):
        atypes = ['C', 'HMO']
        ahash = {'HMO':1, 'Caco2 cells with no B. infantis and no HMO':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHead2011 = getHead2011


def getChen2017(self, tn=1):
    self.prepareData("COV135")
    h = self.h
    atype = h.getSurvName('c infection status')
    atypes = ['C', 'Zika']
    ahash = {'mock':0, 'ZIKA-infected':1}
    if (tn == 2):
        atype = h.getSurvName('c src1')
        atypes = ['DMSO', 'AQ', 'HH']
        ahash = {'ZIKA infection, AQ drug':1,
                'ZIKA infection, HH drug':2,
                'ZIKA infection, DMSO vehicle':0}
    if (tn == 3):
        atype = h.getSurvName('c src1')
        atypes = ['DMSO', 'AQ']
        ahash = {'ZIKA infection, AQ drug':1,
                'ZIKA infection, DMSO vehicle':0}
    if (tn == 4):
        atype = h.getSurvName('c src1')
        atypes = ['DMSO', 'HH']
        ahash = {'ZIKA infection, HH drug':1,
                'ZIKA infection, DMSO vehicle':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChen2017 = getChen2017


def getRialdi2016(self, tn=1):
    self.prepareData("COV136")
    h = self.h
    atype = h.getSurvName('c infection')
    ahash = {'no infection':0, 'A/PR/8/34(\xce\x94NS1) Infection':1}
    ahash = asciiNorm(ahash)
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c treatment')
    atypes = ['C', 'iTop']
    ahash = {'Top1 siRNA':1, 'Control siRNA':0, 'no siRNA':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRialdi2016 = getRialdi2016


def getTheken2019(self, tn=1):
    self.prepareData("COV137")
    h = self.h
    atype = h.getSurvName('c drug treatment')
    ahash = {'ibuprofen sodium':1, 'Placebo':0}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c timepoint')
    ahash = {'baseline':0, 'post-surgery 1':1, 'post-surgery 2':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c response group')
    atypes = ['P', 'PR', 'CR']
    ahash = {'Partial responder':1, 'Placebo':0, 'Full responder':2}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 1 or tval[i] == 2
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = h.getSurvName('c drug treatment')
        atypes = ['P', 'IBU']
        ahash = {'ibuprofen sodium':1, 'Placebo':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTheken2019 = getTheken2019


def getFerretti2018(self, tn=1):
    self.prepareData("COV138")
    h = self.h
    atype = h.getSurvName('c Title')
    atype = [re.sub("H.* ileum (.*)_.*", "\\1", str(k)) for k in atype]
    atypes = ['control', 'IBU']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFerretti2018 = getFerretti2018


def getJabbari2015(self, tn=1):
    self.prepareData("COV139")
    h = self.h
    atype = h.getSurvName('c treatment')
    atypes = ['C', 'B']
    ahash = {'vehiclecontrol':0, 'baricitinib':1, 'vehicle control':0}
    if (tn == 2):
        atype = h.getSurvName('c src1')
        ahash = {'skin, topical baricitinib after disease establishment, week 0':0,
                'skin, topical baricitinib after disease establishment, week 12':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJabbari2015 = getJabbari2015


def getDrawnel2017(self, tn=1):
    self.prepareData("COV140")
    h = self.h
    atype = h.getSurvName('c cell type')
    atypes = ['CM']
    ahash = {'iPS-derived cardiomyocytes':0}
    if (tn == 2):
        atype = h.getSurvName('c compound name')
        atypes = ['C', 'R', 'Q']
        ahash = {'media_Ctrl':0, 'Resveratrol':1, 'Quinine-HCl-2H2O':2}
    if (tn == 3):
        atype = h.getSurvName('c compound name')
        atypes = ['C', 'Q']
        ahash = {'media_Ctrl':0, 'Quinine-HCl-2H2O':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDrawnel2017 = getDrawnel2017


def getTakeshita2019(self, tn=1):
    self.prepareData("COV141")
    h = self.h
    atype = h.getSurvName('c disease status')
    atypes = ['H', 'SN', 'Non', 'TCZ', 'MTX', 'IFX']
    ahash = {'RA TCZ treatment':3, 'RA non treatment':2, 'RA MTX treatment':4,
            'RA IFX treatment':5, 'healthy':0, 'RA synovial fluid':1}
    if (tn == 2):
        atypes = ['C', 'TCZ']
        ahash = {'RA non treatment':0, 'RA TCZ treatment':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTakeshita2019 = getTakeshita2019


def getNakamura2016(self, tn=1):
    self.prepareData("COV142")
    h = self.h
    atype = h.getSurvName('c sampling point')
    ahash = {'Before abatacept administration':0,
            'Before infliximab administration':1,
            'Before tocilizumab administration':2}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c clinical outcome')
    atypes = ['R', 'NR']
    ahash = {'Remission':0, 'Non-remission':1}
    if (tn == 2):
        atype = [atype[i] if dval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNakamura2016 = getNakamura2016


def getNishimoto2014(self, tn=1):
    self.prepareData("COV143")
    h = self.h
    atype = h.getSurvName('c treatment')
    atypes = ['bT', 'aT', 'bM', 'aM']
    ahash = {'after Tocilizumab/MRA':1, 'before methotrexate/MTX':2,
            'after methotrexate/MTX':3, 'before Tocilizumab/MRA':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNishimoto2014 = getNishimoto2014


def getGaertner2012(self, tn=1):
    self.prepareData("COV144")
    h = self.h
    atype = h.getSurvName('c ventricle')
    ahash = {'right':1, 'left':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c indication')
    atypes = ['NF', 'DCM', 'ARVC']
    ahash = {'Dilated cardiomyopathy':1,
            'Arrhythmogenic right ventricular cardiomyopathy':2,
            'Non-Failing':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = h.getSurvName('c ventricle')
        atype = [atype[i] if aval[i] == 2
                else None for i in range(len(atype))]
        atypes = ['left', 'right']
        ahash = {}
    if (tn == 3):
        atypes = ['NF', 'HF']
        ahash = {'Dilated cardiomyopathy':1,
                'Arrhythmogenic right ventricular cardiomyopathy':1,
                'Non-Failing':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGaertner2012 = getGaertner2012


def getvandenBerg2018(self, tn=1):
    self.prepareData("COV145")
    h = self.h
    atype = h.getSurvName('c tissue')
    ahash = {'Peripheral Blood Mononuclear Cells (PBMC)':0,
            'Whole Blood (WB)':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c sample day')
    atypes = ['D0', 'D40', 'D44', 'D47', 'D31', 'D37', 'D30']
    ahash = {}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['D0', 'D31', 'D44']
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getvandenBerg2018 = getvandenBerg2018


def getLemoine2018(self, tn=1):
    self.prepareData("COV146")
    h = self.h
    atype = h.getSurvName('c tissue')
    ahash = {'Peripheral blood':0, 'Cord Blood':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c stimulation')
    atypes = ['U', 'BCG', 'TLR', 'RSV']
    ahash = {'Unstimulated':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLemoine2018 = getLemoine2018


def getBerry2010(self, tn=1):
    self.prepareData("COV147")
    h = self.h
    atype = h.getSurvName('c src1')
    atype = [re.sub("Human ", "", str(k)) for k in atype]
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'whole':0, 'Whole':0, 'Neutrophils':2, 'CD8+':4,
            'Monocytes':1, 'CD4+':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c illness')
    atypes = ['C', 'LTB', 'PTB']
    ahash = {'LATENT TB':1, 'Control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = h.getSurvName('c bcg vaccinated')
        atypes = ['Yes', 'No']
        ahash = {}
    if (tn == 3):
        atype = h.getSurvName('c ethnicity')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'Caucasian':0, 'African American':1, 'Hispanic':3,
                'South Asian':2, 'Asian Other':2, 'Asian':2, 'White':0,
                'Native American':3, 'Black':2, 'Asian other':2,
                'Afican American':1, 'Caucasian/Asian':0, 'Other':3}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerry2010 = getBerry2010


def getMatsumiya2014(self, tn=1):
    self.prepareData("COV148")
    h = self.h
    atype = h.getSurvName('c time')
    ahash = {'0':0, '14':14, '2':2, '7':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c group')
    atypes = ['noBCG', 'BCG']
    ahash = {'A':0, 'C':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatsumiya2014 = getMatsumiya2014


def getLoxton2017(self, tn=1):
    self.prepareData("COV149")
    h = self.h
    atype = h.getSurvName('c timepoint')
    ahash = {'6':6, '0':0, '2':2, '26':26, '12':12, '18':18}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c group')
    atypes = ['VPM1002', 'BCG']
    ahash = {'VPM1002':0, 'BCG':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = h.getSurvName('c timepoint')
        atypes = ['0', '2-6-12', '18-26']
        ahash = {'2':1, '6':1, '12':1, '18':2, '26':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLoxton2017 = getLoxton2017


def getKaushal2015(self, tn=1):
    self.prepareData("COV150")
    h = self.h
    atype = h.getSurvName('c vaccinated with')
    atypes = ['MtbDsigH', 'BCG']
    ahash = {'aerosols with attenuated MtbDsigH mutant':0,
            'aerosols with BCG':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKaushal2015 = getKaushal2015


def getDarrah2020(self, tn=1):
    self.prepareData("COV151")
    h = self.h
    atype = h.getSurvName("c time after bcg vaccination")
    ahash = {'Week 13':13, 'Week 25':25}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName("c stimulation")
    ahash = {'Unstimulated cells':0, 'PPD-Stimulated Cells':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c vaccination route')
    atypes = ['noBCG', 'BCG']
    ahash = {'Intradermal':1, 'Naive':0, 'Aerosol':1, 'Intravenous':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 25 and gval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atype = [atype[i] if gval[i] == 0
                else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDarrah2020 = getDarrah2020


def getDarrah2020II(self, tn=1, ta=13, tb=0):
    self.prepareData("COV151.2")
    h = self.h
    atype = h.getSurvName("c Cell Type")
    ahash = {'Epithelial':0, 'Proliferating':1, 'Mac':2, 'Eo':3,
            'Neutrophils':4, 'Mast':5, 'B':6, 'T':7, 'NA':8}
    bval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName("c time after bcg vaccination")
    ahash = {'Week 13':13, 'Week 25':25}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName("c stimulation")
    ahash = {'Unstimulated cells':0, 'PPD-Stimulated Cells':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName('c vaccination route')
    atypes = ['noBCG', 'BCG']
    ahash = {'Intradermal':1, 'Naive':0, 'Aerosol':1, 'Intravenous':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if gval[i] == 0 and tval[i] == 25 and bval[i] == 8
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = gval
        phash = {6:1, 7:1}
        atype = [atype[i] if tval[i] == ta and bval[i] in phash and aval[i] == tb
                else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    if (tn == 4):
        atype = gval
        phash = {2:1, 5:1}
        atype = [atype[i] if tval[i] == ta and bval[i] in phash and aval[i] == tb
                else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDarrah2020II = getDarrah2020II


def getHatae2020(self, tn=1):
    self.prepareData("LU14")
    h = self.h
    atype = h.getSurvName("c treatment")
    atypes = ['Pre', 'Post']
    ahash = {'Pre-treatment':0, 'Post-treatment (nivolumab)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHatae2020 = getHatae2020


def getKadara2013(self, tn=1):
    self.prepareData("LU18")
    h = self.h
    atype = h.getSurvName("c airway site")
    ahash = {'Contralateral':4, 'Adjacent':3, 'Non-adjacent':2,
            'Main carina':1, 'NA':0}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = h.getSurvName("c time point")
    ahash = {'0':0, '12':12, '24':24, '36':36}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['0', '12', '24', '36']
    ahash = {}
    if (tn == 2):
        atype = [gval[i] if tval[i] == 12 or tval[i] == 36
                else None for i in range(len(atype))]
        atypes = ['A', 'Na', 'C', 'Mc']
        ahash = {3:0, 2:1, 1:3, 4:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKadara2013 = getKadara2013


def getNorenHooten2019(self, tn=1):
    self.prepareData("COV152")
    h = self.h
    atype = h.getSurvName("c age group")
    atypes = ['Y', 'O']
    ahash = {}
    if (tn == 2):
        atype = h.getSurvName("c race")
        atypes = ['W', 'AA']
        ahash = {'White':0, 'African American':1}
    if (tn == 3):
        atype = h.getSurvName("c poverty")
        atypes = ['A', 'B']
        ahash = {'Above':0, 'Below':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNorenHooten2019 = getNorenHooten2019


def getPrince2019(self, tn=1):
    self.prepareData("COV153")
    h = self.h
    atype = h.getSurvName("c race")
    atypes = ['W', 'AA']
    ahash = {'White':0, 'African American':1}
    if (tn == 2):
        atype = h.getSurvName("c frailty status")
        ahash = {'Non-Frail':0, 'Frail':1}
        atypes = ['NF', 'F']
    if (tn == 3):
        group = h.getSurvName("c frailty status")
        atype = [atype[i] if group[i] == 'Frail'
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPrince2019 = getPrince2019


def getDluzen2016(self, tn=1):
    self.prepareData("COV154")
    h = self.h
    atype = h.getSurvName("c race")
    atypes = ['W', 'AA']
    ahash = {'African American':1, 'Caucasian':0}
    if (tn == 2):
        atype = h.getSurvName("c src1")
        atypes = ['WN', 'WH', 'AAN', 'AAH']
        ahash = {'African American_hypertensives':3,
                'African American_normotensives':2,
                'white normotensives':0, 'white hypertensives':1,
                'white_hypertensives':1}
    if (tn == 3):
        atype = h.getSurvName("c src1")
        atypes = ['WN', 'WH']
        ahash = {'white normotensives':0, 'white hypertensives':1,
                'white_hypertensives':1}
    if (tn == 4):
        atype = h.getSurvName("c src1")
        atypes = ['AAN', 'AAH']
        ahash = {'African American_hypertensives':1,
                'African American_normotensives':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDluzen2016 = getDluzen2016


def getHubal2017(self, tn=1):
    self.prepareData("COV155")
    h = self.h
    group = h.getSurvName("c group")
    src1 = h.getSurvName("c src1")
    atype = [" ".join([str(group[i]), str(src1[i])])
            for i in range(len(group))]
    atypes = ['HB', 'PreDB', 'HP', 'PreDP']
    ahash = {'Healthy PBMC_baseline':0, 'Healthy PBMC_post-ex':2,
            'Prediabetic PBMC_baseline':1, 'Prediabetic PBMC_post-ex':3}
    if (tn == 2):
        atypes = ['HB', 'PreDB']
        ahash = {'Healthy PBMC_baseline':0, 'Prediabetic PBMC_baseline':1}
    if (tn == 3):
        atypes = ['B', 'Ex']
        ahash = {'Healthy PBMC_baseline':0, 'Healthy PBMC_post-ex':1,
                'Prediabetic PBMC_baseline':0, 'Prediabetic PBMC_post-ex':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHubal2017 = getHubal2017


def getSilva2018(self, tn=1):
    self.prepareData("COV156")
    h = self.h
    atype = h.getSurvName("c race")
    atypes = ['brown', 'white', 'black']
    ahash = {}
    if (tn == 2):
        atypes = ['white', 'black']
    if (tn == 3):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub(".*, ", "", str(k)) for k in atype]
        atype = [re.sub(" .*", "", str(k)) for k in atype]
        atypes = ['C', 'UnC']
        ahash = {'uncontrolled':1, 'controlled':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSilva2018 = getSilva2018


def getLee2014I(self, tn=1):
    self.prepareData("COV157")
    h = self.h
    atype = h.getSurvName("c ethnicity")
    atypes = ['W', 'AA', 'A', 'M']
    ahash = {'Caucasian':0, 'African American':1, 'Asian':2,
            'African-American':1, 'East Asian':2, 'MULTI-RACIAL':3}
    if (tn == 2):
        group = h.getSurvName("c stimulation")
        atype = [atype[i] if group[i] == 'Unstim'
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = h.getSurvName("c stimulation")
        atypes = ['Unstim', 'LPS', 'dNS1']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLee2014I = getLee2014I


def getLee2014II(self, tn=1):
    self.prepareData("COV157.2")
    h = self.h
    atype = h.getSurvName("c ethnicity")
    atypes = ['W', 'AA', 'A', 'M']
    ahash = {'Caucasian':0, 'African American':1, 'Asian':2,
            'African-American':1, 'East Asian':2, 'MULTI-RACIAL':3}
    if (tn == 2):
        group = h.getSurvName("c stimulation")
        atype = [atype[i] if group[i] == 'unstim'
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLee2014II = getLee2014II


def getHuang2015(self, tn=1):
    self.prepareData("LU17")
    atype = self.h.getSurvName("c study")
    ahash = {'ACCURACY STUDY':0, 'EQUIVALENCE STUDY':1, 'PRECISION STUDY':2,
            'SPECIFICITY STUDY':3, 'SENSITIVITY STUDY':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c risk")
    atypes = ['low', 'high']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 2
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHuang2015 = getHuang2015


def getFavreau2008(self, tn=1):
    self.prepareData("COV164")
    atype = self.h.getSurvName("c src1")
    time = [re.sub(".*, (.*)hou.*", "\\1", str(k)) for k in atype]
    ahash = {'48':48, '72':72, '24':24}
    tval = [ahash[i] if i in ahash else None for i in time]
    atype = [re.sub(".*, (.*) infe.*", "\\1", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'HCoV-OC43':1, 'mock':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 72
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFavreau2008 = getFavreau2008


def getMatualatupauw2019(self, tn=1):
    self.prepareData("MACV124")
    atype = self.h.getSurvName("c day")
    ahash = {'day29':29, 'day1':1}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c time point")
    ahash = {'120min':2, '360min':6, '0min':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c disease")
    atypes = ['H', 'MetS']
    ahash = {'Healty':0, 'Metabolic Syndrome':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if dval[i] == 1 and tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = [0, 2, 6]
        ahash = {}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = tval
        atypes = [0, 2, 6]
        ahash = {}
        atype = [atype[i] if aval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if tval[i] == 2
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatualatupauw2019 = getMatualatupauw2019


def getPaczkowskaAbdulsalam2020(self, tn=1):
    self.prepareData("MACV125")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*blood, (.*), sam.*", "\\1", str(k)) for k in atype]
    atypes = ['H, O', 'H, L', 'D, L', 'D, O']
    ahash = {}
    if (tn == 2):
        atypes = ['H, O', 'D, O']
    if (tn == 3):
        atypes = ['H, L', 'D, L']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPaczkowskaAbdulsalam2020 = getPaczkowskaAbdulsalam2020


def getRouet2018(self, tn=1):
    self.prepareData("HRT1")
    atype = self.h.getSurvName("c subject status")
    atypes = ['C', 'HN', 'HR']
    ahash = {'hypertensive patient with left ventricular remodeling':2,
            'hypertensive patient with normal left ventricular size':1,
             'control individual':0}
    if (tn == 2):
        atypes = ['C', 'HN']
        ahash = {'hypertensive patient with normal left ventricular size':1,
                 'control individual':0}
    if (tn == 3):
        atypes = ['C', 'HR']
        ahash = {'hypertensive patient with left ventricular remodeling':1,
                 'control individual':0}
    if (tn == 4):
        atypes = ['HN', 'HR']
        ahash = {'hypertensive patient with left ventricular remodeling':1,
                'hypertensive patient with normal left ventricular size':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRouet2018 = getRouet2018


def getEsser2018(self, tn=1):
    self.prepareData("HRT2")
    atype = self.h.getSurvName("c treatment")
    ahash = {'placebo':0, 'epicatechin (100mg/d)':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c before/after supplementation")
    atypes = ['B', 'A']
    ahash = {'after':1, 'before':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEsser2018 = getEsser2018


def getLi2016Mm(self, tn=1):
    self.prepareData("HRT3")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("D.*", "", str(k)) for k in atype]
    atypes = ['0', 'AngII7', 'AngII28']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2016Mm = getLi2016Mm


def getNelson2018Mm(self, tn=1):
    self.prepareData("HRT4")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("BP.*", "", str(k)) for k in atype]
    atypes = ['young', 'old']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub(".*(BP.).*", "\\1", str(k)) for k in atype]
        atypes = ['BPN', 'BPH']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNelson2018Mm = getNelson2018Mm


def getMahata2018(self, tn=1):
    self.prepareData("HRT5")
    atype = self.h.getSurvName("c genotype")
    atypes = ['Wt', 'Chga', 'CST']
    ahash = {'CST knockout':2, 'Wild-type':0, 'CHGA knockout':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMahata2018 = getMahata2018


def getPeng2020(self, tn=1):
    self.prepareData("HRT6")
    atype = self.h.getSurvName("c group")
    atypes = ['N', 'H', 'Q']
    ahash = {'Hypertension':1, 'QDG Treatment':2, 'Normal':0}
    if (tn == 2):
        atypes = ['N', 'H']
        ahash = {'Hypertension':1, 'Normal':0}
    if (tn == 3):
        atypes = ['H', 'QDG']
        ahash = {'Hypertension':0, 'QDG Treatment':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeng2020 = getPeng2020


def getSweeney2015(self, tn=1):
    self.prepareData("MACV130")
    atype = self.h.getSurvName("c disease")
    atypes = ['C', 'S', 'SS', 'SIRS']
    ahash = {'SepticShock':2, 'SIRS':3, 'Sepsis':1, 'Control':0}
    if (tn == 2):
        atypes = ['C', 'SIRS']
        ahash = {'SIRS':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSweeney2015 = getSweeney2015


def getMcHugh2015(self, tn=1):
    self.prepareData("MACV131.2")
    atype = self.h.getSurvName("c group")
    atypes = ['C', 'S']
    ahash = {'post-surgical':0, 'Sepsis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMcHugh2015 = getMcHugh2015


def getLance2017(self, tn=1):
    self.prepareData("LP3")
    atype = self.h.getSurvName('c Time               (Weeks)')
    ahash = {'0.142857':0, '1':1, '8':8, '3':3, '6':6, '2':2, '7':7,
            '10':10, '22':22, '23':23, '4':4, '5':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c BPD severity                             ')
    ahash = {'Severe':3, 'Mild':1, 'N/A':4, 'Moderate':2, 'None':0}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Treatment")
    atypes = ['CTRL', 'LPS']
    ahash = {'CTRL':0, 'LPS':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if dval[i] == 2
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c BPD severity                             ')
        atypes = ['N', 'S']
        ahash = {'Severe':1, 'Mild':0, 'N/A':0, 'Moderate':0, 'None':0}
        atype = [atype[i] if aval[i] == 0 and tval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLance2017 = getLance2017


def getCohen2007(self, tn=1):
    self.prepareData("LP4")
    atype = self.h.getSurvName("c Disease")
    atypes = ['nobpd', 'bpd']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCohen2007 = getCohen2007


def getCai2018(self, tn=1):
    self.prepareData("LP5")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['N', 'B']
    ahash = {'normal':0, 'Bronchopulmonary dysplasia (BPD)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCai2018 = getCai2018


def getDavidson2013(self, tn=1):
    self.prepareData("LP6")
    atype = self.h.getSurvName("c Treatment")
    atypes = ['L', 'LI', 'LD']
    ahash = {' LPS+IL-10':1, ' LPS+DEX':2, ' LPS':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDavidson2013 = getDavidson2013


def getPietrzyk2013(self, tn=1, ta=0):
    self.prepareData("LP7")
    atype = self.h.getSurvName('c Title')
    atype = [hu.re.sub(",.*", "", str(k)) for k in atype]
    ahash = {'Infant at timeB':1, 'Infant at timeA':0, 'Infant at timeC':2}
    tval = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName('c bronchopulmonary dysplasia (bpd) group')
    atypes = ['N', 'Mi', 'Mo', 'S', 'A']
    ahash = {'#N/A!':4, '0 (no BPD)':0, '1 (mild BPD)':1,
            '3 (severe BPD)':3, '2 (moderate BPD)':2}
    if (tn == 2):
        atype = self.h.getSurvName('c retinopathy of prematurity (rop) group')
        atypes = ['N', 'NT', 'LT', 'NA']
        ahash = {'#N/A!':3, '0 (no ROP)':0,
                '2 (ROP which need laser therapy)':2,
                '1 (ROP not requiring treatment)':1}
    if tn == 3:
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPietrzyk2013 = getPietrzyk2013


def getCho2012Mm(self, tn=1):
    self.prepareData("LP8")
    atype = self.h.getSurvName("c treatment")
    ahash = {'air':0, 'hyperoxia':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c genotype variation")
    atypes = ['W', 'Nrf']
    ahash = {'Nrf+/+':0, 'Nrf-/-':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCho2012Mm = getCho2012Mm


def getRittirsch2016(self, tn=1, tb=0):
    self.prepareData("MACV160")
    atype = self.h.getSurvName("c time point")
    atype = [re.sub(".*day (.*)[^0-9]", "\\1", str(k)) for k in atype]
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'0':0, '1':1, '2':2, '3':3, '5':5, '7':7, '10':10, '14':14,
            '21':21, '26':26, '28':28}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c subject subgroup")
    atypes = ['C', 'S']
    ahash = {'patients with systemic inflammation without infection':0,
            'patients with secondary sepsis after trauma':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRittirsch2016 = getRittirsch2016


def getLoren2014(self, tn=1):
    self.prepareData("MACV161")
    atype = self.h.getSurvName("c patient prognosis")
    atypes = ['G', 'M', 'P']
    ahash = {'Good prognosis':0,
            'Poor prognosis':2,
            'Moderate prognosis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLoren2014 = getLoren2014


def getStapels2018(self, tn=1):
    self.prepareData("MACV162")
    atype = self.h.getSurvName("c macrophage subtype")
    atypes = ['U', 'iHK', 'By', 'iG', 'iNG']
    ahash = {'Infected, with growing Salmonella':3,
            'Infected, with non-growing Salmonella':4,
            'Bystander':2,
            'Uninfected':0,
            'Infected, with non-viable Salmonella':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getStapels2018 = getStapels2018


def getStammet2012(self, tn=1):
    self.prepareData("MACV163")
    atype = self.h.getSurvName("c cpc")
    atypes = ['1', '2', '3', '4', '5']
    ahash = {}
    if (tn == 2):
        atypes = ['G', 'B']
        ahash = {'1':0, '2':0, '3':1, '4':1, '5':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getStammet2012 = getStammet2012


def getMoyer2010(self, tn=1):
    self.prepareData("MACV164")
    atype = self.h.getSurvName("c molecular group")
    atypes = ['U', 'I', 'F']
    ahash = {'Fibrosis':2, 'Inflammation':1, 'Unclassified':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMoyer2010 = getMoyer2010


def getHunter2010(self, tn=1):
    self.prepareData("MACV165")
    atype = self.h.getSurvName("c outcome at one year")
    atypes = ['P', 'E']
    ahash = {'persistent':0, 'extended':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHunter2010 = getHunter2010


def getKlapper2008(self, tn=1):
    self.prepareData("MACV166")
    atype = self.h.getSurvName("c Ann Arbour Stage")
    atypes = ['NA', 'I', 'II', 'III', 'IV']
    atype = self.h.getSurvName("c Gene expression")
    atype = [re.sub(".*: ", "", str(k)) for k in atype]
    atypes = ['N', 'I', 'M']
    ahash = {'intermediate':1, 'mBL':2, 'non-mBL':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKlapper2008 = getKlapper2008


def getGuthridge2020(self, tn=1):
    self.prepareData("MACV169")
    atype = self.h.getSurvName("c case/control")
    atypes = ['C', 'SLE']
    ahash = {'SLE Case':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuthridge2020 = getGuthridge2020


def getBohne2012I(self, tn=1):
    self.prepareData("MACV170")
    atype = self.h.getSurvName("c immunotolerance group")
    atypes = ['NT', 'TP', 'T']
    ahash = {'TOL POST':1, 'Non TOL':0, 'TOL':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBohne2012I = getBohne2012I


def getBohne2012II(self, tn=1):
    self.prepareData("MACV171")
    atype = self.h.getSurvName("c liver sample group")
    atypes = ['Cont', 'Non-TOL', 'HEPC', 'Non-TOL REJ', 'TOL', 'Cont-Tx', 'REJ']
    ahash = {}
    if (tn == 2):
        atypes = ['Cont-Tx', 'REJ']
    if (tn == 3):
        atypes = ['TOL', 'REJ']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBohne2012II = getBohne2012II


def getZander2011I(self, tn=1):
    self.prepareData("MACV172")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*: ", "", str(k)) for k in atype]
    atypes = ['C', 'BC']
    ahash = {'bronchial carcinoma':1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZander2011I = getZander2011I


def getZander2011II(self, tn=1):
    self.prepareData("MACV173")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(".*: ", "", str(k)) for k in atype]
    atype = [re.sub("X.*", "X", str(k)) for k in atype]
    atypes = ['C', 'BC', 'X']
    ahash = {'bronchial carcinoma':1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZander2011II = getZander2011II


def getPennycuick2020(self, tn=1):
    self.prepareData("COV165")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['PA', 'AA', 'PF', 'AF', 'PC', 'AC']
    ahash = {'Paediatric_Airway':0, 'Adult_Airway':1, 'Paediatric_FACS':2,
            'Adult_FACS':3, 'Paediatric_Cultured':4, 'Adult_Cultured':5}
    if (tn == 2):
        atypes = ['P', 'A']
        ahash = {'Paediatric_FACS':0, 'Adult_FACS':1,
                'Paediatric_Cultured':0, 'Adult_Cultured':1}
    if (tn == 3):
        atypes = ['P', 'A']
        ahash = {'Paediatric_Airway':0, 'Adult_Airway':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPennycuick2020 = getPennycuick2020


def getAulicino2018(self, tn=1):
    self.prepareData("MACV174")
    atype = self.h.getSurvName("c infection")
    atypes = ['C', 'LT', 'D']
    ahash = {'STM-D23580':2, 'STM-LT2':1, 'Mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAulicino2018 = getAulicino2018


def getAulicino2018Sc(self, tn=1):
    self.prepareData("MACV174.3")
    atype = self.h.getSurvName("c infection")
    atypes = ['M', 'LT', 'D']
    ahash = {'D23580':2, 'Mock':0, 'LT2':1, 'Blank':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAulicino2018Sc = getAulicino2018Sc


def getAvital2017(self, tn=1):
    self.prepareData("MACV175")
    atype = self.h.getSurvName("c agent")
    atypes = ['N', 'U', 'S', 'D']
    ahash = {'':3, 'none':0, 'Salmonella typhimurium SL1344':2, 'unexposed':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAvital2017 = getAvital2017


def getYe2018(self, tn=1):
    self.prepareData("MACV181")
    atype = self.h.getSurvName("c condition")
    atypes = ['B', 'InfA', 'IFNB']
    ahash = {'baseline':0, 'influenza stimulated':1, 'IFN-beta stimulated':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYe2018 = getYe2018


def getKash2017(self, tn=1):
    self.prepareData("COV166")
    atype = self.h.getSurvName("c virus strain")
    atypes = ['none', 'Ebola']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKash2017 = getKash2017


def getTang2019(self, tn=1):
    self.prepareData("COV167")
    sex = self.h.getSurvName("c Sex")
    ahash = {'f':0, 'm':1}
    sex = [ahash[i] if i in ahash else None for i in sex]
    age = self.h.getSurvName("c age")
    age = [int(age[i]) if i > 1 and age[i] != 'NA'
            else None for i in range(len(age))]
    atype = self.h.getSurvName("c severity")
    atypes = ['C', 'M', 'S']
    ahash = {'flu_mod':1, 'flu_svre':2, 'hlty_ctrl':0}
    if (tn == 2):
        atypes = ['M', 'S']
        ahash = {'flu_mod':0, 'flu_svre':1}
    if (tn == 3):
        atype = [atype[i] if age[i] > 50
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTang2019 = getTang2019


def getZhang2020(self, tn=1):
    self.prepareData("COV169")
    atype = self.h.getSurvName("c patient group")
    atypes = ['H', 'CoV']
    ahash = {'mild':1, 'severe':1, 'healthy control':0,
            'severe COVID-19 patient':1}
    if (tn == 2):
        atypes = ['H', 'M', 'S']
        ahash = {'mild':1, 'severe':2, 'healthy control':0,
                'severe COVID-19 patient':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2020 = getZhang2020


def getZhang2020Mac(self, tn=1):
    self.prepareData("COV169.2")
    atype = self.h.getSurvName("c patient group")
    atypes = ['H', 'CoV']
    ahash = {'mild':1, 'severe':1, 'healthy control':0,
            'severe COVID-19 patient':1}
    if (tn == 2):
        atypes = ['H', 'M', 'S']
        ahash = {'mild':1, 'severe':2, 'healthy control':0,
                'severe COVID-19 patient':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2020Mac = getZhang2020Mac


def getZhang2020CB(self, tn=1, tb=0):
    self.prepareData("COV169.4")
    atype = self.h.getSurvName("c Cell Type")
    atype = [re.sub("[_-][c1].*", "", str(k)) for k in atype]
    ahash = {'B':0, 'T':1, 'CD4_T':2, 'CD8_T':3, 'Natural_killer':4,
            'Macs_Monos_DCs':5, 'Epithelial':6}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c patient group")
    atypes = ['H', 'CoV']
    ahash = {'mild':1, 'severe':1, 'healthy control':0,
            'severe COVID-19 patient':1}
    if (tn == 2):
        atypes = ['H', 'M', 'S']
        ahash = {'mild':1, 'severe':2, 'healthy control':0,
                'severe COVID-19 patient':2}
    if (tn == 3):
        atypes = ['H', 'M', 'S']
        ahash = {'mild':1, 'severe':2, 'healthy control':0,
                'severe COVID-19 patient':2}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['H', 'CoV']
        ahash = {'mild':1, 'severe':1, 'healthy control':0,
                'severe COVID-19 patient':1}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2020CB = getZhang2020CB


def getZhang2020Epi(self, tn=1):
    self.prepareData("COV169.5")
    atype = self.h.getSurvName("c patient group")
    atypes = ['H', 'CoV']
    ahash = {'mild':1, 'severe':1, 'healthy control':0,
            'severe COVID-19 patient':1}
    if (tn == 2):
        atypes = ['H', 'M', 'S']
        ahash = {'mild':1, 'severe':2, 'healthy control':0,
                'severe COVID-19 patient':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2020Epi = getZhang2020Epi


def getZhang2020CD8(self, tn=1):
    self.prepareData("COV169.6")
    atype = self.h.getSurvName("c patient group")
    atypes = ['H', 'CoV']
    ahash = {'mild':1, 'severe':1, 'healthy control':0,
            'severe COVID-19 patient':1}
    if (tn == 2):
        atypes = ['H', 'M', 'S']
        ahash = {'mild':1, 'severe':2, 'healthy control':0,
                'severe COVID-19 patient':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2020CD8 = getZhang2020CD8


def getButler2011(self, tn=1, ta=0):
    self.prepareData("COV170")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'alveolar':0, 'airway':1}
    self.tissue = [ahash[i] if i in ahash else None for i in atype]
    v1 = self.h.getSurvName("c Age")
    v2 = self.h.getSurvName("c age")
    atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                            for i in range(len(v1))]
    self.age = [int(atype[i]) if i > 1 else None for i in range(len(atype))]
    v1 = self.h.getSurvName("c Ancestry")
    v2 = self.h.getSurvName("c ancestry")
    v3 = self.h.getSurvName("c Ethnic group")
    v4 = self.h.getSurvName("c ethnic group")
    atype = [ "".join([str(k) for k in [v1[i], v2[i], v3[i], v4[i]]])
                                    for i in range(len(v1))]
    ahash = {'African':1, 'hispanic':2, 'white':0,
            'black':1, 'European':0, 'Hispanic':2}
    self.race = [ahash[i] if i in ahash else None for i in atype]
    v1 = self.h.getSurvName("c Sex")
    v2 = self.h.getSurvName("c sex")
    self.sex = [ "".join([str(k) for k in [v1[i], v2[i]]])
                            for i in range(len(v1))]
    v1 = self.h.getSurvName("c smoking status")
    v2 = self.h.getSurvName("c Smoking Status")
    v3 = self.h.getSurvName("c Smoking status")
    atype = [ "".join([str(k) for k in [v1[i], v2[i], v3[i]]])
                            for i in range(len(v1))]
    atype = [re.sub(".*, ", "", str(k)) for k in atype]
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atype = [re.sub("non-smoker", "0", str(k)) for k in atype]
    self.packs = [float(atype[i]) if i > 1 else None for i in range(len(atype))]
    atype = [ "".join([str(k) for k in [v1[i], v2[i], v3[i]]])
                            for i in range(len(v1))]
    self.ss = [re.sub(",.*", "", str(k)) for k in atype]
    atype = self.ss
    atypes = ['NS', 'S']
    ahash = {'non-smoker':0, 'smoker':1}
    if (tn == 2):
        atype = ['Y' if self.age[i] is not None and 
                self.age[i] < 40 else 'O'
                for i in range(len(atype))]
        atype = [atype[i] if self.age[i] is not None
                else None for i in range(len(atype))]
        atypes = ['Y', 'O']
        ahash = {}
    if (tn == 3):
        atype = self.race
        atypes = ['W', 'B', 'H']
        ahash = {0:0,1:1, 2:2}
    if (tn == 4):
        atype = self.sex
        atypes = ['F', 'M']
        ahash = {}
    if (tn == 5):
        atype = self.race
        atypes = ['W', 'B']
        ahash = {0:0, 1:1}
    atype = [atype[i] if self.tissue[i] == ta
            else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getButler2011 = getButler2011


def getTilley(self, tn=1):
    if ("c ethnicity" in self.h.survhdrs):
        v1 = self.h.getSurvName("c ethnicity")
        v2 = self.h.getSurvName("c ethnic group")
        self.eg = [ "".join([str(k) for k in [v1[i], v2[i]]])
                for i in range(len(v1))]
    else:
        self.eg = self.h.getSurvName("c ethnic group")
    atype = self.h.getSurvName("c smoking status")
    atype = [re.sub(", .*", "", str(k)) for k in atype]
    atypes = ['NS', 'S']
    ahash = {'smoker':1, 'S':1, 'nonsmoker':0, 'NS':0, 'non-smoker':0}
    self.ss = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 6):
        atype = self.h.getSurvName("c copd status")
        atypes = ['H', 'COPD']
        ahash = {'':0, 'yes':1}
    if (tn == 5):
        atype = self.eg
        atypes = ['W', 'B']
        ahash = {'Afr':1, 'Eur':0, 'black':1, 'white':0}
    if (tn == 3):
        atype = self.eg
        atypes = ['His', 'W', 'B', 'As']
        ahash = {'Afr':2, 'Eur':1, 'black':2, 'hispanic':0, 'white':1,
                'asian':3}
    if (tn == 2):
        atype = self.h.getSurvName("c age")
        self.age = [int(atype[i]) if i > 1 and atype[i] != ''
                else None for i in range(len(atype))]
        atype = ['Y' if  self.age[i] is not None and self.age[i] < 40
                else 'O' for i in range(len(atype))]
        atype = [atype[i] if self.age[i] is not None
                else None for i in range(len(atype))]
        atypes = ['Y', 'O']
        ahash = {}
    if (tn == 4):
        atype = self.h.getSurvName("c sex")
        atypes = ['F', 'M']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTilley = getTilley

def getTilley2016(self, tn=1):
    self.prepareData("COV171")
    self.getTilley(tn)
bone.IBDAnalysis.getTilley2016 = getTilley2016

def getWang2010(self, tn=1):
    atype = self.h.getSurvName("c Age")
    if ("c age" in self.h.survhdrs):
        v1 = atype
        v2 = self.h.getSurvName("c age")
        atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    self.age = [int(atype[i]) if i > 1 and atype[i] != ''
            else None for i in range(len(atype))]
    atype = self.h.getSurvName("c Ethnic group")
    if ("c ethnic group" in self.h.survhdrs):
        v1 = atype
        v2 = self.h.getSurvName("c ethnic group")
        atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    if ("c Ethnicity" in self.h.survhdrs):
        v1 = atype
        v2 = self.h.getSurvName("c Ethnicity")
        atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    ahash = {'white':1, 'hispanic':0, 'black/hispanic':2,
            'black':2, 'hispnaic':0, 'asian':3}
    self.eg = [ahash[i] if i in ahash else None for i in atype]
    self.sex = self.h.getSurvName("c Sex")
    if ("c sex" in self.h.survhdrs):
        v1 = self.sex
        v2 = self.h.getSurvName("c sex")
        self.sex = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    if ("c Gender" in self.h.survhdrs):
        v1 = self.sex
        v2 = self.h.getSurvName("c Gender")
        ahash = {'':'', 'Male':'M', 'male':'M'}
        v2 = [ahash[i] if i in ahash else None for i in v2]
        self.sex = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    atype = ["" for k in atype]
    if "c Smoking Status" in self.h.survhdrs:
        atype = self.h.getSurvName("c Smoking Status")
    if "c smoking status" in self.h.survhdrs:
        v1 = atype
        v2 = self.h.getSurvName("c smoking status")
        atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    if "c Smoking status" in self.h.survhdrs:
        v1 = atype
        v2 = self.h.getSurvName("c Smoking status")
        atype = [ "".join([str(k) for k in [v1[i], v2[i]]])
                                for i in range(len(v1))]
    self.ss = atype
    atype = [re.sub(".*, ", "", str(k)) for k in atype]
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atype = [re.sub("non-smoker", "0", str(k)) for k in atype]
    self.packs = [float(atype[i]) if i > 1 and atype[i] != ''
            else None for i in range(len(atype))]
    atype = self.ss
    self.ss = [re.sub(",.*", "", str(k)) for k in atype]
    atype = self.ss
    atypes = ['NS', 'S']
    ahash = {'non-smoker':0, 'smoker':1}
    if (tn == 5):
        atype = self.eg
        atypes = ['W', 'B']
        ahash = {1:0, 2:1}
    if (tn == 4):
        atype = self.sex
        atypes = ['F', 'M']
        ahash = {}
    if (tn == 3):
        atype = self.eg
        atypes = ['His', 'W', 'B', 'As']
        ahash = {0:0, 1:1, 2:2, 3:3}
    if (tn == 2):
        atype = ['Y' if self.age[i] is not None and self.age[i] < 40
                else 'O' for i in range(len(atype))]
        atype = [atype[i] if self.age[i] is not None
                else None for i in range(len(atype))]
        atypes = ['Y', 'O']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2010 = getWang2010


def getWang2010I(self, tn=1):
    self.prepareData("COV172")
    self.getWang2010(tn)
bone.IBDAnalysis.getWang2010I = getWang2010I
def getWang2010II(self, tn=1):
    self.prepareData("COV173")
    self.getWang2010(tn)
bone.IBDAnalysis.getWang2010II = getWang2010II
def getTilley2011(self, tn=1):
    self.prepareData("COV174")
    self.getTilley(tn)
bone.IBDAnalysis.getTilley2011 = getTilley2011
def getShaykhiev2011(self, tn=1):
    self.prepareData("COV175")
    self.getWang2010(tn)
bone.IBDAnalysis.getShaykhiev2011 = getShaykhiev2011
def getStruloviciBarel2010(self, tn=1):
    self.prepareData("COV176")
    self.getWang2010(tn)
bone.IBDAnalysis.getStruloviciBarel2010 = getStruloviciBarel2010
def getTuretz2009(self, tn=1):
    self.prepareData("COV177")
    self.getTilley(tn)
bone.IBDAnalysis.getTuretz2009 = getTuretz2009
def getCarolan2008(self, tn=1):
    self.prepareData("COV178")
    self.getWang2010(tn)
bone.IBDAnalysis.getCarolan2008 = getCarolan2008
def getTilley2009(self, tn=1):
    self.prepareData("COV179")
    self.getWang2010(tn)
bone.IBDAnalysis.getTilley2009 = getTilley2009
def getCarolan2006I(self, tn=1):
    self.prepareData("COV180")
    self.getWang2010(tn)
bone.IBDAnalysis.getCarolan2006I = getCarolan2006I
def getCarolan2006II(self, tn=1):
    self.prepareData("COV181")
    self.getWang2010(tn)
bone.IBDAnalysis.getCarolan2006II = getCarolan2006II
def getCarolan2006III(self, tn=1):
    self.prepareData("COV182")
    self.getWang2010(tn)
bone.IBDAnalysis.getCarolan2006III = getCarolan2006III

def getAlmansa2012(self, tn=1):
    self.prepareData("MACV151")
    atype = self.h.getSurvName('c Characteristics[DiseaseState]')
    atypes = ['H', 'COPD']
    ahash = {'critical chronic obstructive pulmonary disease':1,
            'normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAlmansa2012 = getAlmansa2012


def getBigler2017(self, tn=1):
    self.prepareData("COV162")
    atype = self.h.getSurvName('c cohort')
    atypes = ['H', 'Asthma']
    ahash = {'Healthy, non-smoking':0, 'Severe asthma, non-smoking':1,
            'Severe asthma, smoking':1, 'Moderate asthma, non-smoking':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['NS', 'S']
        ahash = {'Healthy, non-smoking':0, 'Severe asthma, non-smoking':0,
                'Severe asthma, smoking':1, 'Moderate asthma, non-smoking':0}
    if (tn == 3):
        atype = self.h.getSurvName('c gender')
        atypes = ['F', 'M']
        ahash = {'male':1, 'female':0}
    if (tn == 4):
        atype = self.h.getSurvName('c race')
        atypes = ['W', 'AA']
        ahash = {'white_caucasian':0, 'black_african':1}
    if (tn == 5):
        atype = self.h.getSurvName('c race')
        atypes = ['W', 'AA', 'A', 'O']
        ahash = {'white_caucasian':0, 'south_asian':2, 'other':3,
                'black_african':1, 'arabic_north_heritage':3,
                'south_east_asian':2, 'multiple_races':3, 'east_asian':2,
                'central_asian':2}
        #atype = [atype[i] if aval[i] == 0
        #        else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBigler2017 = getBigler2017


def getKlenerman2020(self, tn=1):
    self.prepareData("COV183")
    atype = self.h.getSurvName('c cirrhosis present')
    atypes = ['C', 'Cir']
    ahash = {'Yes':1, 'No':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKlenerman2020 = getKlenerman2020


def getWyler2020(self, tn=1, ta = 0, tb = 4):
    self.prepareData("COV185")
    atype = self.h.getSurvName('c cell line')
    ahash = {'Caco2':0, 'Calu3':1, 'H1299':2}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c molecule subtype")
    ahash = {'polyA RNA extracted from whole cells':0,
            'total RNA extracted from whole cells':1}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c time point')
    ahash = {'4h':4, '12h':12, '24h':24, '8h':8, '36h':36}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c infection')
    atypes = ['U', 'M', 'CoV1', 'CoV2']
    ahash = {'SARS-CoV-1':2, 'SARS-CoV-2':3, 'mock':1, 'untreated':0}
    if (tn == 2):
        atypes = ['C', 'CoV2']
        ahash = {'SARS-CoV-2':1, 'mock':0, 'untreated':0}
        atype = [atype[i] if gval[i] == ta and tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'CoV2']
        ahash = {'SARS-CoV-2':1, 'mock':0, 'untreated':0}
        atype = [atype[i] if gval[i] == ta
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['C', 'CoV2']
        ahash = {'SARS-CoV-2':1, 'mock':0, 'untreated':0}
        atype = [atype[i] if gval[i] == 1 and mval[i] == 0
                and tval[i] >= 12
                else None for i in range(len(atype))]
    if (tn == 5):
        atypes = ['C', 'CoV2']
        ahash = {'SARS-CoV-2':1, 'mock':0, 'untreated':0}
        atype = [atype[i] if gval[i] == 1 and mval[i] == 0
                and tval[i] <= 12
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWyler2020 = getWyler2020


def getLamers2020(self, tn=1):
    self.prepareData("COV186")
    atype = self.h.getSurvName('c desc')
    atype = [re.sub("Bulk.*d, in ", "", str(k)) for k in atype]
    medium = [re.sub("ion.*", "ion", str(k)) for k in atype]
    ahash = {'differentiation':0, 'expansion':1}
    mval = [ahash[i] if i in ahash else None for i in medium]
    atype = [re.sub("Biological.*", "", str(k)) for k in atype]
    time = [k.split(" ")[2] if len(k.split()) > 2 else '' for k in atype]
    ahash = {'':0, '24':24, '60':60, '72':72}
    tval = [ahash[i] if i in ahash else None for i in time]
    time = [k.split(" ")[2] if len(k.split()) > 2 else '' for k in time]
    atype = [k.split(" ")[5] if len(k.split()) > 5 else '' for k in atype]
    atypes = ['U', 'CoV1', 'CoV2']
    ahash = {'':0, 'SARS-CoV2':2, 'SARS-CoV':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['U', 'CoV2']
        ahash = {'':0, 'SARS-CoV2':1}
    if (tn == 3):
        atypes = ['U', 'CoV1']
        ahash = {'':0, 'SARS-CoV':1}
    if (tn == 4):
        atypes = ['U', 'CoV2']
        ahash = {'':0, 'SARS-CoV2':1}
        atype = [None if tval[i] != 72 and aval[i] == 2
                else atype[i] for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLamers2020 = getLamers2020


def getVerhoeven2014(self, tn=1):
    self.prepareData("COV187")
    atype = self.h.getSurvName('c tissue')
    ahash = {'lung mucosa':0, 'colonic mucosa':1, 'jejunal mucosa':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c desc')
    atype = [re.sub("Gene.*from ", "", str(k)) for k in atype]
    atype = [re.sub(" mac.*", "", str(k)) for k in atype]
    atype = [k.split(" ")[-1] for k in atype]
    atypes = ['H', 'U', 'T']
    ahash = {'treated':2, 'untreated':1, 'healthy':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVerhoeven2014 = getVerhoeven2014


def getHosmillo2020(self, tn=1):
    self.prepareData("COV188")
    atype = self.h.getSurvName('c tissue')
    ahash = {'ileum organoid TI006':0, 'ileum organoid TI365':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['M', 'NoV', 'uNoV']
    ahash = {'HuNoV GII.4 infection 48h':1, 'Mock 48h':0,
            'UV-treated HuNoV GII.4 48h':2}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['M', 'NoV']
        ahash = {'HuNoV GII.4 infection 48h':1, 'Mock 48h':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHosmillo2020 = getHosmillo2020


def getCuadras2002I(self, tn=1):
    self.prepareData("COV189")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'T1', 'T6', 'T12', 'T24']
    ahash = {'Tc 1h':1, 'Tc control':0, 'Tc 6h':2, 'Tc 24h':4, 'Tc 12h':3}
    if (tn == 2):
        atypes = ['C', 'I']
        ahash = {'Tc 1h':0, 'Tc control':0, 'Tc 6h':1, 'Tc 24h':1, 'Tc 12h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCuadras2002I = getCuadras2002I


def getCuadras2002II(self, tn=1):
    self.prepareData("COV189.2")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" .*[^0-9]([0-9]+)h", " \\1h", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'Infection 1h':0, 'Control 16h':0, 'Control 1h':0, 'Infection 16h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCuadras2002II = getCuadras2002II


def getMedigeshi2020(self, tn=1):
    self.prepareData("COV191")
    atype = self.h.getSurvName('c infection status')
    atypes = ['M', 'JEV']
    ahash = {'Japanese encephalitis virus':1, 'mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMedigeshi2020 = getMedigeshi2020


def getHu2013(self, tn=1):
    self.prepareData("MACV146")
    atype = self.h.getSurvName('c group')
    atypes = ['Ctrl', 'RNA-virus','DNA-virus']
    ahash = {'Adenovirus':2, 'Virus-negative Control':0,'HHV6':2,
            'Enterovirus':1,'Rhinovirus':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'I']
        ahash = {'Adenovirus':1, 'Virus-negative Control':0,'HHV6':1,
                'Enterovirus':1,'Rhinovirus':1}
    if (tn == 3):
        atype = self.h.getSurvName("c ethnicity")
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'White':0, 'Black':1, 'Other':3}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'Rhinovirus':1}
    if (tn == 5):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'Enterovirus':1}
    if (tn == 6):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'Adenovirus':1}
    if (tn == 7):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'HHV6':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHu2013 = getHu2013


def getAltman2019I(self, tn=1):
    self.prepareData("COV192")
    atype = self.h.getSurvName('c analysis.visit')
    ahash = {'Visit 0':0, 'Visit 1b':1, 'Visit 1a':2, 'Visit 2a':3, 'Visit 2b':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c csteroid.start.relative.to.visit')
    ahash = {'Before':0, 'After':1, '':2}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c case.or.control.status.original')
    ahash = {'':2, 'Control':0, 'Case':1}
    bval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c viral.type.at.visit')
    ahash = {'Viral':1, 'Non-viral':0, '':2, 'NA':2} 
    mval = [ahash[i] if i in ahash else None for i in atype]
    atypes = [ 'Non-viral', 'Viral']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c GLI_RACE')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'Black':1, 'Other':3, 'White':0, 'SE Asian':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAltman2019I = getAltman2019I


def getAltman2019II(self, tn=1):
    self.prepareData("COV192.2")
    atype = self.h.getSurvName('c analysis.visit')
    ahash = {'Visit 0':0, 'Visit 1b':1, 'Visit 1a':2, 'Visit 2a':3, 'Visit 2b':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c csteroid.start.relative.to.visit')
    ahash = {'Before':0, 'After':1, '':2}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c case.or.control.status.v0')
    ahash = {'':2, 'Control':0, 'Case':1}
    bval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c GLI_RACE')
    ahash = {'Black':1, 'Other':3, 'White':0, 'SE Asian':2}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c viral.type.at.visit')
    ahash = {'Viral':1, 'Non-viral':0, '':2, 'NA':2} 
    mval = [ahash[i] if i in ahash else None for i in atype]
    atypes = [ 'Non-viral', 'Viral']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 0 and sval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c GLI_RACE')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'Black':1, 'Other':3, 'White':0, 'SE Asian':2}
        atype = [atype[i] if sval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAltman2019II = getAltman2019II


def getHealy2020(self, tn=1):
    self.prepareData("COV193")
    atype = self.h.getSurvName('c Title')
    tissue = [re.sub("-.*", "", str(k)) for k in atype]
    ahash = {'Macrophage':0, 'Microglia':1}
    tval = [ahash[i] if i in ahash else None for i in tissue]
    atype = self.h.getSurvName('c treatment')
    atypes = ['Ctrl', 'VitD']
    ahash = {'':0, '100nM calcitriol':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHealy2020 = getHealy2020


def getCasella2019(self, tn=1):
    self.prepareData("COV194")
    atype = self.h.getSurvName('c desc')
    cells = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'WI-38':0, 'IMR-90':1, 'HAECs':2, 'HUVEC':3}
    tval = [ahash[i] if i in ahash else None for i in cells]
    atype = [re.sub("[^ ]*\s+(.*)", "\\1", str(k)) for k in atype]
    atypes = ['Ctrl', 'Sen']
    ahash = {'doxorubicin-induced senescence':1,
            'proliferating control':0, 'control':0,
            'replicative senescence':1, 'IR-induced senescence':1,
            'empty vector control':0, 'oncogene-induced senescence':1,
            'IR treatment':1}
    if (tn == 2):
        atype = self.h.getSurvName('c desc')
        atypes = ['C', 'S']
        ahash = {'WI-38 doxorubicin-induced senescence':1,
                'WI-38  proliferating control':0,
                'WI-38 control':0,
                'WI-38  replicative senescence':1,
                'WI-38 IR-induced senescence':1,
                'WI-38 empty vector control':0,
                'WI-38 oncogene-induced senescence':1}
    if (tn == 3):
        atype = self.h.getSurvName('c desc')
        atypes = ['PC', 'IR', 'S']
        ahash = {'IMR-90 proliferating control':0,
                'IMR-90 IR treatment':1,
                'IMR-90  replicative senescence':2}
    if (tn == 4):
        atype = self.h.getSurvName('c desc')
        atypes = ['C', 'S']
        ahash = {'HAECs IR-induced senescence':1,
                'HAECs control':0}
    if (tn == 5):
        atype = self.h.getSurvName('c desc')
        atypes = ['C', 'S']
        ahash = {'HUVEC control':0,
                'HUVEC IR-induced senescence':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCasella2019 = getCasella2019


def getNurminen2015(self, tn=1):
    self.prepareData("COV195")
    atype = self.h.getSurvName('c treatment')
    atypes = ['Ctrl', 'VitD']
    ahash = {'vehicle':0, '1,25D':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNurminen2015 = getNurminen2015


def getSalamon2014(self, tn=1):
    self.prepareData("COV196")
    atype = self.h.getSurvName('c infection')
    ahash = {'M. tuberculosis H37Rv':1, 'uninfected':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['Ctrl', 'VitD']
    ahash = {'100 nM of 1α,25-dihydroxyvitamin D3':1, 'untreated':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSalamon2014 = getSalamon2014


def getHeikkinen2011(self, tn=1):
    self.prepareData("COV197")
    atype = self.h.getSurvName('c treatment')
    atypes = ['Ctrl', 'VitD']
    ahash = {'vehicle':0, '1,25D':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHeikkinen2011 = getHeikkinen2011


def getCosta2019(self, tn=1):
    self.prepareData("COV198")
    atype = self.h.getSurvName('c desc')
    atype = [re.sub("G.*with a ", "", str(k)) for k in atype]
    atype = [re.sub("VDR.*with ", "", str(k)) for k in atype]
    atypes = ['Ctrl', 'VitD', 'Mut', 'MutVitD']
    ahash = {'wild-type 1,25-vitamin D':1,
            'wild-type vehicle':0,
            'null p.Arg30* vehicle':2,
            'null p.Arg30* 1,25-vitamin D':3}
    if (tn == 2):
        atypes = ['Ctrl', 'VitD']
        ahash = {'wild-type 1,25-vitamin D':1, 'wild-type vehicle':0}
    if (tn == 3):
        atypes = ['Mut', 'MutVitD']
        ahash = {'null p.Arg30* vehicle':0, 'null p.Arg30* 1,25-vitamin D':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCosta2019 = getCosta2019


def getWang2011(self, tn=1):
    self.prepareData("COV199")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'VitD', 'tC', 'tVitD']
    ahash = {'No-Testosterone':0,
            '5nM-Testosterone-with100nM-VitD':3,
            'No-Testosterone-withVitD':1,
            '5nM-Testosterone':2}
    if (tn == 2):
        atypes = ['C', 'VitD']
        ahash = {'No-Testosterone':0,
                'No-Testosterone-withVitD':1}
    if (tn == 3):
        atypes = ['tC', 'tVitD']
        ahash = {'5nM-Testosterone-with100nM-VitD':1, '5nM-Testosterone':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2011 = getWang2011


def getCosta2009(self, tn=1):
    self.prepareData("COV200")
    atype = self.h.getSurvName('c Title')
    atypes = ['DR1', 'DR2']
    ahash = {'MCF7 VDr (2) X MCF7 (2)':0, 'MCF7 DR X MCF7':0,
            'MCF7 VDr X MCF7':0, 'MCF7 DRA (2) X MCF7 P38 (2)':1,
            'MCF7 DRA X MCF7 P38':1, 'MCF7 D3RE (2) X MCF7 P38 (2)':1,
            'MCF7 D3RE X MCF7 P38':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCosta2009 = getCosta2009


def getMartinezSena2020(self, tn=1):
    self.prepareData("COV201")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'C+VDR', 'VitD']
    ahash = {'control Ad-C + Vehicle':0, 'Ad-VDR + Vehicle':1, 'Ad-VDR + 10nMVitD':2}
    if (tn == 2):
        atypes = ['C', 'VitD']
        ahash = {'Ad-VDR + Vehicle':0, 'Ad-VDR + 10nMVitD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinezSena2020 = getMartinezSena2020


def getFerrerMayorga(self, tn=1):
    self.prepareData("COV202")
    atype = self.h.getSurvName('c tissue')
    ahash = {'Colon normal fibroblasts':0, 'Colon tumor fibroblasts':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'VitD']
    ahash = {'Vehicle':0, '1,25(OH)2D3':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFerrerMayorga = getFerrerMayorga


def getLiu2015(self, tn=1):
    self.prepareData("COV203")
    atype = self.h.getSurvName('c cell line')
    ahash = {'LNCaP':0, 'MCF-7':1, '':2} #Prostasphere
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'VitD']
    ahash = {'control for calcitriol':0, '100nM calcitriol':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 2
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiu2015 = getLiu2015


def getOda2015Mm(self, tn=1):
    self.prepareData("COV204")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'VDR_KO']
    ahash = {'VDR Control (CON) wounded skin':0,
            'VDR knockout (KO) wounded skin':1,
            'VDR knockout (KO) keratinocytes':1,
            'Control (CON) keratinocytes':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOda2015Mm = getOda2015Mm


def getHangelbroek2019(self, tn=1):
    self.prepareData("COV205")
    atype = self.h.getSurvName('c time of sampling')
    ahash = {'after intervention':1, 'before intervention (baseline)':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c intervention group')
    atypes = ['C', 'VitD']
    ahash = {'Placebo':0, '25-hydroxycholecalciferol (25(OH)D3)':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['C', 'VitD']
        ahash = {0:0, 1:1}
        atype = [tval[i] if aval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHangelbroek2019 = getHangelbroek2019


def getmeugnier2018rat(self, tn=1):
    self.prepareData("COV206")
    atype = self.h.getSurvName('c diet')
    atypes = ['C', 'VitDd', 'VitD']
    ahash = {'control diet (1000UI/Kg vit D)':0,
            'Vitamin D depletion diet':1,
            'depletion and repletion of vitamin D in the diet':2}
    if (tn == 2):
        atypes = ['C', 'VitDd']
        ahash = {'control diet (1000UI/Kg vit D)':0,
                'Vitamin D depletion diet':1}
    if (tn == 3):
        atypes = ['C', 'VitD']
        ahash = {'control diet (1000UI/Kg vit D)':0,
                'depletion and repletion of vitamin D in the diet':1}
    if (tn == 4):
        atypes = ['C', 'VitD']
        ahash = {'Vitamin D depletion diet':0,
                'depletion and repletion of vitamin D in the diet':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getmeugnier2018rat = getmeugnier2018rat


def getCummings2017(self, tn=1):
    self.prepareData("COV207")
    atype = self.h.getSurvName('c tissue')
    ahash = {"Barrett's esophagus segment":1, 'Normal esophageal squamous mucosa':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c arm')
    ahash = {'Arm A':0, 'Arm B':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c timepoint')
    atypes = ['C', 'VitD']
    ahash = {'T1':1, 'T0':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 and gval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 0 and gval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if tval[i] == 1 and gval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if tval[i] == 1 and gval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCummings2017 = getCummings2017


def getBae2011rat(self, tn=1):
    self.prepareData("COV208")
    atype = self.h.getSurvName('c treatment')
    atypes = ['U', 'EV', 'VitD', 'V', 'EP']
    ahash = {'untreated':0, 'enalapril and vehicle':1, 'paricalcitol':2,
            'vehicle':3, 'enalapril and paricalcitol':4}
    if (tn == 2):
        atypes = ['V', 'VitD']
        ahash = {'paricalcitol':1, 'vehicle':0}
    if (tn == 3):
        atypes = ['V', 'VitD']
        ahash = {'enalapril and vehicle':0, 'enalapril and paricalcitol':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBae2011rat = getBae2011rat


def getDankers2019(self, tn=1):
    self.prepareData("COV209")
    atype = self.h.getSurvName('c treatment')
    atypes = ['V', 'VitD']
    ahash = {'VitD':1, 'vehicle':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDankers2019 = getDankers2019


def getColdren2006II(self, tn=1):
    self.prepareData("COV211")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['A', 'B']
    ahash = {'H358':0, 'Calu3':0, 'A427':0, 'A549':0, 'H1703':1, 'H1299':1, 'HCC827':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getColdren2006II = getColdren2006II


def getLoboda2011II(self, tn=1):
    self.prepareData("COV213")
    atype = self.h.getSurvName('c cell line')
    atypes = ['A', 'B']
    ahash = {'H358':0, 'CALU3':0, 'A427':0, 'A549':0, 'H1703':1, 'H1299':1, 'HCC827':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLoboda2011II = getLoboda2011II


def getByers2013(self, tn=1):
    self.prepareData("COV215")
    atype = self.h.getSurvName('c cell line')
    atypes = ['A', 'B']
    ahash = {'H358':0, 'Calu-3':0, 'A427':0, 'A549':0, 'H1703':1, 'H1299':1, 'HCC827':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getByers2013 = getByers2013


def getAstraZeneca2014(self, tn=1):
    self.prepareData("COV218")
    atype = self.h.getSurvName('c cell line')
    atypes = ['A', 'B']
    ahash = {'H358':0, 'CALU3':0, 'A427':0, 'A549':0, 'CACO2':1, 'H1299':1, 'HCC827':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAstraZeneca2014 = getAstraZeneca2014


def getTan2019(self, tn=1):
    self.prepareData("COV219")
    atype = self.h.getSurvName('c cell type')
    ahash = {'CD3-CD56+ NK cells':0, 'Whole PBMC':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['U', 'I']
    ahash = {'12 hour exposure to influenza-infected A549 cells':1,
            '12 hour exposure to non-infected A549 cells':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTan2019 = getTan2019


def getBoeijen2019(self, tn=1):
    self.prepareData("COV220")
    atype = self.h.getSurvName('c hbv clinical phase')
    atypes = ['C', 'HCV', 'HIV', 'M0', 'M1', 'M2', 'Hep']
    ahash = {'HCV':1, 'HIV-1':2, 'Healthy':0, 'Immune Active':4,
            'HBeAg- Hepatitis':6, 'Immune Control':3, 'Immune Tolerant':5}
    if (tn == 2):
        atypes = ['C', 'HCV']
        ahash = {'HCV':1, 'Healthy':0}
    if (tn == 3):
        atypes = ['C', 'HIV']
        ahash = {'HIV-1':1, 'Healthy':0}
    if (tn == 4):
        atypes = ['C', 'HBV']
        ahash = {'Healthy':0, 'Immune Active':1, 'HBeAg- Hepatitis':1,
                'Immune Control':1, 'Immune Tolerant':1}
    if (tn == 5):
        atypes = ['C', 'Hep']
        ahash = {'Healthy':0, 'HBeAg- Hepatitis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBoeijen2019 = getBoeijen2019


def getWagstaffe2019(self, tn=1):
    self.prepareData("COV221")
    atype = self.h.getSurvName('c cell fraction')
    ahash = {'NK cell depleted':0, 'NK cell enriched':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c time point')
    atypes = ['B', 'pV']
    ahash = {'30 days post-vaccination':1, 'Baseline (pre-vaccination)':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWagstaffe2019 = getWagstaffe2019


def getSantosa2020(self, tn=1):
    self.prepareData("COV222")
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['Wt', 'Mut']
    ahash = {'LDHA-deficient':1, 'wild-type or LDHA-deficient':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSantosa2020 = getSantosa2020


def getTing2020(self, tn=1, tb = 0):
    self.prepareData("COV223")
    atype = self.h.getSurvName('c src1')
    ahash = {'lung':0, 'heart':1, 'jejunum':2, 'liver':3, 'kidney':4, 'bowel':5,
            'fat':6, 'skin':7, 'marrow':8}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c sample case')
    atypes = ['C', 'I']
    ahash = {'1':1, '2':1, '3':1, '4':1, '5':1, 'Control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [tval[i] if  tval[i] is not None
                else None for i in range(len(atype))]
        atype = [-1 if  aval[i] is not None and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTing2020 = getTing2020


def getHoek2015(self, tn=1, tb = 0, tc=1):
    self.prepareData("COV224")
    atype = self.h.getSurvName('c src1')
    ahash = {'PBMC':0, 'myeloid DC':1, 'Monocytes':2, 'Neutrophils':3,
            'B cells':4, 'T cells':5, 'NK cells':6}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c time')
    ahash = {'0 d':0, '1 d':1, '3 d':3, '7 d':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['0 d', '1 d', '3 d', '7 d']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
        atypes = ['C', 'pV']
        ahash = {'0 d':0, '1 d':1, '3 d':1, '7 d':1}
    if (tn == 4):
        atype = [atype[i] if sval[i] == tb and
                (tval[i] == 0 or tval[i] == tc)
                else None for i in range(len(atype))]
        atypes = ['C', 'pV']
        ahash = {'0 d':0, '1 d':1, '3 d':1, '7 d':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHoek2015 = getHoek2015


def getXing2014mm(self, tn=1):
    self.prepareData("COV225")
    atype = self.h.getSurvName("c treatment")
    atypes = ['PBS', 'IL15']
    ahash = {'IL15':1, 'PBS- ip':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXing2014mm = getXing2014mm


def getCribbs2018(self, tn=1):
    self.prepareData("COV226")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_.*", "", str(k)) for k in atype]
    atypes = ['DMSO', 'GSK-J4']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCribbs2018 = getCribbs2018


def getBayo2019(self, tn=1):
    self.prepareData("COV228")
    atype = self.h.getSurvName("c treatment")
    atypes = ['DMSO', 'GSK-J4', 'JIB-04', 'SD-70']
    ahash = {'24 h with DMSO':0,
            '24 h with JIB-04 (1mM)':2,
            '24 h with GSK-J4 (6.2 mM)':1,
            '24 h with SD-70 (2.2 mM)':3}
    if (tn == 2):
        atypes = ['DMSO', 'GSK-J4']
        ahash = {'24 h with DMSO':0, '24 h with GSK-J4 (6.2 mM)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBayo2019 = getBayo2019


def getNarang2018(self, tn=1):
    self.prepareData("COV229")
    atype = self.h.getSurvName("c gender")
    ahash = {'Female':0, 'Male':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c time point")
    ahash = {'Day 0':0, 'Day 2':2, 'Day 7':7, 'Day 28':28}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['Day 0', 'Day 2', 'Day 7', 'Day 28']
    ahash = {}
    if (tn == 2):
        atype = [sval[i] if tval[i] == 28
                else None for i in range(len(atype))]
        atypes = ['F', 'M']
        ahash = {0:0, 1:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNarang2018 = getNarang2018


def getHoang2014(self, tn=1, tb=2):
    self.prepareData("COV238")
    atype = self.h.getSurvName("c timepoint")
    ahash = {'Acute':0, 'Follow Up':1, 'day_28':1, 'day_0':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c severity")
    atypes = ['OFI', 'Mild', 'Severe']
    ahash = {'OFI':0, 'Mild':1, 'Moderate':1, 'Severe':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = tval
        atypes = ['CV', 'AV']
        ahash = {0:1, 1:0}
        atype = [atype[i] if aval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHoang2014 = getHoang2014


def getKim2015(self, tn=1, tb=2):
    self.prepareData("COV239")
    atype = self.h.getSurvName("c time point")
    ahash = {'48hrs':48, '2hrs':2, '12hrs':12, '8hrs':8, '36hrs':36,
            '4hrs':4, '72hrs':72, '60hrs':60, '24hrs':24, 'baseline':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c virus infection")
    atypes = ['U', 'IAV', 'RSV', 'IAV+RSV']
    ahash = {'Influenza & Rhino virus infected':3, 'Influenza virus infected':1,
            'Rhino virus infected':2, 'uninfected':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [0 if aval[i] == 0
                else tval[i] for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKim2015 = getKim2015


def getMuramoto2014(self, tn=1, tb=2):
    self.prepareData("COV240")
    atype = self.h.getSurvName("c time")
    ahash = {'Pre':0, '3':3, '1':1, '5':5, '7':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c influenza strain")
    atypes = ['0', '8', '9']
    ahash = {'VN30259':2, 'VN3040':0, 'VN3028II':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [tval[i] if aval[i] == tb
                else None for i in range(len(atype))]
        atypes = sorted(hu.uniq([i for i in atype if i is not None]))
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMuramoto2014 = getMuramoto2014


def getDavenport2015(self, tn=1, tb=2):
    self.prepareData("COV241")
    atype = self.h.getSurvName("c timepoint")
    ahash = {'Pre-challenge':0, '48 hours post-challenge':48,
            '12 hours post-challenge':12, '24 hours post-challenge':24}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c vaccination status")
    ahash = {'Control':0, 'Vaccinee':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c symptom severity")
    atypes = ['H', 'M', 'S']
    ahash = {'Moderate/severe':2, 'None':0, 'Mild':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 24
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDavenport2015 = getDavenport2015


def getChang2020(self, tn=1, tb=2):
    self.prepareData("COV242")
    atype = self.h.getSurvName("c clinical _w")
    ahash = {'Yes':1, 'No':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c virus_treatment_status")
    atypes = ['Control', 'RVA', 'RVC']
    ahash = {'Control':0, 'RVA':1, 'RVC':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [tval[i] if aval[i] == 2
                else None for i in range(len(atype))]
        atypes = ['No', 'Yes']
        ahash = {0:0, 1:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChang2020 = getChang2020


def getHuang2018(self, tn=1):
    self.prepareData("COV243")
    atype = self.h.getSurvName("c disease state")
    atypes = ['C', 'CP', 'AP']
    ahash = {'acute-phase KD':2, 'convalescent-phase KD':1, 'control':0}
    if (tn == 2):
        atypes = ['C', 'AP']
        ahash = {'acute-phase KD':1, 'control':0}
    if (tn == 3):
        atypes = ['CV', 'AV']
        ahash = {'acute-phase KD':1, 'convalescent-phase KD':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHuang2018 = getHuang2018


def getPopper2007I(self, tn=1):
    self.prepareData("COV244")
    atype = self.h.getSurvName("c Disease State")
    atypes = ['SA', 'TX', 'A', 'C', 'O']
    ahash = {'':4}
    if (tn == 2):
        atypes = ['Conv', 'SubAcute', 'Acute']
        ahash = {'C':0, 'SA':1, 'A':2}
    if (tn == 3):
        atypes = ['Conv', 'Acute']
        ahash = {'C':0, 'A':1}
    if (tn == 4):
        atypes = ['SubAcute', 'Acute']
        ahash = {'SA':0, 'A':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPopper2007I = getPopper2007I


def getPopper2007II(self, tn=1):
    self.prepareData("COV245")
    atype = self.h.getSurvName("c Phenotype")
    atypes = ['A', 'R', 'D', 'NR', 'N']
    ahash = {}
    if (tn == 2):
        atypes = ['R', 'NR']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPopper2007II = getPopper2007II


def getWright2018I(self, tn=1):
    self.prepareData("COV246")
    atype = self.h.getSurvName("c category")
    atypes = ['C', 'V', 'B', 'KD', 'U']
    ahash = {'Definite Viral':1, 'Control':0, 'Uncertain':4,
            'Definite Bacterial':2, 'Kawasaki Disease':3}
    if (tn == 2):
        atypes = ['C', 'KD']
        ahash = {'Control':0, 'Kawasaki Disease':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWright2018I = getWright2018I

def getWright2018Disc(self, tn=1):
    self.prepareData("COV246.2")
    atype = self.h.getSurvName("c category")
    atypes = ['C', 'V', 'B', 'KD', 'U']
    ahash = {'JIA':4, 'HSP':4, 'Uncertain':4,
            'Control':0, 'KD':3,
            'Definite Viral':1, 'Definite Bacterial':2}
    if (tn == 2):
        atypes = ['C', 'KD']
        ahash = {'Control':0, 'KD':1}
    if (tn == 3):
        atypes = ['V', 'B', 'KD']
        ahash = {'Definite Viral':0, 'Definite Bacterial':1, 'KD':2}
    if (tn == 4):
        atypes = ['B', 'KD']
        ahash = {'Definite Bacterial':0, 'KD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWright2018Disc = getWright2018Disc

def getWright2018II(self, tn=1):
    self.prepareData("COV247")
    atype = self.h.getSurvName("c category")
    atypes = ['C', 'V', 'B', 'KD', 'U', 'JIA', 'HSP']
    ahash = {'Definite Bacterial':2, 'Control':0, 'Uncertain':4,
            'Kawasaki Disease':3, 'Definite Viral':1}
    if (tn == 2):
        atypes = ['C', 'KD']
        ahash = {'Control':0, 'Kawasaki Disease':1}
    if (tn == 3):
        atypes = ['B', 'KD']
        ahash = {'Definite Bacterial':0, 'Kawasaki Disease':1}
    if (tn == 4):
        atypes = ['V', 'KD']
        ahash = {'Definite Viral':0, 'Kawasaki Disease':1}
    if (tn == 5):
        atypes = ['C', 'V']
        ahash = {'Control':0, 'Definite Viral':1}
    if (tn == 6):
        atypes = ['C', 'B']
        ahash = {'Control':0, 'Definite Bacterial':1}
    if (tn == 7):
        atypes = ['B', 'V']
        ahash = {'Definite Bacterial':0, 'Definite Viral':1}
    if (tn == 8):
        atypes = ['NS', 'S']
        atype = self.h.getSurvName('c Shock: 1=yes, 2=no')
        ahash = {'2.0':0, '1.0':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWright2018II = getWright2018II


def getPopper2009(self, tn=1):
    self.prepareData("COV248")
    atype = self.h.getSurvName("c disease state")
    atypes = ['C-dr', 'C-ai', 'C-sf', 'KD']
    ahash = {}
    if (tn == 2):
        atypes = ['C', 'KD']
        ahash = {'C-dr':0}
    if (tn == 3):
        atypes = ['C', 'KD']
        ahash = {'C-sf':0}
    if (tn == 4):
        atypes = ['C', 'KD']
        ahash = {'C-sf':0, 'C-dr':0, 'C-ai':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPopper2009 = getPopper2009


def getOkuzaki2017(self, tn=1):
    self.prepareData("COV249")
    atype = self.h.getSurvName("c ivig treatment")
    atypes = ['Pre', 'Post']
    ahash = {'before IVIg treatment':0, 'after IVIg treatment':1}
    if (tn == 2):
        atype = self.h.getSurvName("c gender")
        atypes = ['female', 'male']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOkuzaki2017 = getOkuzaki2017


def getJaggi2018(self, tn=1):
    self.prepareData("COV250")
    atype = self.h.getSurvName("c final condition")
    atypes = ['healthy', 'cKD', 'HAdV', 'inKD', 'GAS', 'GAS/SF']
    ahash = {'healthy':0, 'cKD':1, 'HAdV':2, 'inKD':3, 'GAS':4, 'GAS/SF':5}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['healthy', 'KD']
        ahash = {'inKD':1}
    if (tn == 3):
        atypes = ['healthy', 'cKD']
        ahash = {}
    if (tn == 4):
        atypes = ['healthy', 'inKD']
        ahash = {}
    if (tn == 5):
        atype = self.h.getSurvName("c race")
        atypes = ['W', 'B', 'A', 'H', 'O']
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJaggi2018 = getJaggi2018


def getHoang2014(self, tn=1):
    self.prepareData("COV251")
    atype = self.h.getSurvName("c ivig")
    ahash = {'Responsive':0, 'Resistant':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c phase")
    ahash = {'Convalescent':0, 'Acute':1}
    hval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c phenotype")
    atypes = ['OFI', 'Mild', 'Moderate', 'Severe']
    ahash = {'OFI':0, 'Mild':1, 'Moderate':2, 'Severe':3}
    pval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = hval
        atypes = ['CV', 'AV']
        ahash = {0:0, 1:1}
    if (tn == 3):
        atype = [atype[i] if hval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = tval
        atypes = ['R', 'NR']
        ahash = {0:0, 1:1}
        atype = [atype[i] if pval[i] is not None and pval[i] > 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHoang2014 = getHoang2014


def getOgata2009(self, tn=1):
    self.prepareData("COV252")
    atype = self.h.getSurvName("c time")
    ahash = {'post-treatment':1, 'pre-treatment':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c patient")
    ahash = {'IVIG-resistant patient':1, 'IVIG-responsive patient':0}
    pval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    ahash = {'IVIG':0, 'methylprednisolone (IVMP) + IVIG':1}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['IVIG', 'IVMP+IVIG']
    if (tn == 2):
        atype = pval
        atypes = ['R', 'NR']
        ahash = {0:0, 1:1}
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = pval
        atypes = ['R', 'NR']
        ahash = {0:0, 1:1}
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = tval
        atypes = ['Pre', 'Post']
        ahash = {0:0, 1:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOgata2009 = getOgata2009


def getFury2009(self, tn=1):
    self.prepareData("COV253")
    atype = self.h.getSurvName("c treatment_category")
    atypes = ['HC', 'R_A', 'R_IVIG', 'NR_A', 'NR_IVIG']
    ahash = {'IVIG responder_Acute':1, 'healthy controls':0,
            'IVIG non responder_Acute':3, 'IVIG non responder_after IVIG':4,
            'IVIG responder_after IVIG':2}
    if (tn == 2):
        atypes = ['HC', 'KD']
        ahash = {'healthy controls':0, 'IVIG responder_Acute':1,
                'IVIG non responder_Acute':1}
    if (tn == 3):
        atypes = ['Pre', 'Post']
        ahash = {'IVIG responder_Acute':0, 'IVIG responder_after IVIG':1,
                'IVIG non responder_Acute':0, 'IVIG non responder_after IVIG':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFury2009 = getFury2009


def getRowley2015(self, tn=1):
    self.prepareData("COV254")
    atype = self.h.getSurvName("c disease")
    atypes = ['C', 'tKD', 'uKD']
    ahash = {'untreated Kawasaki Disease':2, 'treated Kawasaki Disease':1, 'Control':0}
    if (tn == 2):
        atypes = ['HC', 'KD']
        ahash = {'Control':0, 'untreated Kawasaki Disease':1}
    if (tn == 3):
        atypes = ['Pre', 'Post']
        ahash = {'untreated Kawasaki Disease':0, 'treated Kawasaki Disease':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRowley2015 = getRowley2015


def getOCarroll2014(self, tn=1):
    self.prepareData("MACV182")
    atype = self.h.getSurvName("c cell state")
    atypes = ['C', 'A', 'T', 'R']
    ahash = {'Recovery':3, 'control':0, 'LPS Tolerance':2,
            'Acute response to LPS':1}
    if (tn == 2):
        atypes = ['C', 'A', 'T']
        ahash = {'control':0, 'LPS Tolerance':2, 'Acute response to LPS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOCarroll2014 = getOCarroll2014


def getMages2007(self, tn=1):
    self.prepareData("MACV183")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(".*with (.*) stim.*", "\\1", str(k)) for k in atype]
    atypes = ['1_0', '0_0', '1_1', '0_1']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMages2007 = getMages2007


def getFoster2007(self, tn=1):
    self.prepareData("MACV184")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(".*ages, ", "", str(k)) for k in atype]
    atypes = ['U', 'T', 'TR']
    ahash = {'treated with LPS for 24hours, then restimulated':2,
            'treated with LPS':1,
             'untreated':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFoster2007 = getFoster2007


def getUtz2020MmII(self, tn=1):
    self.prepareData("MACV185.2")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(".*- ", "", str(k)) for k in atype]
    atypes = ['Microglia', 'BAM']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getUtz2020MmII = getUtz2020MmII


def getHagemeyer2016Mm(self, tn=1):
    self.prepareData("MACV186")
    atype = self.h.getSurvName("c genotype")
    atypes = ['Wt', 'Irf8']
    ahash = {'Irf8 knockout':1, 'wildtype':0}
    if (tn == 2):
        atype = self.h.getSurvName("c tissue")
        atypes = ['liver', 'brain', 'skin', 'kidney', 'yolk sac']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHagemeyer2016Mm = getHagemeyer2016Mm


def getPoidinger2018Mm(self, tn=1):
    self.prepareData("MACV187")
    atype = self.h.getSurvName("c src1")
    atypes = ['M6h', 'M6l', 'S']
    ahash = {'microglia Ly6Chigh':0, 'microglia Ly6CLow':1, 'spleen monocytes':2}
    if (tn == 2):
        atypes = ['Ly6C+', 'Ly6C-']
        ahash = {'microglia Ly6Chigh':0, 'microglia Ly6CLow':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPoidinger2018Mm = getPoidinger2018Mm


def getAvraham2015Mm(self, tn=1, tb=0):
    self.prepareData("MACV188")
    atype = self.h.getSurvName("c time after salmonella exposure")
    ahash = {'0':0, '4':4, '2.5':2.5, '8':8}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c phrodo positive")
    ahash = {'Yes':1, 'No':0}
    pval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c gfp positive")
    atypes = ['GFP-', 'GFP+']
    ahash = {'Yes':1, 'No':0}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = [2 if gval[i] == 1
            else pval[i] for i in range(len(atype))]
    atypes = ['P-', 'P+', 'P+GFP+']
    ahash = {0:0, 1:1, 2:2}
    mval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = tval
        atypes = [0, 2.5, 4, 8]
        ahash = {}
    if (tn == 3):
        atype = tval
        atypes = [0, 2.5, 4, 8]
        ahash = {}
        atype = [atype[i] if mval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAvraham2015Mm = getAvraham2015Mm


def getDelorey2018Mm(self, tn=1):
    self.prepareData("MACV189")
    atype = self.h.getSurvName("c src1")
    atypes = ['M', 'M-', 'M+']
    ahash = {'macrophage':1, 'macrophage and C. albicans':2, '':0}
    if (tn == 2):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub("_.*", "", str(k)) for k in atype]
        atypes = ['Ct', 'U', 'M', 'CA', 'E']
        ahash = {'Ct':0, 'Uninfected':1, '':2, 'CandidaOnly':3, 'Exposed':4}
    if (tn == 3):
        atypes = ['M-', 'M+']
        ahash = {'macrophage':0, 'macrophage and C. albicans':1}
    if (tn == 4):
        atype = self.h.getSurvName("c Title")
        atype = [re.sub("_.*", "", str(k)) for k in atype]
        atypes = ['Ct', 'U']
        ahash = {'Ct':0, 'Uninfected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelorey2018Mm = getDelorey2018Mm


def getWestermann2016Hs(self, tn=1):
    self.prepareData("MACV190")
    atype = self.h.getSurvName("c infection agent")
    atypes = ['M', 'M+1', 'M+2']
    ahash = {'Salmonella typhimurium SL1344':2, 'Salmonella Typhimurium SL1344':2,
            'Salmonella Typhimurium SL1344 delta-pinT':1, '':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWestermann2016Hs = getWestermann2016Hs


def getWestermann2016Pig(self, tn=1):
    self.prepareData("MACV190.2")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(" h repl.*", "", str(k)) for k in atype]
    atype = [re.sub(".*pinT", "pinT", str(k)) for k in atype]
    atype = [re.sub(".*31 ", "", str(k)) for k in atype]
    atypes = ['pinT 00', 'pinT 01', 'pinT 02', 'pinT 06', 'pinT 16',
            'WT 00', 'WT 01', 'WT 02', 'WT 06', 'WT 16', 'mock 01']
    ahash = {}
    if (tn == 2):
        atypes = ['pinT 00', 'pinT 01', 'pinT 02', 'pinT 06', 'pinT 16']
    if (tn == 3):
        atypes = ['WT 00', 'WT 01', 'WT 02', 'WT 06', 'WT 16', 'mock 01']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWestermann2016Pig = getWestermann2016Pig


def getSanchezCabo2018Mm(self, tn=1):
    self.prepareData("MACV193")
    atype = self.h.getSurvName("c condition")
    atypes = ['Uh', '3', '7', '30']
    ahash = {'uninfarcted heart':0, '3 days post cryoinjury':1,
            '7 days post cryoinjury':2, '30 days post cryoinjury':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSanchezCabo2018Mm = getSanchezCabo2018Mm


def getSu2016tamMm(self, tn=1):
    self.prepareData("MACV194")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_\s*23.*", "", str(k)) for k in atype]
    atypes = ['mWT', 'mK', 'tWt', 'tK']
    ahash = {'KOM_CD206hi':1, 'WT_tumor':2, 'KOM_tumor':3, 'WT_CD206hi':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSu2016tamMm = getSu2016tamMm


def getPisu2020Mm(self, tn=1):
    self.prepareData("MACV195")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("M.*", "M", str(k)) for k in atype]
    atypes = ['AM', 'IM']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPisu2020Mm = getPisu2020Mm


def getChevalier2015Mm(self, tn=1):
    self.prepareData("MACV159")
    atype = self.h.getSurvName("c src1")
    atypes = ['rt', '6C', 'rt+anti', '6C+anti']
    ahash = {'room temp 30 days':0, '6 C 30 days':1,
            'room temp + antibiotics 30 days':2, '6 C + antibiotics 30 days':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChevalier2015Mm = getChevalier2015Mm


def getBlish2020PseudoBulk(self, tn=1):
    self.prepareData("COV255")
    atype = self.h.getSurvName("c sample origin")
    atypes = ['C', 'CoV']
    ahash = {'Patient with confirmed COVID-19':1, 'Healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlish2020PseudoBulk = getBlish2020PseudoBulk


def getBlish2020Mac(self, tn=1):
    self.prepareData("COV255.2")
    atype = self.h.getSurvName("c sample origin")
    atypes = ['C', 'CoV']
    ahash = {'Patient with confirmed COVID-19':1, 'Healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlish2020Mac = getBlish2020Mac


def getBlish2020B(self, tn=1):
    self.prepareData("COV255.3")
    atype = self.h.getSurvName("c sample origin")
    atypes = ['C', 'CoV']
    ahash = {'Patient with confirmed COVID-19':1, 'Healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlish2020B = getBlish2020B


def getBlish2020T(self, tn=1):
    self.prepareData("COV255.4")
    atype = self.h.getSurvName("c sample origin")
    atypes = ['C', 'CoV']
    ahash = {'Patient with confirmed COVID-19':1, 'Healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlish2020T = getBlish2020T


def getBlish2020NK(self, tn=1):
    self.prepareData("COV255.5")
    atype = self.h.getSurvName("c sample origin")
    atypes = ['C', 'CoV']
    ahash = {'Patient with confirmed COVID-19':1, 'Healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlish2020NK = getBlish2020NK


def getBlish2020Epi(self, tn=1):
    self.prepareData("COV255.6")
    atype = self.h.getSurvName("c sample origin")
    atypes = ['C', 'CoV']
    ahash = {'Patient with confirmed COVID-19':1, 'Healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlish2020Epi = getBlish2020Epi


def getPG2020Gut(self, tn=1):
    self.prepareData("COV256")
    atype = self.h.getSurvName("c Group")
    atypes = ['C', 'CoV']
    ahash = {'colon-48h':1, 'ileum-72h':1, 'ileum-un':0, 'colon-un':0, 'ileum-48h':1}
    if (tn == 2):
        ahash = {'colon-48h':1, 'colon-un':0}
    if (tn == 3):
        atypes = ['C', '48', '72']
        ahash = {'ileum-72h':2, 'ileum-un':0, 'ileum-48h':1}
    if (tn == 4):
        atypes = ['C', 'CoV']
        ahash = {'ileum-72h':1, 'ileum-un':0, 'ileum-48h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020Gut = getPG2020Gut


def getShalova2015(self, tn=1):
    self.prepareData("MACV196")
    atype = self.h.getSurvName("c treatment")
    ahash = {'None':0, 'LPS':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c src1")
    atypes = ['H', 'R', 'S']
    ahash = {'healthy donor':0,
            'patient recovering from sepsis':1,
            'patient with sepsis':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShalova2015 = getShalova2015


def getAGonzalez2017mm(self, tn=1, tb=0):
    self.prepareData("MACV197")
    atype = self.h.getSurvName("c tissue")
    ahash = {'Spleen':0, 'Bone marrow':1, 'Intestine':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c phenotype")
    atypes = ['NE', 'E']
    ahash = {'Non engulfing':0, 'Engulfing':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAGonzalez2017mm = getAGonzalez2017mm


def getButcher2018mm(self, tn=1, tb=0):
    self.prepareData("MACV200")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(",.*", "", str(k)) for k in atype]
    atypes = ['U', 'PA', 'PT', 'CA', 'CT', 'LA', 'LT', 'IA', 'IT']
    ahash = {'Untreated':0,
            'Pam3csk4 tolerance':2, 'Pam3csk4 acute':1, 'CpG acute':3, 'CpG tolerance':4,
            'LPS acute':5, 'LPS tolerance':6, 'Poly IC acute':7, 'Poly IC tolerance':8}
    if (tn == 2):
        atypes = ['U', 'LA', 'LT']
        ahash = {'Untreated':0, 'LPS acute':1, 'LPS tolerance':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getButcher2018mm = getButcher2018mm


def getBurns2020KD(self, tn=1, tb=0):
    self.prepareData("COV257")
    atype = self.h.getSurvName("c ethnicity")
    ahash = {'0': 'Unknown', '1': 'Asian', '2':'Black/ African American',
            '3': 'Caucasian', '4': 'Hispanic', '6': 'More than one race',
             '7': 'American Indian/Alaska Native',
            '8': 'Native Hawaiian or Other Pacific Islander', 9: 'Other'}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c disease_phase")
    atypes = ['C', 'A']
    ahash = {'Convalescent':0, 'Acute':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    idhash = {'UCSD-3815_C','UCSD-3450','UCSD-3929','UCSD-3291_A','UCSD-3318',
            'UCSD-3573','UCSD-3838','UCSD-3877_C','UCSD-3868_C','UCSD-3826_A',
            'UCSD-3836'}
    if (tn == 2):
        atype = [atype[i] if self.h.headers[i] not in idhash
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if self.h.headers[i] in idhash
                else None for i in range(len(atype))]
    if (tn == 4):
        btype = self.h.getSurvName("c with_matching_pair")
        atype = [atype[i] if btype[i] == 'yes'
                else None for i in range(len(atype))]
    if (tn == 5):
        atypes = ['C', 'A1', 'A2', 'A4']
        btype = self.h.getSurvName("c CA status")
        ahash = {'4':3, '2':2, '1':1}
        atype = [ahash[btype[i]] if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        ahash = {'Convalescent':0, 'Acute':1, 0:0, 1:1, 2:2, 3:3}
    if (tn == 6):
        atypes = ['W', 'B', 'A', 'H']
        atype = rval
        ahash = {'Caucasian':0, 'Black/ African American':1, 'Asian':2, 
                'Hispanic':3}
    if (tn == 7):
        atype = self.h.getSurvName("c illday")
        ahash = {'10':1}
        dval = [ahash[i] if i in ahash else None for i in atype]
        atypes = ['A2', 'A4']
        btype = self.h.getSurvName("c CA status")
        ahash = {'4':3, '2':2, '1':1}
        atype = [ahash[btype[i]] if aval[i] == 1 and dval[i] is None
                else atype[i] for i in range(len(atype))]
        ahash = {2:0, 3:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurns2020KD = getBurns2020KD


def getVallveJuanico2019(self, tn=1, tb=0):
    self.prepareData("MACV202")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_[0-9]*$", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'M1_Ctrl':1, 'M1_Endo':1, 'M2_Ctrl':2, 'M2_Endo':2}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1_Ctrl':1, 'M2_Ctrl':2}
    if (tn == 3):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1_Endo':1, 'M2_Endo':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVallveJuanico2019 = getVallveJuanico2019


def getCader2020Mm(self, tn=1, tb=0):
    self.prepareData("MACV203")
    atype = self.h.getSurvName("c phenotype")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'M0':0, 'M1 polarised':1, 'M2 polarised':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCader2020Mm = getCader2020Mm


def getWang2019(self, tn=1, tb=0):
    self.prepareData("MACV204")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub(", do.*", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {'IFNy/LPS':1, 'IL4':2, 'control media':0}
    if (tn == 2):
        ahash = {'UTD, IL13':2, 'UTD, IL4':2, 'UTD, media':0}
    if (tn == 3):
        ahash = {'CAR, IL13':2, 'CAR, IL4':2, 'CAR, media':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2019 = getWang2019


def getOConnell2019Mm(self, tn=1, tb=0):
    self.prepareData("MACV205")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2']
    ahash = {'RAW media':0, 'RAW M1':1, 'RAW M2':2,
            'BMM media':0, 'BMM M1':1, 'BMM M2':2}
    if (tn == 2):
        ahash = {'RAW media':0, 'RAW M1':1, 'RAW M2':2}
    if (tn == 3):
        ahash = {'BMM media':0, 'BMM M1':1, 'BMM M2':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOConnell2019Mm = getOConnell2019Mm


def getDiazJimenez2020I(self, tn=1, tb=0):
    self.prepareData("MACV206")
    atype = self.h.getSurvName("c cell type")
    atypes = ['Mono', 'Mac']
    ahash = {'Monocytes':0, 'macrophages-derived monocyte':1}
    if (tn == 2):
        atype = self.h.getSurvName("c agent")
        atypes = ['V', 'D']
        ahash = {'Vehicle':0, 'Dex':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDiazJimenez2020I = getDiazJimenez2020I


def getDiazJimenez2020II(self, tn=1, tb=0):
    self.prepareData("MACV206.2")
    atype = self.h.getSurvName("c cell subtype")
    atypes = ['Mono', 'Mac']
    ahash = {'Cell line':0, 'Cell line diffferentited':1}
    if (tn == 2):
        atype = self.h.getSurvName("c treatment")
        atypes = ['V', 'D']
        ahash = {'Vehicle':0, 'Dex treatment':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDiazJimenez2020II = getDiazJimenez2020II


def getGrossVered2020Mm(self, tn=1, tb=0):
    self.prepareData("MACV207")
    atype = self.h.getSurvName("c src1")
    ahash = {'Blood':0, 'Colon':1, 'Ileum':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c genotype")
    ahash = {'B6.FVB-Tg[Itgax-DTR/GFP]57Lan/J':0, 'Cx3cr1-DTR':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c selection marker")
    ahash = {'CD45+CD11b+CD115+Ly6C+/-':0,
            'DAPI-CD45+CD11b+CD64+Ly6C-MHCII+':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    ahash = {'N/A':0, '9-18 ng DT/gr bodyweight':1}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atypes = ['U', 'T']
    if (tn == 2):
        atype = tval
        atypes = ['B', 'C', 'I']
        ahash = {0:0, 1:1, 2:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGrossVered2020Mm = getGrossVered2020Mm


def getGopinathan2017(self, tn=1, tb=0):
    self.prepareData("MACV208")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_3h.*$", "", str(k)) for k in atype]
    atypes = ['U', 'Nm', 'IL10', 'Nm+Il10', 'm-Il10', 'p']
    ahash = {'Monocytes_unstimulated':0, 'Monocytes_IL-10':2, 'Monocytes_Nm':1,
            'Monocytes__mild_meningococcemia_plasma_immunodepleted for IL-10':4,
            'Monocytes_Nm+IL-10':3,
            'Monocytes__meningococcal_sepsis_plasma_immunodepleted_for_IL-10':4,
            'Monocytes__mild_meningococcemia_plasma':5,
            'Monocytes__meningitis_plasma_immunodepleted for IL-10':4,
            'Monocytes__meningococcal_sepsis_plasma':5,
            'Monocytes__meningitis_plasma':5}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'Monocytes_unstimulated':0, 'Monocytes_IL-10':2, 'Monocytes_Nm':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGopinathan2017 = getGopinathan2017


def getSerazin2020(self, tn=1, tb=0):
    self.prepareData("MACV209")
    atype = self.h.getSurvName("c cell type")
    atypes = ['Mono', 'Mac']
    ahash = {'monocytes':0, 'macrophages':1}
    if (tn == 2):
        atype = self.h.getSurvName("c src1")
        atypes = ['M0', 'M1', 'M2', 'M-IL34', 'M-MCSF']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSerazin2020 = getSerazin2020


def getBlack2018(self, tn=1, tb=0):
    self.prepareData("MACV210")
    atype = self.h.getSurvName("c cell type")
    atypes = ['M', 'Mc', 'Mnc']
    ahash = {'Monocyte':0,
            'MonocyteNonclassical':2,
            'MonocyteClassical':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlack2018 = getBlack2018


def getMuhitch2019(self, tn=1, tb=0):
    self.prepareData("MACV211")
    atype = self.h.getSurvName("c src1")
    atypes = ['Mc', 'Mnc', 'DCc', 'DCnc']
    ahash = {'non-classical monocyte derived dendritic cell':3,
            'classical monocyte derived dendritic cell':2,
            'non-classical monocytes':1,
            'classical monocytes':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMuhitch2019 = getMuhitch2019


def getMuhitch2019MonoI(self, tn=1, tb=0):
    self.prepareData("MACV212")
    atype = self.h.getSurvName("c disease state")
    atypes = ['C', 'RCC']
    ahash = {'control':0, 'renal cell carcinoma (RCC)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMuhitch2019MonoI = getMuhitch2019MonoI


def getMuhitch2019MonoII(self, tn=1, tb=0):
    self.prepareData("MACV213")
    atype = self.h.getSurvName("c monocyte subset")
    ahash = {'intermediate':1, 'classical':0, 'non-classical':2}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c disease state")
    atypes = ['C', 'RCC']
    ahash = {'renal cell carcinoma':1, 'healthy':0}
    if (tn == 2):
        atype = [atype[i] if mval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMuhitch2019MonoII = getMuhitch2019MonoII


def getXu2019(self, tn=1, tb=0):
    self.prepareData("MACV215")
    atype = self.h.getSurvName("c cell type")
    atypes = ['Mc', 'Mi', 'Mnc']
    ahash = {'Classical monocytes':0,
            'Intermediate monocytes':1,
            'Non classical monocytes':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXu2019 = getXu2019


def getBouchlaka2017(self, tn=1, tb=0):
    self.prepareData("MACV216")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(" Mac.*", "", str(k)) for k in atype]
    atypes = ['BM', 'BL', 'MSC']
    ahash = {'Blood Derived':1, 'MSC Educated':2, 'BM Derived':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBouchlaka2017 = getBouchlaka2017


def getTausendschon2015(self, tn=1, tb=0):
    self.prepareData("MACV217")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub("hu.* \(", "", str(k)) for k in atype]
    atype = [re.sub("\) .*IL-10", " IL10", str(k)) for k in atype]
    atype = [re.sub("\)\s*16.*, ", " -IL10 ", str(k)) for k in atype]
    atypes = ['C1', 'C21', 'IL10-1', 'IL10-21']
    ahash = {'si control IL10, 4h 1% oxygen':2,
            'si control IL10, 4h 21% oxygen':3,
            'si control -IL10 4h 1% oxygen':0,
            'si control -IL10 4h 21% oxygen':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTausendschon2015 = getTausendschon2015


def getWhite2012(self, tn=1, tb=0):
    self.prepareData("MACV218")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'LPS', '4F', '4F+LPS']
    ahash = {'control':0, '4F':2, '4F + lipopolysaccharides (LPS)':3,
            'control + lipopolysaccharides (LPS)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWhite2012 = getWhite2012


def getJardine2019(self, tn=1, tb=0):
    self.prepareData("MACV219")
    atype = self.h.getSurvName("c cell type")
    tissue = [re.sub(",.*", "", str(k)) for k in atype]
    ahash = {'Blood':0, 'BAL':1}
    tval = [ahash[i] if i in ahash else None for i in tissue]
    ctype = [re.sub(".*, ", "", str(k)) for k in atype]
    ahash = {'Classical Monocyte':0, 'int. Monocyte':1,
            'Non-classical Monocyte':2, 'monocyte-derived DC':3,
            'AM':4, 'DC1':5, 'DC2':6, 'DC3':7, 'DC2/3':8, 'PDC':9}
    sval = [ahash[i] if i in ahash else None for i in ctype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'LPS']
    ahash = {'Healthy':0, 'Saline':0, 'LPS':1}
    if (tn == 2):
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [sval[i] if tval[i] == 0
                else None for i in range(len(atype))]
        atypes = ['C', 'I', 'NC']
        ahash = {0:0, 1:1, 2:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJardine2019 = getJardine2019


def getSansom2020(self, tn=1, tb=0):
    self.prepareData("MACV221")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_hh.*$", "", str(k)) for k in atype]
    atype = [re.sub("cLP_", "", str(k)) for k in atype]
    atypes = ['mac_ko', 'mac_wt', 'p1mono_ko', 'p1mono_wt', 'p2mono_ko', 'p2mono_wt']
    ahash = {}
    if (tn == 2):
        atypes = ['mac_wt', 'p1mono_wt', 'p2mono_wt']
    if (tn == 3):
        atypes = ['mac_wt', 'mac_ko']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSansom2020 = getSansom2020


def getLi2019Mm(self, tn=1, tb=0):
    self.prepareData("MACV222")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_.*$", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2', 'lnATM', 'obATM']
    ahash = {}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
    if (tn == 3):
        atypes = ['lnATM', 'obATM']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2019Mm = getLi2019Mm


def getLi2019MmBlk(self, tn=1, tb=0):
    self.prepareData("MACV222.2")
    atype = self.h.getSurvName("c Title")
    atypes = ['M0', 'M1', 'M2', 'P']
    ahash = {'M012':3,
            'A0.0':1, 'A0.1':1, 'A0.2':1, 'A0.3':1, 'A0.4':1,
            'A0.5':2, 'A0.6':2, 'A0.7':2, 'A0.8':2, 'A0.9':2, 'A1.0':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2019MmBlk = getLi2019MmBlk


def getRealegeno2016(self, tn=1, tb=0):
    self.prepareData("MACV223")
    atype = self.h.getSurvName("c stimulation")
    atypes = ['C', 'IFNG', 'TLR2/1']
    ahash = {'media alone':0, 'TLR2/1 ligand':2, 'interferon gamma':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRealegeno2016 = getRealegeno2016


def getOkuzaki2020(self, tn=1, tb=0):
    self.prepareData("COV258")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("_.*", "", str(k)) for k in atype]
    atypes = ['A549', 'lungNHBE', 'lungORG', 'V', 'VC', 'Ctl']
    ahash = {}
    if (tn == 2):
        atypes = ['C', 'CoV']
        ahash = {'V':1, 'Ctl':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOkuzaki2020 = getOkuzaki2020


def getBrykczynska2020(self, tn=1, tb=0):
    self.prepareData("MACV224")
    atype = self.h.getSurvName("c tissue")
    ahash = {'adipose tissue macrophages':0, 'monocytes':1, 'microglia':2,
            'colonic macrophages':3, 'islet macrophages':4, 'Kupffer cells':5,
            'peritoneal macrophages':6}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['N', 'HFD', 'STZ', 'F', 'RE']
    ahash = {'fasted':3, 'refed':4, 'normal_diet':0, 'high_fat_diet_with_STZ':2,
            'high_fat_diet':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBrykczynska2020 = getBrykczynska2020


def getLavin2014(self, tn=1, tb=0):
    self.prepareData("MACV225")
    atype = self.h.getSurvName("c cell type")
    atypes = ['Nu', 'Mono', 'B', 'PC', 'Sp', 'Li', 'Co', 'Il', 'Lu']
    ahash = {'Neutrophils':0, 'Monocytes':1, 'Brain microglia':2,
            'Peritoneal cavity macrophages':3, 'Spleen red pulp macrophages':4,
            'Kupffer cells':5, 'Large intestine macrophages':6,
            'Small intestine macrophages':7, 'Lung macrophages':8}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLavin2014 = getLavin2014


def getHaney2018(self, tn=1, tb=0):
    self.prepareData("MACV226")
    atype = self.h.getSurvName("c protocol")
    atypes = ['U', 'D', 'Uc', 'Dc']
    ahash = {'Undifferentiated GFP':0, 'Undifferentiated Ctrl sgRNA':2,
            'Differentiated GFP':1, 'Differentiated Ctrl sgRNA':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHaney2018 = getHaney2018


def getDutertre2019(self, tn=1, tb=0):
    self.prepareData("MACV227.2")
    atype = self.h.getSurvName("c disease state")
    atypes = ['H', 'SLE', 'Ss', 'NA']
    ahash = {'Healthy':0, 'NA':3, 'SLE':1, 'Systemic sclerosis':2}
    if (tn == 2):
        atype = self.h.getSurvName("c cell type")
        atypes = ['5+', '5-163-', '163+14-', '163+14+']
        ahash = {'cDC2_CD5+':0, 'cDC2_CD5-CD163-':1,
                'cDC2_CD5-CD163+CD14-':2, 'cDC2_CD5-CD163+CD14+':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDutertre2019 = getDutertre2019


def getBeins2016I(self, tn=1, tb=0):
    self.prepareData("MACV228")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2', 'TGFb']
    ahash = {'TGFb':3, 'none':0, 'LPS/IFNg':1, 'IL-4':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBeins2016I = getBeins2016I


def getBeins2016II(self, tn=1, tb=0):
    self.prepareData("MACV228.2")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1', 'M2', 'TGFb']
    ahash = {'TGFb':3, 'none':0, 'LPS/IFNg':1, 'IL-4':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBeins2016II = getBeins2016II


def getBeins2016III(self, tn=1, tb=0):
    self.prepareData("MACV228.3")
    atype = self.h.getSurvName("c stimulation")
    atypes = ['M0', 'M1', 'M2', 'TGFb']
    ahash = {'unstimulated':0, 'LPS+IFNg':1, 'IL4':2, 'TGFb':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBeins2016III = getBeins2016III


def getRossi2018Mm(self, tn=1, tb=0):
    self.prepareData("MACV229")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['GFP', 'IL4']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRossi2018Mm = getRossi2018Mm


def getWel2020(self, tn=1, tb=0):
    self.prepareData("MACV230")
    atype = self.h.getSurvName("c genotype/variation")
    atypes = ['WT', 'FES']
    ahash = {'Wild-type':0, 'FES_S700C mutant':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWel2020 = getWel2020


def getSingh2016(self, tn=1, tb=0):
    self.prepareData("MACV231")
    atype = self.h.getSurvName("c src1")
    atypes = ['CD15-', 'CD15+']
    ahash = {'CD15-':0, 'CD15+ cancer stem like cells':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSingh2016 = getSingh2016


def getJoshi2014(self, tn=1, tb=0):
    self.prepareData("MACV232")
    atype = self.h.getSurvName("c genotype/variation")
    atypes = ['WT', 'Rac2']
    ahash = {'Rac2 -/-':1, 'WT mice':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJoshi2014 = getJoshi2014


def getJoshi2020MmI(self, tn=1, tb=0):
    self.prepareData("MACV233")
    atype = self.h.getSurvName("c src1")
    atypes = ['Un', 'LPS', 'IL4', 'TAM', 'TAMSyk']
    ahash = {'BMDMs, IL4-stimulated, Syk flox':2,
            'BMDMs, IL4-stimulated, Syk cre':2,
            'BMDMs, not stimulated, Syk cre':0,
            'BMDMs, LPS-stimulated, Syk flox':1,
            'BMDMs, LPS-stimulated, Syk cre':1,
            'BMDMs, not stimulated, Syk flox':0,
            'Tumor-associated macrophages, not stimulated, Syk flox':4,
            'Tumor-associated macrophages, not stimulated, Syk cre':3}
    if (tn == 2):
        atypes = ['C', 'T']
        ahash = {'Tumor cells, vehicle treated, Syk WT':0,
                'Tumor cells, SRX3207 treated, Syk WT':1,
                'Tumor cells, not stimulated, Syk cre':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJoshi2020MmI = getJoshi2020MmI


def getJoshi2020MmII(self, tn=1, tb=0):
    self.prepareData("MACV233.2")
    atype = self.h.getSurvName("c Group")
    atypes = ['WT', 'Syk']
    ahash = {'KO':1, 'WT':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJoshi2020MmII = getJoshi2020MmII


def getJoshi2020MmIII(self, tn=1, tb=0):
    self.prepareData("MACV233.3")
    atype = self.h.getSurvName("c Title")
    atypes = ['WT', 'Syk']
    ahash = {'CD11b_Syk_WT':0, 'CD11b_Syk_KO':1, 'CD45_Syk_KO':1,
            'CD45_Syk_WT':0, 'Tumor_Syk_KO':1, 'Tumor_Syk_WT':0,
            'CD11b_Syk_WT_mac':0, 'CD11b_Syk_KO_mac':1, 'CD45_Syk_KO_mac':1,
            'CD45_Syk_WT_mac':0, 'Tumor_Syk_KO_mac':1, 'Tumor_Syk_WT_mac':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJoshi2020MmIII = getJoshi2020MmIII


def getRedmond2020(self, tn=1, tb=0):
    self.prepareData("COV259")
    atype = self.h.getSurvName("c tissue")
    ahash = {'hESC Pancreas':0, 'Lung':1, 'Liver Organoid':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infected")
    atypes = ['C', 'CoV']
    ahash = {'Mock':0, 'sars-Cov2':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c Title')
        atype = [re.sub(".*SARS", "SARS", str(k)) for k in atype]
        atype = [re.sub(".*Mock", "Mock", str(k)) for k in atype]
        atypes = ['C', 'CoV']
        ahash = {'Mock 1':0, 'Mock 2':0, 'SARS-CoV-2 1':1, 'SARS-CoV-2 2':1}
    if (tn == 4):
        atype = self.h.getSurvName('c Title')
        atype = [re.sub(".*SARS", "SARS", str(k)) for k in atype]
        atype = [re.sub(".*Mock", "Mock", str(k)) for k in atype]
        atypes = ['C', 'CoV']
        ahash = {'Mock 1 R2':0, 'Mock 2 R2':0, 'SARS-CoV-2 1 R2':1, 'SARS-CoV-2 2 R2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRedmond2020 = getRedmond2020


def getHe2020Mm(self, tn=1, tb=0):
    self.prepareData("COV260")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'hACE2']
    ahash = {'control':0, 'Ad5-hACE2-transduced':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHe2020Mm = getHe2020Mm


def getWang2020Cov(self, tn=1, tb=0):
    self.prepareData("COV262")
    atype = self.h.getSurvName("c infection")
    atypes = ['M', 'CoV']
    ahash = {'SARS-CoV-2':1, 'mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2020Cov = getWang2020Cov


def getPG2020HSAE(self, tn=1, tb=0):
    self.prepareData("COV263")
    atype = self.h.getSurvName("c Group")
    atypes = ['M', 'CoV']
    ahash = {'CoV-72':1, 'U':0, 'CoV-48':0}
    if (tn == 2):
        atype = self.h.getSurvName("c Protocol")
        atypes = ['M', 'CoV', '96well']
        ahash = {'72':1, '0':0, '96 well 48':2, '48':1, '96 well 72':2}
    if (tn == 3):
        atype = self.h.getSurvName("c Protocol")
        atypes = ['M', 'CoV']
        ahash = {'72':1, '0':0, '96 well 48':1, '48':1, '96 well 72':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020HSAE = getPG2020HSAE


def getJulia2020(self, tn=1, tb=0):
    self.prepareData("COV264")
    atype = self.h.getSurvName("c treatment")
    atypes = ['control', 'abatacept']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJulia2020 = getJulia2020


def getRuppin2020(self, tn=1, tb=0):
    self.prepareData("COV265")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'I']
    ahash = {'Control (mock infection)':0, 'SARS-CoV-2 infection':1}
    if (tn == 2):
        atype = self.h.getSurvName("c src1")
        ahash = {'Vero E6 cells':0, 'SARS-CoV-2-infected Vero E6 cells':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRuppin2020 = getRuppin2020


def getChua2020Blk(self, tn=1, tb=0):
    self.prepareData("COV267")
    atype = self.h.getSurvName("c Type")
    ahash = {'NS':0, 'BL':1, 'PB':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c COVID-19 severity')
    atypes = ['C', 'M', 'S']
    ahash = {'':0, 'critical':2, 'moderate':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'CoV']
        ahash = {'':0, 'critical':1, 'moderate':1}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChua2020Blk = getChua2020Blk


def getPG2020LungHs(self, tn=1, tb=0):
    #self.prepareData("covirus1",
    #        cfile="/Users/sataheri/public_html/Hegemon/explore.conf")
    self.prepareData("COV372")
    atype = self.h.getSurvName('c title')
    atype = [str(k).split("_")[2] if len(str(k).split("_")) > 2
                     else None for k in atype]
    ahash = {'48h':48, 'UN':0, '72h':72, 'Un':0}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c type")
    ahash = {'monolayer':0, 'ALI':1, 'organoid_P1':2, 'organoid_P4':2,
            'Tissue_a':3, 'Tissue_b':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    ahash = {'infected':1, 'Uninfected':0}
    atypes = ['Un', 'I']
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['U', '48', '72']
        atype = sval
        ahash = {0:0, 48:1, 72:2}
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020LungHs = getPG2020LungHs


def getPG2020LungHam(self, tn=1, tb=0):
    #self.prepareData("covirus1.2",
    #        cfile="/Users/sataheri/public_html/Hegemon/explore.conf")
    self.prepareData("COV373")
    atype = self.h.getSurvName("c title")
    atype = [re.sub("[-_].*", "", str(k)) for k in atype]
    ahash = {}
    atypes = ['UN', '3', '4']
    if (tn == 2):
        atypes = ['UN', '3']
    if (tn == 3):
        atypes = ['3', '4']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020LungHam = getPG2020LungHam


def getHan2020(self, tn=1, tb=0):
    self.prepareData("COV269")
    atype = self.h.getSurvName("c tissue/cell type")
    ahash = {'hPSC_Lung organoid':0, 'adult lung autopsy':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c subject status')
    atypes = ['H', 'CoV']
    ahash = {'healthy':0, 'COVID-19':1}
    if (tn == 2):
        atype = self.h.getSurvName("c treatment")
        atypes = ['M', 'CoV', 'DMSO', 'Im']
        ahash = {'Mock':0, 'SARS-CoV-2':1, 'SARS-CoV-2+imatinib':3,
                'SARS-CoV-2+DMSO':2}
    if (tn == 3):
        atype = self.h.getSurvName("c treatment")
        atypes = ['M', 'CoV']
        ahash = {'Mock':0, 'SARS-CoV-2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHan2020 = getHan2020


def getMaayan2020(self, tn=1, tb=0):
    self.prepareData("COV270")
    atype = self.h.getSurvName("c src1")
    ahash = {'A549-ACE2':0, 'Pancreatic organoids':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c drug")
    ahash = {'mock':0, 'DMSO':1, 'Amlodipine':2,
            'Terfenadine':2, 'Loperamide':2, 'Berbamine':2,
            'Trifluoperazine':2, 'RS504395':2, 'RS504393':2, 'RS504394':2}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c sars-cov-2 infected')
    atypes = ['H', 'CoV']
    ahash = {'Yes':1, 'No':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['H', 'CoV', 'T']
        atype = [dval[i] if tval[i] == 1
                else None for i in range(len(atype))]
        ahash = {0:0, 1:1, 2:2}
    if (tn == 5):
        atype = [atype[i] if tval[i] == 0 and (dval[i] == 0 or dval[i] == 1)
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaayan2020 = getMaayan2020


def getChen2020(self, tn=1, tb=0):
    self.prepareData("COV271")
    atype = self.h.getSurvName('c Title')
    atypes = ['H', 'CoV']
    ahash = {'GSM4451223':0, 'GSM4451224':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChen2020 = getChen2020


def getLeem2020Mm(self, tn=1, tb=0):
    self.prepareData("COV272")
    atype = self.h.getSurvName('c src1')
    atypes = ['N', 'Il15', 'NI', 'NIL15', 'TC', 'IFN']
    ahash = {'Naive mouse without IL-15':0,
            'Naive mouse with IL-15':1,
            'MCMV-infected mouse without IL-15':2,
            'MCMV-infected mouse with IL-15':3,
            'TC_1_Ct':4,
            'TC_1_IFNr':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLeem2020Mm = getLeem2020Mm


def getLangelier2020(self, tn=1, tb=0):
    self.prepareData("COV273")
    atype = self.h.getSurvName('c disease state')
    atypes = ['C', 'CoV', 'V']
    ahash = {'no virus':0, 'SC2':1, 'other virus':2}
    if (tn == 2):
        atypes = ['C', 'CoV']
        ahash = {'no virus':0, 'SC2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLangelier2020 = getLangelier2020


def getJaitovich2020(self, tn=1, tb=0):
    self.prepareData("COV274")
    hfd = self.h.getSurvName('c hospital-free days post 45 day followup (days)')
    atype = self.h.getSurvName('c icu')
    ahash = {'yes':1, 'no':0}
    ival = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c mechanical ventilation')
    ahash = {'yes':1, 'no':0}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state')
    atypes = ['C', 'CoV']
    ahash = {'COVID-19':1, 'non-COVID-19':0}
    dval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = mval
        atypes = ['C', 'MV']
        ahash = {1:1, 0:0}
    if (tn == 3):
        atype = mval
        atypes = ['C', 'ICU']
        ahash = {1:1, 0:0}
    if (tn == 4):
        st = ["-".join([str(ival[i]), str(mval[i])])
                for i in range(len(atype))]
        atypes = ['0-0', '0-1', '1-0', '1-1']
        atypes = ['1-1']
        ahash = {}
        atype = hfd
        atype[1] = 0
        atype = [0 if int(k) >= 20 else 1 for k in atype]
        atype = [atype[i] if st[i] == '1-1' and dval[i] == 1
                and hfd[i] != '0'
                else None for i in range(len(atype))]
        ahash = {0:0, 1:1}
        atypes = ['G', 'P']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJaitovich2020 = getJaitovich2020


def getGeorge2019Mm(self, tn=1, tb=0):
    self.prepareData("MACV234")
    atype = self.h.getSurvName('c injected with')
    atypes = ['PBS', 'LPS', 'IL4']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGeorge2019Mm = getGeorge2019Mm


def getNugent2020Mm(self, tn=1, tb=0):
    self.prepareData("MACV235")
    atype = self.h.getSurvName('c cell_type')
    ahash = {'microglia':0, 'other':2, 'astrocyte':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c genotype')
    atypes = ['WT', 'Het', 'M']
    ahash = {'Trem2 +/-':1, 'Trem2 -/-':2, 'Trem2 +/+':0,
            'TREM2 +/+':0, 'TREM2 -/-':2}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = [atype[i] if tval[i] == 2
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = tval
        atypes = ['M', 'A', 'O']
        ahash = {0:0, 1:1, 2:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNugent2020Mm = getNugent2020Mm


def getGutbier2020(self, tn=1, tb=0):
    self.prepareData("MACV236")
    atype = self.h.getSurvName('c cell type')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'iPSC-derived M0 macrophage':0,
            'iPSC-derived M1 macrophage':1,
            'iPSC-derived M2 macrophage':2,
            'iPSC-derived macrophage progenitors':0,
            'PBMC monocyte':0,
            'iPSC-derived microglia cell':0,
            'PBMC-derived M0 macrophage':0,
            'PBMC-derived M1 macrophage':1,
            'PBMC-derived M2 macrophage':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGutbier2020 = getGutbier2020


def getDuan2020(self, tn=1, tb=0):
    self.prepareData("MACV237")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'M1', 'M2', 'V']
    ahash = {'GSM4557108':0, 'GSM4557109':1, 'GSM4557110':3,
            'GSM4557111':2, 'GSM4557112':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuan2020 = getDuan2020


def getDuan2020Mac(self, tn=1, tb=0):
    self.prepareData("MACV237.2")
    atype = self.h.getSurvName('c Title')
    atypes = ['V', 'M1', 'M2']
    ahash = {'GSM4557108':0, 'GSM4557109':1, 'GSM4557110':0,
            'GSM4557111':2, 'GSM4557112':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuan2020Mac = getDuan2020Mac


def getDuan2020MacII(self, tn=1, tb=0):
    self.prepareData("MACV237.3")
    atype = self.h.getSurvName('c Title')
    atypes = ['V', 'M1', 'M2']
    ahash = {'M0':0, 'A0.0':1, 'A0.1':1, 'A0.2':1, 'A0.3':1, 'A0.4':1,
            'A0.5':2, 'A0.6':2, 'A0.7':2, 'A0.8':2, 'A0.9':2, 'A1.0':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuan2020MacII = getDuan2020MacII


def getDuan2020MacIII(self, tn=1, tb=0):
    self.prepareData("MACV237.4")
    atype = self.h.getSurvName('c Title')
    atypes = ['M0', 'M1', 'M2']
    ahash = {'M012':0,
            'A0.0':1, 'A0.1':1, 'A0.2':1, 'A0.3':1, 'A0.4':1,
            'A0.5':2, 'A0.6':2, 'A0.7':2, 'A0.8':2, 'A0.9':2, 'A1.0':2,
            'B0.0':1, 'B0.1':1, 'B0.2':1, 'B0.3':1, 'B0.4':1,
            'B0.5':2, 'B0.6':2, 'B0.7':2, 'B0.8':2, 'B0.9':2, 'B1.0':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuan2020MacIII = getDuan2020MacIII


def getDuan2020Lung(self, tn=1, tb=0):
    self.prepareData("MACV237.5")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'M1', 'M2', 'V']
    ahash = {'GSM4557108':0, 'GSM4557109':1, 'GSM4557110':3,
            'GSM4557111':2, 'GSM4557112':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuan2020Lung = getDuan2020Lung


def getDuan2020LungII(self, tn=1, tb=0):
    self.prepareData("MACV237.6")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'CoV']
    ahash = {'L0':0,
            'A0.0':0, 'A0.1':0, 'A0.2':0, 'A0.3':0, 'A0.4':0,
            'A0.5':1, 'A0.6':1, 'A0.7':1, 'A0.8':1, 'A0.9':1, 'A1.0':1,
            'B0.0':0, 'B0.1':0, 'B0.2':0, 'B0.3':0, 'B0.4':0,
            'B0.5':1, 'B0.6':1, 'B0.7':1, 'B0.8':1, 'B0.9':1, 'B1.0':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDuan2020LungII = getDuan2020LungII


def getGurvich2020(self, tn=1, tb=0):
    self.prepareData("MACV238")
    atype = self.h.getSurvName('c cell type')
    atypes = ['M0', 'M1', 'M2', 'O']
    ahash = {'CD14':3, 'M0':0, 'M1':1, 'M2a':2, 'Mreg':3, 'Mreg_UKR':3, 'PCMO':3}
    if (tn == 2):
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M0':0, 'M1':1, 'M2a':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGurvich2020 = getGurvich2020


def getNienhold2020(self, tn=1, tb=0):
    self.prepareData("COV276")
    atype = self.h.getSurvName('c diagnosis')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['C', 'CoV', 'O']
    ahash = {'COVID-19':1, 'Control':0, 'Other':2}
    if (tn == 2):
        atypes = ['C', 'CoV']
        ahash = {'COVID-19':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNienhold2020 = getNienhold2020


def getSuarez2015(self, tn=1):
    self.prepareData("MACV133")
    atype = self.h.getSurvName("c condition")
    atypes = ['C', 'S']
    ahash = {'COINFECTION':1, 'Healthy Control':0, 'VIRUS':1, 'BACTERIA':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'B', 'V', 'CO']
        ahash = {'COINFECTION':3, 'Healthy Control':0, 'VIRUS':2, 'BACTERIA':1}
    if (tn == 3):
        atype = self.h.getSurvName("c race")
        atypes = ['W', 'B', 'A']
        ahash = {'White':0, 'Black':1, 'Asian':2}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSuarez2015 = getSuarez2015


def getAhn2015(self, tn=1):
    self.prepareData("MACV142")
    atype = self.h.getSurvName("c pathogen")
    atypes = ['C', 'S']
    ahash = {'-':0,
            'Staphylococcus aureus':1,
            'Escherichia coli':1,
            'Staphylococcus aureus and Streptococcus pneumoniae':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'Ec', 'Sa', 'S']
        ahash = {'-':0,
                'Staphylococcus aureus':2,
                'Escherichia coli':1,
                'Staphylococcus aureus and Streptococcus pneumoniae':3}
    if (tn == 3):
        atype = self.h.getSurvName("c ethnicity")
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'White':0, 'Unknown':3, 'Black':1, 'unknown':3, 'Asian':2,
                'black':1, 'white':0}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAhn2015 = getAhn2015


def getBloom2013(self, tn=1):
    self.prepareData("MACV155")
    atype = self.h.getSurvName("c disease state")
    atypes = ['C', 'TB', 'B', 'P', 'S', 'C']
    ahash = {'TB':1, 'Sarcoid':4, 'Control':0, 'pneumonia':3,
            'Active sarcoidosis':4, 'Non-active sarcoidosis':4,
            'lung cancer':5, 'Active Sarcoid':4, 'Pneumonia':3,
            'Baseline':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'TB']
        ahash = {'TB':1, 'Control':0}
    if (tn == 3):
        atype = self.h.getSurvName("c ethnicity")
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'Afro-Carribbean':3, 'Indian subcontinent':2,
                '':3, 'Black':1, 'White':0, 'Caucasian':0,
                'Middle Eastern':3, 'SE Asian':2, 'Central Asia':2,
                'None':3}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBloom2013 = getBloom2013


def getYang2017(self, tn=1):
    self.prepareData("MACV252")
    atype = self.h.getSurvName("c asthma")
    atypes = ['C', 'A']
    ahash = {'TRUE':1, 'FALSE':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 3):
        atype1 = self.h.getSurvName('c race_white')
        atype2 = self.h.getSurvName('c race_aa')
        atype3 = self.h.getSurvName('c race_hispanic')
        atype = [" ".join([str(atype1[i]), str(atype2[i]), str(atype3[i])]) \
                for i in range(len(atype1))]
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'No Yes No':1, 'Yes Yes Yes':3, 'No No Yes':3,
           'Yes No Yes':0, 'Yes Yes No':3, 'No No No':3, 'No Yes Yes':1}
        #atype = [atype[i] if aval[i] == 0
        #        else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2017 = getYang2017


def getLewis2019(self, tn=1):
    self.prepareData("COV159")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'V']
    ahash = {'VARILRIX':1, 'ENGERIXB3':1, 'PLACEBOB3':0, 'ENGERIXB1':1,
             'AGRIPPAL':1, 'PLACEBOAB1C':0, 'STAMARIL':1, 'FLUADC':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['C', 'VAR', 'ENG', 'AGR', 'STA', 'FLU']
        ahash = {'VARILRIX':1, 'ENGERIXB3':2, 'PLACEBOB3':0, 'ENGERIXB1':2,
                 'AGRIPPAL':3, 'PLACEBOAB1C':0, 'STAMARIL':4, 'FLUADC':5}
    if (tn == 3):
        atype = self.h.getSurvName('c race')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'white':0, 'asian':2, 'black or african american':1, 'other':3}
        #atype = [atype[i] if aval[i] == 0
        #        else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLewis2019 = getLewis2019


def getObermoser2013(self, tn=1):
    self.prepareData("COV160")
    atype = self.h.getSurvName("c vaccine")
    atypes = ['C', 'I', 'V', 'NS']
    ahash = {'saline':0, 'PNEUM':1, 'Pneumovax':2, 'FLUZONE':2, 'Saline':0,
            'Pneumovax vaccine group':2, 'Influenza vaccine group':2,
            'Flu':1, 'NS':3, 'Influenza':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['C', 'I', 'V']
        ahash = {'saline':0, 'PNEUM':1, 'Pneumovax':2, 'FLUZONE':2, 'Saline':0,
                'Pneumovax vaccine group':2, 'Influenza vaccine group':2,
                'Flu':1, 'NS':0, 'Influenza':1}
    if (tn == 3):
        atype = self.h.getSurvName('c race')
        atypes = ['W', 'B', 'A']
        ahash = {'Caucasian':0, 'Asian':2, 'African American':1}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getObermoser2013 = getObermoser2013


def getMahajan2016(self, tn=1):
    self.prepareData("COV261")
    atype = self.h.getSurvName("c condition")
    atypes = ['C', 'SBI', 'nSBI']
    ahash = {'nonSBI':2, 'SBI':1, 'Healthy Control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['C', 'I']
        ahash = {'nonSBI':1, 'SBI':1, 'Healthy Control':0}
    if (tn == 3):
        atype = self.h.getSurvName('c race')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'Black or African American':1, 'White':0, 'Stated as Unknown':3,
                'Asian':2, 'Other':3, 'American Indian or Alaska Native':3,
                'Native Hawiian or Other Pacific Islander':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMahajan2016 = getMahajan2016


def getChaussabel2008(self, tn=1):
    self.prepareData("COV163")
    atype = self.h.getSurvName("c Illness")
    atypes = ['HC', 'I']
    ahash = {'Healthy':0, 'SLE':1, 'UTI':1, 'Bacteremia':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'I']
        ahash = {'Healthy':0, 'SLE':1}
    if (tn == 3):
        atype = self.h.getSurvName('c Race')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'(Race Not Reported)':3, 'Black or African American':1,
                'White':0, 'Other':3, 'Asian':2, 'Selected More than One Race':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChaussabel2008 = getChaussabel2008


def getMiller2014(self, tn=1):
    self.prepareData("MACV262")
    atype = self.h.getSurvName("c sample group")
    atypes = ['HC', 'S']
    ahash = {'caregiver group':1, 'non-stressed control group':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c education")
        atypes = ['L', 'H', 'U', 'O']
        ahash = {'4':2, '1':0, '2':1, '6':3, '3':1, '5':2, 'NA':3, '0':0}
    if (tn == 3):
        atype = self.h.getSurvName('c ethnicity')
        atypes = ['W', 'O']
        ahash = {'Caucasian':0, 'Non-Caucasian':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMiller2014 = getMiller2014


def getYang2020(self, tn=1):
    self.prepareData("COV279")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['HC', 'IPF', 'CHP']
    ahash = {'chp':2, 'ipf':1, 'control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c Sex")
        atypes = ['F', 'M']
        ahash = {'male':1, 'female':0}
    if (tn == 3):
        atype = self.h.getSurvName('c race (hispanic;1, black;3,asian;4, white;5, other;6)')
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'5':0, '1':3, '3':1, '4':2, 'UNKNOWN':3, '6':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2020 = getYang2020


def getHuang2020(self, tn=1):
    self.prepareData("COV281")
    atype = self.h.getSurvName("c infection status")
    atypes = ['C', 'CoV']
    ahash = {'uninfected':0, 'infected with SARS-CoV-2 MOI 140':1}
    if (tn == 2):
        atype = self.h.getSurvName("c days post infection")
        ahash = {'NA':0, '1':1}
    if (tn == 3):
        atype = self.h.getSurvName("c days post infection")
        ahash = {'NA':0, '4':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHuang2020 = getHuang2020


def getYouk2020hAO(self, tn=1):
    self.prepareData("COV282")
    atype = self.h.getSurvName("c Timepoint")
    atypes = ['C', 'CoV']
    ahash = {'D0':0, 'D1':1, 'D3':1}
    if (tn == 2):
        ahash = {'D0':0, 'D1':1}
    if (tn == 3):
        ahash = {'D0':0, 'D3':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYouk2020hAO = getYouk2020hAO


def getYouk2020hBO(self, tn=1):
    self.prepareData("COV283")
    atype = self.h.getSurvName("c Timepoint")
    atypes = ['C', 'CoV']
    ahash = {'D0':0, 'D1':1, 'D3':1}
    if (tn == 2):
        ahash = {'D0':0, 'D1':1}
    if (tn == 3):
        ahash = {'D0':0, 'D3':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYouk2020hBO = getYouk2020hBO


def getVlasovaStLouis2018(self, tn=1):
    self.prepareData("COV290")
    atype = self.h.getSurvName("c group description")
    atypes = ['C', 'E', 'L']
    ahash = {'earlyIRIS':1, 'na':0, 'LateIRIS':2}
    if (tn == 2):
        atype = self.h.getSurvName("c group")
        ahash = {'IRIS':1, 'control':0}
        atypes = ['C', 'I']
    if (tn == 3):
        ahash = {'earlyIRIS':1, 'na':0}
        atypes = ['C', 'I']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVlasovaStLouis2018 = getVlasovaStLouis2018


def getZaas2010Mm(self, tn=1):
    self.prepareData("COV294")
    atype = self.h.getSurvName("c infection duration")
    ahash = {'1':1, '2':2, '3':3, '4':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection status")
    atypes = ['H', 'B', 'F']
    ahash = {'candida':2, 'healthy':0, 'staph':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['H', 'F']
        ahash = {'candida':1, 'healthy':0}
        atype = [atype[i] if tval[i] == 2 or aval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['H', 'B']
        ahash = {'staph':1, 'healthy':0}
        atype = [atype[i] if tval[i] == 2 or aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZaas2010Mm = getZaas2010Mm


def getLi2019I(self, tn=1):
    self.prepareData("COV295")
    atype = self.h.getSurvName("c disease state")
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['Healthy', 'Fungal']
    ahash = {'healthy':0, 'patient':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2019I = getLi2019I


def getDix2015(self, tn=1):
    self.prepareData("COV296")
    atype = self.h.getSurvName("c infection")
    atypes = ['M', 'B', 'F']
    ahash = {'Escherichia coli':1, 'mock':0, 'Staphylococcus aureus':1,
            'Candida albicans':2, 'Aspergillus fumigatus':2}
    if (tn == 2):
        atypes = ['Mock', 'Fungal']
        ahash = {'mock':0, 'Candida albicans':1, 'Aspergillus fumigatus':1}
        ahash = {'mock':0, 'Aspergillus fumigatus':1}
        ahash = {'mock':0, 'Candida albicans':1}
    if (tn == 3):
        atypes = ['Mock', 'B']
        ahash = {'mock':0, 'Escherichia coli':1, 'Staphylococcus aureus':1}
        ahash = {'mock':0, 'Escherichia coli':1}
        ahash = {'mock':0, 'Staphylococcus aureus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDix2015 = getDix2015


def getSubramani2020Mm(self, tn=1):
    self.prepareData("MACV266")
    atype = self.h.getSurvName("c infection state")
    atypes = ['Mock', 'Fungal']
    ahash = {'mock-infected':0, 'infected':1, 'not infected':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSubramani2020Mm = getSubramani2020Mm


def getBruno2020(self, tn=1):
    self.prepareData("COV297")
    atype = self.h.getSurvName("c stimulus")
    atypes = ['Ctl', 'Fungal']
    ahash = {'calb_man':1, 'calb_bg':1, 'caur_KTClive':1,
            'calb_live':1, 'RPMI':0, 'caur_KTCbg':1, 'caur_KTCman':1}
    if (tn == 2):
        ahash = {'calb_live':1, 'RPMI':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBruno2020 = getBruno2020


def getDAiuto2015(self, tn=1):
    self.prepareData("COV298")
    atype = self.h.getSurvName("c source_name (ch1)")
    atype = [re.sub("[ -].*", "", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'HSV':1, 'uninfected':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDAiuto2015 = getDAiuto2015


def getTommasi2020(self, tn=1):
    self.prepareData("COV299")
    atype = self.h.getSurvName("c vzv infection (ch1)")
    atype = [re.sub("[ -].*", "", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'VZV':1, 'Uninfected':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTommasi2020 = getTommasi2020


def getHU2019(self, tn=1):
    self.prepareData("COV300")
    atype = self.h.getSurvName("c source_name (ch1)")
    atypes = ['C', 'I']
    ahash = {'CHB WTC':1, 'CHB YTC':1, 'Healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHU2019 = getHU2019


def getRiou2017(self, tn=1):
    self.prepareData("COV301")
    atype = self.h.getSurvName("c hcmv exposition (ch1)")
    atypes = ['C', 'I']
    ahash = {'No HCMV-specific IgG detected':0,
            'HCMV-specific IgG detected':0,
            'During acute primary HCMV infection':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRiou2017 = getRiou2017


def getZilliox2007(self, tn=1):
    self.prepareData("COV302")
    atype = self.h.getSurvName("c title2")
    atypes = ['C', 'I']
    ahash = {'Control':0, 'Patient':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZilliox2007 = getZilliox2007


def getOberstein2017(self, tn=1):
    self.prepareData("COV303")
    atype = self.h.getSurvName("c treatment (ch1)")
    atypes = ['C', 'I']
    ahash = {'mock':0, 'HCMV-infected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOberstein2017 = getOberstein2017


def getOliver2017(self, tn=1):
    self.prepareData("COV304")
    atype = self.h.getSurvName("c status (ch1)")
    atypes = ['C', 'I']
    ahash = {'Uninfected':0, 'infected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOliver2017 = getOliver2017


def getMarkus2014(self, tn=1):
    self.prepareData("COV305")
    atype = self.h.getSurvName("c treatment (ch1)")
    atype = [re.sub("[ -].*", "", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'infection':1, 'none':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMarkus2014 = getMarkus2014


def getYoon2019(self, tn=1):
    self.prepareData("COV306")
    atype = self.h.getSurvName("c ebv status (ch1)")
    atypes = ['C', 'I']
    ahash = {'EBV-negative':0, 'EBV-positive':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYoon2019 = getYoon2019


def getAkinci2020(self, tn=1, ta=0):
    self.prepareData("COV307")
    atype = self.h.getSurvName("c tissue")
    ahash = {'Colon':0, 'Liver':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atype = [re.sub(",.*", "", str(k)) for k in atype]
    atypes = ['C', 'HCQ', 'Rem']
    ahash = {'DMSO':0, 'HCQ':1, 'Remdesivir':2}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'Rem']
        ahash = {'DMSO':0, 'Remdesivir':1}
        atype = [atype[i] if tval[i] == ta
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAkinci2020 = getAkinci2020


def getThompson2017(self, tn=1):
    self.prepareData("COV308")
    atype = self.h.getSurvName('c disease state')
    atypes = ['H', 'TB', 'DxC', 'MTP']
    ahash = {'TB Subjects':1,
            'Healthy Controls':0,
            'Lung Dx Controls':2,
            'MTP Controls':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getThompson2017 = getThompson2017


def getLindeboom2020(self, tn=1):
    self.prepareData("COV317")
    atype = self.h.getSurvName('c icu admission')
    ahash = {'no ICU admission':1, 'ICU admission':2, 'NA':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['U', '0', '5', 'CA', 'CA+HCQ']
    ahash = {'0 days HCQ':1, '5 days HCQ':2, 'untreated':0, 'HKCA':3, 'HKCA+HCQ':4}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['C', 'HCQ']
        ahash = {'0 days HCQ':0, '5 days HCQ':1}
    if (tn == 3):
        atypes = ['U', 'CA', 'CA+HCQ']
        ahash = {'untreated':0, 'HKCA':1, 'HKCA+HCQ':2}
    if (tn == 4):
        atypes = ['U', 'CoV2']
        ahash = {'untreated':0, '0 days HCQ':1}
    if (tn == 5):
        atypes = ['U', 'no ICU', 'ICU']
        atype = [tval[i] if aval[i] == 0 or aval[i] == 1
                else None for i in range(len(atype))]
        ahash = {0:0, 1:1, 2:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLindeboom2020 = getLindeboom2020


def getBossel2019(self, tn=1):
    self.prepareData("COV318")
    atype = self.h.getSurvName('c time (post-infection)')
    atypes = ['0', '4', '8']
    ahash = {'8 hr':2, '0 hr':0, '4 hr':1}
    if (tn == 2):
        atypes = ['C', 'I']
        ahash = {'8 hr':1, '0 hr':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBossel2019 = getBossel2019


def getRobinson2020(self, tn=1):
    self.prepareData("COV319")
    atype = self.h.getSurvName('c infection')
    atypes = ['C', 'I']
    ahash = {'none':0, 'Salmonella typhimurium':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRobinson2020 = getRobinson2020


def getMcNeill2018(self, tn=1):
    self.prepareData("COV309")
    atype = self.h.getSurvName('c day')
    ahash = {'D28':28, 'D0':0}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c cmv status')
    ahash = {'CMVn':0, 'CMVp':1}
    mval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c risk for tb')
    atypes = ['C', 'I']
    ahash = {'Control':0, 'Case':1}
    atype = [atype[i] if dval[i] == 0 and mval[i] == 0
            else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMcNeill2018 = getMcNeill2018


def getKaforou2013(self, tn=1):
    self.prepareData("COV310")
    atype = self.h.getSurvName('c hiv status')
    ahash = {'HIV negative':0, 'HIV positive':1}
    hval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state')
    atypes = ['C', 'Tb', 'lTb']
    ahash = {'latent TB infection':2,
            'other disease':0,
            'active tuberculosis':1}
    atype = [atype[i] if hval[i] == 0
            else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKaforou2013 = getKaforou2013


def getAnderson2014(self, tn=1):
    self.prepareData("COV312")
    atype = self.h.getSurvName('c hiv status')
    ahash = {'HIV negative':0, 'HIV positive':1}
    hval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease status')
    atypes = ['C', 'Tb', 'lTb']
    ahash = {'latent TB infection':2,
            'other disease':0,
            'active tuberculosis':1}
    atype = [atype[i] if hval[i] == 0
            else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAnderson2014 = getAnderson2014


def getBartholomeus2019(self, tn=1):
    self.prepareData("COV313")
    atype = self.h.getSurvName('c infected with/healthy control')
    atypes = ['C', 'I']
    ahash = {'Control':0, 'Streptococcus pneumoniae':1, 
            'Haemophilus influenzae':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBartholomeus2019 = getBartholomeus2019


def getHerberg2016I(self, tn=1):
    self.prepareData("COV168")
    atype = self.h.getSurvName('c src1')
    atype = [re.sub("W.*from ", "", str(k)) for k in atype]
    atype = [re.sub("pa.*with ", "", str(k)) for k in atype]
    atypes = ['C', 'pB', 'pV', 'B', 'V', 'U']
    ahash = {'healthy control':0,
            'Probable Bacterial infection':1,
            'Probable Viral infection':2,
            'Definite Bacterial infection':3,
            'Definite Viral infection':4,
            'infection of uncertain bacterial or viral aetiology':5}
    if (tn == 2):
        atypes = ['C', 'B']
        ahash = {'healthy control':0, 'Definite Bacterial infection':1}
    if (tn == 3):
        atypes = ['V', 'B']
        ahash = {'Definite Viral infection':0,
                'Definite Bacterial infection':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHerberg2016I = getHerberg2016I


def getLieberman2020(self, tn=1):
    self.prepareData("COV275")
    atype = self.h.getSurvName('c src1')
    ahash = {'Nasopharyngeal Swab':0, 'Human Airway Epithelial Cells':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c src1')
    ahash = {'Nasopharyngeal Swab':0, 'Human Airway Epithelial Cells':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c sars-cov-2 infection')
    ahash = {'infected':1, 'uninfected':0}
    atype = self.h.getSurvName('c sars-cov-2 positivity')
    atypes = ['C', 'I']
    ahash = {'pos':1, 'neg':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLieberman2020 = getLieberman2020


def getPrice2020Cov2(self, tn=1):
    self.prepareData("COV278")
    atype = self.h.getSurvName('c day post-infection')
    ahash = {'Day 0':0, 'Day 1':1, 'Day 3':3, 'Day 5':5, 'Day 7':7,
            'Day 10':10, 'Day 12':12, 'Day 15':15, 'Day 17':17}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c agent')
    atypes = ['C', 'CoV2']
    ahash = {'Control':0, 'SARS-CoV-2':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 or tval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = [0, 1, 3, 5, 7, 10, 12, 15, 17]
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPrice2020Cov2 = getPrice2020Cov2


def getFiege2020(self, tn=1):
    self.prepareData("COV320")
    atype = self.h.getSurvName('c treatment')
    ahash = {'untreated':0, 'Remdesivir pretreated':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c infection status')
    atypes = ['C', 'CoV2']
    ahash = {'uninfected':0,
            'SARS-CoV-2 infected, MOI 5, 24 hours post infection':1,
            'SARS-CoV-2 infected, MOI 5, 48 hours post infection':1}
    if (tn == 2):
        atypes = ['C', 'CoV2', 'R']
        atype = ['R' if tval[i] == 1
                else atype[i] for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFiege2020 = getFiege2020


def getGajMSD(self, tn=1):
    self.prepareData("COV257.2")
    atype = self.h.getSurvName('c Disease_stage')
    ahash = {'':0, 'Acute':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Diagnosis')
    ahash = {'':0, 'KD':1, 'MIS-C':2}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'I', 'S']
    ahash = {'1y Convalescent (CAA+)':0, '1y Convalescent (CAA-)':0,
            'Acute (CAA+)':1, 'Acute(CAA-)':1, 'MIS-C':1,
            'Healthy Adult':0, 'Acute COVID19':1, 'Standards':2}
    if (tn == 2):
        atype = dval
        ahash = {1:0, 2:1}
        atypes = ['KD', 'MISC']
        atype = [atype[i] if sval[i] == 1
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['CV', 'AV']
        ahash = {'1y Convalescent (CAA+)':0, '1y Convalescent (CAA-)':0,
                'Acute (CAA+)':1, 'Acute(CAA-)':1}
    if (tn == 4):
        atypes = ['CV', 'AV', 'M']
        ahash = {'1y Convalescent (CAA+)':0, '1y Convalescent (CAA-)':0,
                'Acute (CAA+)':1, 'Acute(CAA-)':1, 'MIS-C':2}
    if (tn == 5):
        atypes = ['AV', 'M']
        ahash = {'Acute (CAA+)':0, 'Acute(CAA-)':0, 'MIS-C':1}
    if (tn == 6):
        atype = self.h.getSurvName('c Group')
        atypes = ['C', 'CoV']
        ahash = {'Healthy Adult':0, 'Acute COVID19':1}
    if (tn == 7):
        atype = self.h.getSurvName('c Gender')
        ahash = {'F':0, 'M':1}
        sval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.h.getSurvName('c Severity')
        atypes = ['NC', 'Critical']
        ahash = {'severe':1, 'critical':1, 'moderate':0,
                'moderate to severe':0, 'fatal':1, 'Asymptomatic':0,
                'Fatal':1, 'Severe':1, 'Critical':1, 'Mod-Severe':0,
                'Moderate':0}
        atype = [atype[i] if sval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGajMSD = getGajMSD


def getBurns2020KDMISC(self, tn=1):
    self.prepareData("COV257.3")
    atype = self.h.getSurvName('c Disease_stage')
    ahash = {'Acute':1, 'Convalescent':0, '':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Diagnosis')
    atypes = ['K', 'M', 'F', 'U']
    ahash = {'KD':0, '':3, 'FC':2, 'MIS-C':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['CV', 'AV', 'M']
        atype = ['CV' if tval[i] == 0 and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if tval[i] == 1 and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 3):
        atypes = ['CAA-', 'CAA+']
        atype = self.h.getSurvName('c CAA_pos_neg')
        ahash = {'pos':1, 'neg':0}
    if (tn == 4):
        atypes = ['C', 'A1', 'A2', 'A3', 'A4']
        btype = self.h.getSurvName("c CA status")
        atype = self.h.getSurvName('c Disease_stage')
        ahash = {'Acute':1, 'Convalescent':0, '':2}
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'4':4, '3':3, '2':2, '1':1, '':-1, 'na': -1}
        atype = [ahash[btype[i]] if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        ahash = {'Convalescent':0, 'Acute':1, 0:0, 1:1, 2:2, 3:3, 4:4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurns2020KDMISC = getBurns2020KDMISC


def getBurns2020KDMISCII(self, tn=1):
    self.prepareData("COV257.4")
    atype = self.h.getSurvName('c Treatment with statin')
    ahash = {'No':0, 'Yes':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease_status')
    atypes = ['K', 'M']
    ahash = {'KD':0, 'MIS-C':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['K', 'M', 'D']
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
        atype = ['D' if str(self.h.headers[i]) == 'S204026'
                else atype[i] for i in range(len(atype))]
    if (tn == 3):
        atypes = ['K', 'M']
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
        atype = [None if str(self.h.headers[i]) == 'S204026'
                else atype[i] for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurns2020KDMISCII = getBurns2020KDMISCII


def getBurns2020KDMISCIII(self, tn=1, ta=0):
    self.prepareData("COV257.5")
    atype = self.h.getSurvName('c zworstever')
    zval = ['Z', 'c Z']
    for i in atype[2:]:
        if i == '' or i == 'na':
            zval.append(0)
        elif float(i) < 2:
            zval.append(1)
        elif float(i) >= 2.5 and float(i) < 10:
            zval.append(2)
        elif float(i) >= 10:
            zval.append(4)
        else:
            zval.append(0)
    iday = self.h.getSurvName("c illday")
    idval = ['i', 'c i']
    for i in iday[2:]:
        if i == '' or i == 'na':
            idval.append(0)
        elif float(i) <= 10:
            idval.append(1)
        elif float(i) > 10:
            idval.append(2)
        else:
            idval.append(0)
    atype1 = self.h.getSurvName('c sex')
    atype2 = self.h.getSurvName('c gender')
    atype = [atype1[i] if atype1[i] != '' else atype2[i]
            for i in range(len(atype1))]
    ahash = {'2':0, '1':1, '':2} #1=male, 2=female
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c batch')
    ahash = {'statin KD WB':2, 'MIS-C_UCSD':3, 'iLess10':0, 'UCL':1, '':4}
    bval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Treatment with statin')
    ahash = {'No':0, 'Yes':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Disease_stage')
    ahash = {'Acute ':2, 'Subacute':1, 'Acute':2, 'Convalescent':0,
            'Acute (Label changed)':2, '':3}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Diagnosis')
    atypes = ['K', 'M', 'FC']
    ahash = {'KD':0, 'MIS-C':1, '':0, 'FC':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['K', 'M', 'D']
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
        atype = ['D' if str(self.h.headers[i]) == 'S204026'
                else atype[i] for i in range(len(atype))]
    if (tn == 3):
        atypes = ['K', 'M']
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
        atype = [None if str(self.h.headers[i]) == 'S204026'
                else atype[i] for i in range(len(atype))]
    if (tn == 4):
        atypes = ['CV', 'SA', 'AV', 'M', 'FC']
        atype = ['CV' if dval[i] == 0 and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['SA' if dval[i] == 1 and aval[i] == 0 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and aval[i] == 0 and bval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        atype = ['FC' if aval[i] == 2
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 5):
        atypes = ['CAA-', 'CAA+']
        atype = self.h.getSurvName('c CAA_pos_neg')
        ahash = {'pos':1, 'neg':0}
    if (tn == 6):
        atypes = ['CV', 'CAA-', 'CAA+S', 'CAA+G']
        btype = self.h.getSurvName("c CA status")
        atype = self.h.getSurvName('c Disease_stage')
        ahash = {'Acute ':2, 'Subacute':1, 'Acute':2, 'Convalescent':0,
                'Acute (Label changed)':2, '':3}
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'4':4, '3':3, '2':2, '1':1, '':-1, 'na': -1}
        atype = [ahash[btype[i]] if aval[i] == 2
                else atype[i] for i in range(len(atype))]
        ahash = {'Convalescent':0, 'Subacute':1, 'Acute (Label changed)':1,
                'Acute':1, 'Acute ':1, 0:0, 1:1, 2:2, 4:3}
    if (tn == 7):
        atypes = ['CV', 'AV']
        atype = ['CV' if dval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and bval[i] == 0
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 8):
        atypes = ['SA', 'AV']
        atype = ['SA' if dval[i] == 1 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and tval[i] is not None and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 9):
        atypes = ['SA', 'AV', 'ST', 'M']
        atype = ['SA' if dval[i] == 1 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and tval[i] is not None and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['ST' if dval[i] == 1 and tval[i] == 1 and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 10):
        atypes = ['CV', 'CAA-', 'CAA+S', 'CAA+G']
        atype = self.h.getSurvName('c Disease_stage')
        ahash = {'Acute ':2, 'Subacute':1, 'Acute':2, 'Convalescent':0,
                'Acute (Label changed)':2, '':3}
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = [zval[i] if aval[i] == 2 and idval[i] == 1 and bval[i] == 0
                else atype[i] for i in range(len(atype))]
        ahash = {'Convalescent':0, 1:1, 2:2, 4:3}
    if (tn == 11):
        atypes = ['SA', 'AV', 'M']
        atype = ['SA' if dval[i] == 1 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and tval[i] is not None
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 12):
        atypes = ['SA', 'AV', 'M']
        atype = ['SA' if dval[i] == 1 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and tval[i] is not None
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        atype = [None if str(self.h.headers[i]) == 'S204026'
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 13):
        atypes = ['K', 'M']
        atype = ['K' if dval[i] == 2 and tval[i] is not None
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 14):
        atypes = ['CAA+S', 'CAA+G']
        atype = self.h.getSurvName('c Disease_stage')
        ahash = {'Acute ':2, 'Subacute':1, 'Acute':2, 'Convalescent':0,
                'Acute (Label changed)':2, '':3}
        aval = [ahash[i] if i in ahash else None for i in atype]
        atype = [zval[i] if aval[i] == 2 and idval[i] == 1 and bval[i] == 0
                else atype[i] for i in range(len(atype))]
        ahash = {2:0, 4:1}
    if (tn == 15):
        atypes = ['SA', 'AV', 'M', 'FC']
        atype = ['SA' if dval[i] == 1 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and aval[i] == 0 and bval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and tval[i] is not None
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        atype = ['FC' if aval[i] == 2
                else atype[i] for i in range(len(atype))]
        ahash = {}
    if (tn == 16):
        atypes = ['CV', 'SA', 'AV', 'M', 'FC']
        atype = ['CV' if dval[i] == 0 and aval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['SA' if dval[i] == 1 and aval[i] == 0 and tval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['AV' if dval[i] == 2 and aval[i] == 0 and bval[i] == 0
                else atype[i] for i in range(len(atype))]
        atype = ['M' if aval[i] == 1
                else atype[i] for i in range(len(atype))]
        atype = ['FC' if aval[i] == 2
                else atype[i] for i in range(len(atype))]
        ahash = {}
        atype = [atype[i] if gval[i] == ta
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurns2020KDMISCIII = getBurns2020KDMISCIII


def getBurns2020KDMISCIV(self, tn=1):
    self.prepareData("COV257.7")
    atype = self.h.getSurvName('c Diagnosis')
    atypes = ['K', 'M']
    ahash = {'KD':0, 'MIS-C':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBurns2020KDMISCIV = getBurns2020KDMISCIV


def getManne2020(self, tn=1):
    self.prepareData("COV321")
    atype = self.h.getSurvName('c disease_stage')
    atypes = ['C', 'I', 'S']
    ahash = {'ICU':2, 'Non-ICU':1, 'Healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getManne2020 = getManne2020


def getLacerdaMariano2020Mm(self, tn=1):
    self.prepareData("MACV267")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("-R.:.*", "", str(k)) for k in atype]
    atypes = ['N-L', 'N-H', 'I-L', 'I-H']
    ahash = {'naive-L':0, 'naive-H':1, 'inf-L':2, 'inf-H':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLacerdaMariano2020Mm = getLacerdaMariano2020Mm


def getMatkovich2017(self, tn=1):
    self.prepareData("MACV268")
    atype = self.h.getSurvName('c condition')
    atypes = ['N', 'S', 'ICM', 'NICM']
    ahash = {'septic cardiomyopathy':1,
            'ischemic heart disease':2,
            'nonischemic dilated cardiomyopathy':3,
            'nonfailing heart':0}
    if (tn == 2):
        atypes = ['N', 'S']
        ahash = {'septic cardiomyopathy':1,
                'nonfailing heart':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatkovich2017 = getMatkovich2017


def getCoulibaly2019Gr(self, tn=1):
    self.prepareData("MACV269")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['PS', 'SS', 'SIRS']
    ahash = {'SIRS':2, 'septic shock':1, 'presurgical':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCoulibaly2019Gr = getCoulibaly2019Gr


def getCoulibaly2019NK(self, tn=1):
    self.prepareData("MACV270")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['PS', 'SS', 'SIRS']
    ahash = {'SIRS':2, 'septic shock':1, 'presurgical':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCoulibaly2019NK = getCoulibaly2019NK


def getPG2020LungHamNamir(self, tn=1):
    self.prepareData("COV323")
    atype = self.h.getSurvName('c info')
    atypes = ['U', 'V', 'N', '3', '4']
    ahash = {'Merck Drug':2, 'Merck Drug Vehicle Control':1, 'UN':0}
    if (tn == 2):
        atypes = ['U', 'V', 'N']
        ahash = {'Merck Drug':2, 'Merck Drug Vehicle Control':1, 'UN':0}
        ah = {'3', '4'}
        atype = [None if k in ah else k for k in atype]
    if (tn == 3):
        atypes = ['U', 'V']
        ahash = {'Merck Drug Vehicle Control':1, 'UN':0}
        ah = {'3', '4'}
        atype = [None if k in ah else k for k in atype]
    if (tn == 4):
        atypes = ['V', 'N']
        ahash = {'Merck Drug':1, 'Merck Drug Vehicle Control':0}
        ah = {'3', '4'}
        atype = [None if k in ah else k for k in atype]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020LungHamNamir = getPG2020LungHamNamir


def getThair2020(self, tn=1, tb=0):
    self.prepareData("COV327")
    atype = self.h.getSurvName('c disease')
    atypes = ['H', 'CoV']
    ahash = {'Healthy control':0, 'COVID19':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getThair2020 = getThair2020


def getBernardes2020(self, tn=1, tb=0):
    self.prepareData("COV328")
    atype = self.h.getSurvName('c remission')
    atypes = ['H', 'R', 'NR']
    ahash = {'Remission':1, 'Healthy':0, 'No Remission':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBernardes2020 = getBernardes2020


def getMaia2020(self, tn=1, tb=0):
    self.prepareData("COV336")
    atype = self.h.getSurvName('c src1')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'Classical':0, 'Non-classical':1, 'Plasmacytoid':2,
            'Neutrophils':3, 'Basophils':4, 'Myeloid':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c haematological tumor')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    ahash = {'No':0, 'Monoclonal':1, 'Acute':2, 'Diffuse':3}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c covid status')
    atypes = ['R', 'I']
    ahash = {'Infected':1, 'Recovered':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaia2020 = getMaia2020


def getParnell2013(self, tn=1, tb=0):
    self.prepareData("MACV272")
    atype = self.h.getSurvName('c disease status')
    atypes = ['H', 'S', 'NS']
    ahash = {'healthy':0, 'sepsis survivor':1, 'sepsis nonsurvivor':2}
    if (tn == 2):
        atypes = ['H', 'S']
        ahash = {'healthy':0, 'sepsis survivor':1, 'sepsis nonsurvivor':1}
    if (tn == 3):
        atypes = ['S', 'NS']
        ahash = {'sepsis survivor':0, 'sepsis nonsurvivor':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getParnell2013 = getParnell2013


def getParedes2020CRC(self, tn=1, tb=0):
    self.prepareData("CRC148")
    atype = self.h.getSurvName('c individual')
    atype = [re.sub("_.*", "", str(k)) for k in atype]
    ahash = {'Caucasian':0, 'AfricanAmerican':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c tissue')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['N', 'T']
    ahash = {'Tumor':1, 'Adjacent':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'AA']
        atype = tval
        ahash = {0:0, 1:1}
        atype = [atype[i] if aval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getParedes2020CRC = getParedes2020CRC


def getKim2020VV(self, tn=1, tb=0):
    self.prepareData("MACV273")
    atype = self.h.getSurvName('c Title')
    time = [str(k).split("_")[2] if len(str(k).split("_")) > 2
                     else None for k in atype]
    ahash = {'0h':0, '3h':3, '6h':6}
    tval = [ahash[i] if i in ahash else None for i in time]
    host = [str(k).split("_")[0] if len(str(k).split("_")) > 2
                     else None for k in atype]
    ahash = {'dTHP-1':0, 'HT-29':1}
    hval = [ahash[i] if i in ahash else None for i in host]
    atype = self.h.getSurvName('c src1')
    atype = [re.sub(" ce.*", "", str(k)) for k in atype]
    atypes = ['M', 'I']
    ahash = {'Vibrio-infected dTHP-1':1, 'Mock-treated dTHP-1':0,
            'Vibrio-infected HT-29':1, 'Mock-treated HT-29':0}
    ival = [ahash[i] if i in ahash else None for i in atype]
    atypes = [0, 3, 6]
    atype = tval
    atype = [atype[i] if hval[i] == 0
            else None for i in range(len(atype))]
    atype = [0 if ival[i] == 0
            else atype[i] for i in range(len(atype))]
    ahash= {}
    if (tn == 2 or tn == 4):
        atypes = [0, 3, 6]
        atype = tval
        atype = [atype[i] if hval[i] == 1
                else None for i in range(len(atype))]
        atype = [0 if ival[i] == 0
                else atype[i] for i in range(len(atype))]
    if (tn == 3 or tn == 4):
        atypes = [0, 6]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKim2020VV = getKim2020VV


def getStathopoulos2009(self, tn=1, tb=0):
    self.prepareData("MACV274")
    atype = self.h.getSurvName('c src1')
    atypes = ['H', '1', '2']
    ahash = {'blood, two primary malignancies':2,
            'blood, healthy':0, 'blood, one primary malignancy':1}
    if (tn == 2):
        atypes = ['H', 'T']
        ahash = {'blood, two primary malignancies':1,
                'blood, healthy':0, 'blood, one primary malignancy':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getStathopoulos2009 = getStathopoulos2009


def getGarridoMartin2020NSCLC(self, tn=1, tb=0):
    self.prepareData("MACV275")
    atype = self.h.getSurvName('c cells')
    atypes = ['M', 'TAM']
    ahash = {'Macrophages':0, 'Tumour Associated Macrophages':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGarridoMartin2020NSCLC = getGarridoMartin2020NSCLC


def getAltman2020Blood(self, tn=1, tb=0):
    self.prepareData("MACV256")
    atype = self.h.getSurvName('c src1')
    disease = [str(k).split("-")[1] if len(str(k).split("-")) > 1
                     else None for k in atype]
    stype = [str(k).split("-")[2] if len(str(k).split("-")) > 2
                     else None for k in atype]
    atypes = ['C', 'COPD']
    ahash = {'whole blood-COPD-Control':0,
            'whole blood-COPD-COPD':1}
    if (tn == 2):
        atypes = ['C', 'Staph']
        ahash = {'whole blood-Staph-Control':0,
                'whole blood-Staph-Staph':1}
    if (tn == 3):
        atypes = ['C', 'Sepsis']
        ahash = {'whole blood-Sepsis-Control':0,
                'whole blood-Sepsis-melioidosis':1}
    if (tn == 4):
        atypes = ['C', 'TB']
        ahash = {'whole blood-TB-Control':0,
                'whole blood-TB-PTB':1}
    if (tn == 5):
        atypes = ['C', 'Melanoma']
        ahash = {'whole blood-Melanoma-Control':0,
                'whole blood-Melanoma-Melanoma':1}
    if (tn == 6):
        atypes = ['C', 'Bcell def']
        ahash = {'whole blood-B-cell deficiency-Bcell':1,
                'whole blood-B-cell deficiency-Control':0}
    if (tn == 7):
        atypes = ['C', 'Flu']
        ahash = {'whole blood-Flu-Control':0,
                'whole blood-Flu-FLU':1}
    if (tn == 8):
        atypes = ['C', 'HIV']
        ahash = {'whole blood-HIV-Control':0,
                'whole blood-HIV-HIV':1}
    if (tn == 9):
        atypes = ['C', 'JDM']
        ahash = {'whole blood-Juvenile Dermatomyositis-Control':0,
                'whole blood-Juvenile Dermatomyositis-JDM':1}
    if (tn == 10):
        atypes = ['C', 'KD']
        ahash = {'whole blood-Kawasaki-Control':0,
                'whole blood-Kawasaki-Kawasaki':1}
    if (tn == 11):
        atypes = ['C', 'Liver Transplant']
        ahash = {'whole blood-Liver Transplant-Control':0,
                'whole blood-Liver Transplant-Transplant':1}
    if (tn == 12):
        atypes = ['C', 'MS']
        ahash = {'whole blood-MS-Control':0,
                'whole blood-MS-MS Patient':1}
    if (tn == 13):
        atypes = ['C', 'Pregnancy']
        ahash = {'whole blood-Pregnancy-Control':0,
                'whole blood-Pregnancy-Pregnancy':1}
    if (tn == 14):
        atypes = ['C', 'RSV']
        ahash = {'whole blood-RSV-Control':0,
                'whole blood-RSV-RSV':1}
    if (tn == 15):
        atypes = ['C', 'SLE']
        ahash = {'whole blood-SLE-Control':0,
                'whole blood-SLE-SLE':1}
    if (tn == 16):
        atypes = ['C', 'SoJIA']
        ahash = {'whole blood-SoJIA-Control':0,
                'whole blood-SoJIA-SoJIA':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAltman2020Blood = getAltman2020Blood


def getBadea2008Panc(self, tn=1, tb=0):
    self.prepareData("PANC7")
    atype = self.h.getSurvName("c type")
    atypes = ['N', 'T']
    ahash = {'normal':0, 'tumor':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBadea2008Panc = getBadea2008Panc


def getYang2016PDAC(self, tn=1, tb=0):
    self.prepareData("PANC15")
    atype = self.h.getSurvName("c grading")
    ahash = {'G3':3, 'G2':2, 'G4':4, 'Gx':None, 'G1':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c tissue")
    atypes = ['N', 'T']
    ahash = {'Pancreatic tumor':1, 'adjacent pancreatic non-tumor':0}
    if (tn == 2):
        atype = [atype[i] if gval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2016PDAC = getYang2016PDAC


def getKirby2016(self, tn=1, tb=0):
    self.prepareData("PANC14")
    atype = self.h.getSurvName("c src1")
    atypes = ['CL', 'T']
    ahash = {'pancreatic adenocarcinoma cancer tissue':1,
            'Pancreatic cancer cell line':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKirby2016 = getKirby2016


def getTCGAPAAD(self, tn=1, tb=0):
    self.prepareData("PANC13")
    atype = self.h.getSurvName("c Histology")
    atypes = ['N', 'T']
    ahash = {'Primary Tumor':1, 'Solid Tissue Normal':0}
    if tn == 2:
        atypes = ['N', 'T', 'M']
        ahash = {'Primary Tumor':1, 'Solid Tissue Normal':0, 'Metastatic':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTCGAPAAD = getTCGAPAAD


def getZhang2012PDAC(self, tn=1, tb=0):
    self.prepareData("PANC11")
    atype = self.h.getSurvName("c sample type")
    atypes = ['N', 'T']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2012PDAC = getZhang2012PDAC


def getPei2009PDAC(self, tn=1, tb=0):
    self.prepareData("PANC19")
    atype = self.h.getSurvName("c tissue")
    atypes = ['N', 'T']
    ahash = {'Tumor Tissue in Pancreatic Cancer Sample':1,
            'Normal Tissue in Pancreatic Cancer Sample':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPei2009PDAC = getPei2009PDAC


def getBalagurunathan2008PDAC(self, tn=1, tb=0):
    self.prepareData("PANC4")
    atype = self.h.getSurvName("c sample type")
    atypes = ['N', 'T']
    ahash = {'normal tissue':0, 'primary pancreatic tumor':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBalagurunathan2008PDAC = getBalagurunathan2008PDAC


def getJimeno2008PDAC(self, tn=1, tb=0):
    self.prepareData("PANC6")
    atype = self.h.getSurvName("c Classification")
    atypes = ['S', 'R']
    ahash = {'Sensitive':0, 'Resistant':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJimeno2008PDAC = getJimeno2008PDAC


def getIshikawa2005PDAC(self, tn=1, tb=0):
    self.prepareData("PANC8.1")
    atype = self.h.getSurvName("c atypicalCellProportion")
    atypes = ['L', 'H']
    ahash = {}
    if (tn == 2):
        atypes = ['L', 'M', 'H']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getIshikawa2005PDAC = getIshikawa2005PDAC


def getPeral2021PDAC(self, tn=1, tb=0):
    self.prepareData("PANC17")
    atype = self.h.getSurvName("c src1")
    atypes = ['WT', 'G', 'Het', 'HetG']
    ahash = {'Wildtype control_PDAC':0,
            'CDH11 Heterozygous control_PDAC':2,
            'Wildtype Gemcitabine_PDAC':1,
            'CDH11 Heterozygous Gemcitabine_PDAC':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPeral2021PDAC = getPeral2021PDAC


def getYu2020PDACblood(self, tn=1, tb=0):
    self.prepareData("PANC21")
    atype = self.h.getSurvName("c disease state")
    atypes = ['N', 'T']
    ahash = {'healthy':0, 'PDAC':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYu2020PDACblood = getYu2020PDACblood


def getMoffitt2015PDAC(self, tn=1, tb=0):
    self.prepareData("PANC22")
    atype = self.h.getSurvName('c cell line/tissue')
    ahash = {'Pancreas':0, 'Lung':1, 'Spleen':2, 'Liver':3,
            'Peritoneal':4, 'Colon':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c tissue type")
    atypes = ['N', 'T']
    ahash = {'Primary':1, 'Normal':0}
    atype = [atype[i] if tval[i] == 0
            else None for i in range(len(atype))]
    if (tn == 2):
        atypes = ['N', 'T', 'M']
        ahash = {'Primary':1, 'Metastasis':2, 'Normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMoffitt2015PDAC = getMoffitt2015PDAC


def getMaurer2019PDAC(self, tn=1, tb=0):
    self.prepareData("PANC23")
    atype = self.h.getSurvName("c compartment")
    atypes = ['S', 'E']
    ahash = {'Epithelium':1, 'Stroma':0}
    if (tn == 2):
        atypes = ['S', 'E', 'B']
        ahash = {'Epithelium':1, 'Stroma':0, 'Bulk':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaurer2019PDAC = getMaurer2019PDAC


def getPommier2018PDACmm(self, tn=1, tb=0):
    self.prepareData("PANC24")
    atype = self.h.getSurvName("c cell subpopulation")
    atypes = ['E-', 'E+']
    ahash = {'Ecad+':1, 'Ecad-':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPommier2018PDACmm = getPommier2018PDACmm


def getJanky2016PDAC(self, tn=1, tb=0):
    self.prepareData("PANC25")
    atype = self.h.getSurvName("c tissue")
    atypes = ['N', 'T']
    ahash = {'pancreatic tumor':1, 'non-tumoral pancreatic tissue':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJanky2016PDAC = getJanky2016PDAC


def getRashid2020PDAC(self, tn=1, tb=0):
    self.prepareData("PANC28")
    atype = self.h.getSurvName("c sample type")
    atypes = ['FNA', 'FFPE', 'FF', 'TB']
    ahash = {'FNA':0, 'FFPE':1, 'FF':2, 'tumor biopsies':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRashid2020PDAC = getRashid2020PDAC


def getRamaswamy2021MISC(self, tn=1, tb=0):
    self.prepareData("COV337")
    atype = self.h.getSurvName("c disease subtype")
    ahash = {'Severe':4, 'Severe; Recovered':3, '--':0, 'Moderate':2,
             'Moderate; Recovered':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c disease state")
    atypes = ['H', 'M']
    ahash = {'Multisystem inflammatory syndrome in children (MIS-C)':1, 'healthy':0}
    if (tn == 2):
        atype = tval
        atypes = ['H', 'M', 'S']
        ahash = {0:0, 2:1, 4:2}
    if (tn == 3):
        atype = tval
        atypes = ['H', 'M', 'S']
        ahash = {0:0, 1:0, 3:0, 2:1, 4:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRamaswamy2021MISC = getRamaswamy2021MISC


def getBrunetta2021CoV2(self, tn=1, tb=0):
    self.prepareData("COV338")
    atype = self.h.getSurvName("c condition")
    atypes = ['H', 'CoV']
    ahash = {'Healthy Control individual':0, 'COVID-19 hospitalized patient':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBrunetta2021CoV2 = getBrunetta2021CoV2


def getXu2020CoV2(self, tn=1, tb=0):
    self.prepareData("COV339")
    atype = self.h.getSurvName("c disease condition")
    atypes = ['H', 'CoV', 'Hp', 'IPF', 'Ma', 'Ssa']
    ahash = {'Hypersensitivity pneumonitis':2,
            'Donor':0,
            'Idiopathic pulmonary fibrosis':3,
            'Myositis-associated interstitial lng disease':4,
            'Systemic slcerosis-associated interstitial lung disease':5, '':1}
    if (tn == 2):
        atypes = ['H', 'CoV', 'IPF']
        ahash = {'Donor':0, 'Idiopathic pulmonary fibrosis':2, '':1}
    if (tn == 3):
        atypes = ['H', 'CoV']
        ahash = {'Donor':0, '':1}
    if (tn == 4):
        atypes = ['H', 'IPF']
        ahash = {'Donor':0, 'Idiopathic pulmonary fibrosis':1}
    if (tn == 5):
        atypes = ['IPF', 'CoV']
        ahash = {'Idiopathic pulmonary fibrosis':0, '':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXu2020CoV2 = getXu2020CoV2


def getGeng2019IPF(self, tn=1, tb=0):
    self.prepareData("COV341")
    atype = self.h.getSurvName("c cell type")
    atypes = ['NI', 'I']
    ahash = {'Invasive cells':1, 'non-Invasive cells':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGeng2019IPF = getGeng2019IPF


def getYao2021IPF(self, tn=1, tb=0):
    self.prepareData("COV343")
    atype = self.h.getSurvName("c group")
    atypes = ['H', 'IPF']
    ahash = {'Donor':0, 'IPF':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2021IPF = getYao2021IPF


def getBharat2020CoV2(self, tn=1, tb=1):
    dbid = 'COV342.1'
    if tn == 2:
        dbid = 'COV342.2'
    elif tn == 3:
        dbid = 'COV342.3'
    elif tn == 4:
        dbid = 'COV342.4'
    self.prepareData(dbid)
    atype = self.h.getSurvName("c Tissue Type")
    ahash = {'Parenchyma':0, 'Biopsy':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Diagnosis")
    atypes = ['H', 'CoV', 'IPF']
    ahash = {'COVID-19':1, 'Control (B.)':0, 'Control (H.)':0,
            'IPF':2, 'Other PF':2}
    if (tb == 2):
        atypes = ['H', 'CoV']
        atype = [atype[i] if tval[i] == 0
                else None for i in range(len(atype))]
    if (tb == 3):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBharat2020CoV2 = getBharat2020CoV2


def getBharat2020CoV2scblk(self, tn=1, tb=1):
    self.prepareData('COV342.5')
    atype = self.h.getSurvName("c cell type")
    ahash = {'Sorted Stromal population':2,
            'Sorted Myeloid population':1,
            'Sorted CD31 population':3,
            'Sorted Epithelial population':0,
            'Sorted Stromal and Myeloid populations':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c covid-19")
    atypes = ['H', 'CoV']
    ahash = {'Negative':0, 'Positive':1}
    if (tn == 2):
        atypes = ['H', 'CoV']
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBharat2020CoV2scblk = getBharat2020CoV2scblk


def getBoesch2020ipf(self, tn=1, tb=1):
    self.prepareData('COV344')
    atype = self.h.getSurvName("c condition")
    atypes = ['C', 'IPF']
    ahash = {'idiopathic pulmonary fibrosis (IPF)':1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBoesch2020ipf = getBoesch2020ipf


def getLuo2021ildRat(self, tn=1, tb=1):
    self.prepareData('COV345')
    atype = self.h.getSurvName("c bleomycin")
    atypes = ['C', 'Bleomycin']
    ahash = {'with':1, 'without':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLuo2021ildRat = getLuo2021ildRat


def getBauer2015ipfRat(self, tn=1, tb=1):
    self.prepareData('COV346')
    atype = self.h.getSurvName("c time")
    ahash = {'8 WEEKS':8, '2 WEEKS':2, '6 WEEKS':6, '3 DAYS':0,
            '1 WEEKS':1, '4 WEEKS':4, '3 WEEKS':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'Bleomycin']
    ahash = {'VEHICLE':0, 'BLEOMYCIN':1, 'NAIVE (UNTREATED)':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBauer2015ipfRat = getBauer2015ipfRat


def getGuillotin2021ipf(self, tn=1, tb=1):
    self.prepareData('COV347')
    atype = self.h.getSurvName("c source")
    atypes = ['Biopsy', 'transplant']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuillotin2021ipf = getGuillotin2021ipf


def getMeltzer2011ipf(self, tn=1, tb=1):
    self.prepareData('COV348')
    atype = self.h.getSurvName("c gender")
    ahash = {'female':0, 'male':2, 'unknown':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c phenotype")
    atypes = ['H', 'E', 'A']
    ahash = {'advanced idiopathic pulmonary fibrosis (IPF)':2,
            'early idiopathic pulmonary fibrosis (IPF)':1, 'healthy':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 2 or aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeltzer2011ipf = getMeltzer2011ipf


def getYao2020IPF(self, tn=1, tb=1):
    self.prepareData('COV268')
    atype = self.h.getSurvName("c group")
    atypes = ['control', 'IPF']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2020IPF = getYao2020IPF


def getYao2020IPFmm(self, tn=1, tb=1):
    self.prepareData('COV268.2')
    atype = self.h.getSurvName("c treatment")
    atypes = ['control', 'tamoxifen']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2020IPFmm = getYao2020IPFmm


def getYao2020Lung(self, tn=1, tb=1):
    self.prepareData('COV349')
    atype = self.h.getSurvName("c cd66 status")
    atypes = ['CD66-', 'CD66+']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2020Lung = getYao2020Lung


def getYao2020scblkIPFII(self, tn=1, tb=1):
    self.prepareData('COV349.2')
    atype = self.h.getSurvName("c group")
    atypes = ['Donor', 'IPF']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2020scblkIPFII = getYao2020scblkIPFII


def getSullivan2021a549(self, tn=1, tb=1):
    self.prepareData('COV350')
    atype = self.h.getSurvName('c treatment')
    ahash = {'none':0, 'Doxicycline':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c genotype")
    atypes = ['control', 'senescent']
    ahash = {'WT':0, 'TRF2-DN':1}
    atype = ['WT' if tval[i] == 0
            else atype[i] for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSullivan2021a549 = getSullivan2021a549


def getYao2020AT2CoV2(self, tn=1, tb=1):
    self.prepareData('COV367')
    atype = self.h.getSurvName("c group (ch1)")
    atypes = ['C', 'CoV2']
    ahash = {'Infected':1, 'mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2020AT2CoV2 = getYao2020AT2CoV2


def getYao2020AT2CoV2II(self, tn=1, tb=1):
    self.prepareData('COV369')
    atype = self.h.getSurvName("c group")
    atypes = ['C', 'CoV2']
    ahash = {'infected with SARS-CoV2':1, 'Mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYao2020AT2CoV2II = getYao2020AT2CoV2II


def getLamers2020CoV2(self, tn=1, tb=1):
    self.prepareData('COV326')
    atype = self.h.getSurvName("c cell type")
    ahash = {'Bronchioalveolar':0, 'Small airway':1, 'Lung bud tip':2,
            'Differentiating lung bud tip':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'CoV2']
    ahash = {'Mock':0, 'SARS-CoV-2':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLamers2020CoV2 = getLamers2020CoV2


def getKatsura2020at2CoV2(self, tn=1, tb=1):
    self.prepareData('COV370')
    atype = self.h.getSurvName("c infection status")
    atypes = ['C', 'CoV2']
    ahash = {'without infection':0, 'infected SARS-CoV-2':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKatsura2020at2CoV2 = getKatsura2020at2CoV2


def getPG2020iAT2(self, tn=1, tb=1):
    self.prepareData('COV371')
    atype = self.h.getSurvName("c type")
    atypes = ['C', '48', '72']
    ahash = {'un':0, '48h':1, '72h':2}
    if (tn == 2):
        atypes = ['C', 'CoV2']
        ahash = {'un':0, '48h':1, '72h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2020iAT2 = getPG2020iAT2


def getWinkler2021CoV2hACE2mm(self, tn=1, tb=1):
    self.prepareData('COV374')
    atype = self.h.getSurvName("c treatment")
    ahash = {'Naïve':0, 'LALA-PG D+2':2, 'LALA-PG D+1':3,
            '2050 D+2':4, '2050 D+1':5, 'Isotype':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c infection")
    atypes = ['C', 'CoV2']
    ahash = {'Naïve':0, 'SARS2':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0 or tval[i] == 1
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWinkler2021CoV2hACE2mm = getWinkler2021CoV2hACE2mm


def getZarrinpar2021aging(self, tn=1, tb=1):
    self.prepareData('AGE1')
    atype = self.h.getSurvName("c Sample code")
    ahash = {'D':0, 'BAT':1, 'TI':2, 'LIVER':3, 'WAT':4, 'SKM':5, 'CB':6, '':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c Sample code.1")
    atypes = ['NCD', 'HFD']
    ahash = {}
    if (tn == 2):
        atypes = ['ZT6-8', 'ZT9-12', 'ZT18-20', 'ZT21-24']
        ahash = {'ZT6':0, 'ZT7':0, 'ZT8':0, 'ZT9':1, 'ZT10':1, 'ZT11':1, 'ZT12':1,
                'ZT18':2, 'ZT19':2, 'ZT20':2, 'ZT21':3, 'ZT22':3, 'ZT23':3, 
                'ZT24':3}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZarrinpar2021aging = getZarrinpar2021aging


def getMelms2021CoV2snblk(self, tn=1, tb=1):
    self.prepareData('COV376')
    atype = self.h.getSurvName("c disease")
    atypes = ['C', 'CoV2']
    ahash = {'COVID-19':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMelms2021CoV2snblk = getMelms2021CoV2snblk


def getDelorey2021CoV2I(self, tn=1, tb=1):
    self.prepareData('COV377')
    atype = self.h.getSurvName("c src1")
    ahash = {'Parenchyma':0, 'LUL':1, 'blank sample':2, 'Trachea':3,
            'HeartLV':4, 'RUL':5, 'LLL':6}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c segment")
    atypes = ['C', 'CoV2']
    ahash = {'COVID-':0, 'COVID+':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelorey2021CoV2I = getDelorey2021CoV2I


def getDelorey2021CoV2II(self, tn=1, tb=1):
    self.prepareData('COV377.2')
    atype = self.h.getSurvName("c morphology")
    atypes = ['NA', 'IA', 'BE', 'A']
    ahash = {'Inflamed Alveoli':1, 'Artery':3,
            'Bronchial Epithelium':2, 'Normal Alveoli':0}
    if (tn == 2):
        atypes = ['NA', 'IA']
        ahash = {'Inflamed Alveoli':1, 'Normal Alveoli':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelorey2021CoV2II = getDelorey2021CoV2II


def getDelorey2021CoV2IV(self, tn=1, tb=1):
    self.prepareData('COV378')
    atype = self.h.getSurvName("c tissue")
    atypes = ['Lung', 'Brain', 'Heart']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelorey2021CoV2IV = getDelorey2021CoV2IV


def getDelorey2021CoV2scblk(self, tn=1, tb=1):
    self.prepareData('COV378.2')
    atype = self.h.getSurvName("c tissue")
    atypes = ['lung', 'kidney', 'heart', 'liver',
            'LN', 'spleen', 'airway', 'trachea', 'brain']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDelorey2021CoV2scblk = getDelorey2021CoV2scblk


def getMelillo2021CoV2scblk(self, tn=1, tb=1):
    self.prepareData('COV379')
    atype = self.h.getSurvName("c disease severity")
    atypes = ['Stable', 'Progressive']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMelillo2021CoV2scblk = getMelillo2021CoV2scblk


def getLangelier2021CoV2(self, tn=1, tb=1):
    self.prepareData('COV380')
    atype = self.h.getSurvName("c study group")
    atypes = ['E', 'L', 'T']
    ahash = {'No VAP - longitudinal':1, 'No VAP - late':2, 'No VAP - early':0,
            'VAP - early':0, 'VAP - longitudinal':1, 'VAP - late':2}
    if (tn == 2):
        ahash = {'No VAP - longitudinal':1, 'No VAP - late':2, 'No VAP - early':0}
    if (tn == 3):
        ahash = {'VAP - early':0, 'VAP - longitudinal':1, 'VAP - late':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLangelier2021CoV2 = getLangelier2021CoV2


def getZhao2021CoV2scblk (self, tn=1, tb=1):
    self.prepareData('COV381')
    atype = self.h.getSurvName("c cell type")
    ahash = {'CD3 positive':0, 'CD3 negative':1, 'CD45 negative':2,
            'EpCAM_positive':3, 'CD45 positive':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c disease")
    atypes = ['BN', 'CoV2']
    ahash = {'COVID19':1, 'Bacterial pneumonia':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhao2021CoV2scblk  = getZhao2021CoV2scblk 


def getBost2021CoV2scblk(self, tn=1, tb=1):
    self.prepareData('COV382')
    atype = self.h.getSurvName("c tissue")
    ahash = {'Blood':0, 'BAL':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c clinical outcome")
    atypes = ['H', 'M', 'S']
    ahash = {'NA':0, 'Dead':2, 'Alive':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'CoV2']
        ahash = {'NA':0, 'Dead':1, 'Alive':1}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBost2021CoV2scblk = getBost2021CoV2scblk


def getDesai2020CoV2(self, tn=1, tb=1):
    self.prepareData('COV383')
    atype = self.h.getSurvName("c tissue substructure")
    ahash = {'Alveoli':0, 'Bronchial':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c segment type")
    ahash = {'Geometric':0, 'PanCK_Pos':1, 'PanCK_Neg':2}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c sars-cov-2 rna ish")
    atypes = ['C', 'CoV2']
    ahash = {'Positive':1, 'Negative':0, 'NEARBY_Positive':1, 'Control':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDesai2020CoV2 = getDesai2020CoV2


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

def getWang2020CoV2scblkII(self, tn=1, tb=1):
    self.prepareData('COV384.2')
    cell = self.h.getSurvName("c celltype")
    atype = self.h.getSurvName("c age")
    atypes = ['N', 'C', 'A']
    ahash = {'3 yr':1, '29 wkGA':0, '31 wkGA':0,
            '31 yrs':2, '29 yrs':2, '33 yrs':2}
    if (tn == 2):
        atype = [atype[i] if cell[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        ctypes = ['Alveolar macrophages', 'monocytes']
        atype = [atype[i] if cell[i] in ctypes
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2020CoV2scblkII = getWang2020CoV2scblkII

def getRen2021CoV2scblk(self, tn=1, tb=1):
    self.prepareData('COV385')
    atype = self.h.getSurvName("c sample type")
    ahash = {'frozen PBMC':1, 'fresh PBMC':0,
            'CD19+ B cell sorted from fresh PBMC (FACS)':2,
            'CD3+ T cell sorted from fresh PBMC (FACS)':3,
            'fresh BALF':4, 'fresh PFMC':5, 'fresh Sputum':6,
            'CD3+ T cell and CD19+ B cell sorted from fresh PBMC (FACS)':7,
            'B cells sorted from frozen PBMC (MACS, STEMCELL 19054)':8,
            'CD19+ B cell sorted from fresh PBMC (MACS)':9}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c sample time")
    ahash = {'progression':2, 'convalescence':1, 'control':0}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c covid-19 severity")
    atypes = ['H', 'M', 'S']
    ahash = {'severe/critical':2, 'mild/moderate':1, 'control':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName("c sample time")
        atypes = ['C', 'CV', 'AV']
        ahash = {'progression':2, 'convalescence':1, 'control':0}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRen2021CoV2scblk = getRen2021CoV2scblk


def getSaichi2021cov2(self, tn=1, tb=1):
    self.prepareData('COV386')
    atype = self.h.getSurvName("c disease severity")
    atypes = ['H', 'M', 'S']
    ahash = {'Severe':2, 'Moderate':1, 'Healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSaichi2021cov2 = getSaichi2021cov2


def getFernandez2021iAT2 (self, tn=1, tb=1):
    self.prepareData('COV387')
    atype = self.h.getSurvName("c timepoint")
    ahash = {'Day28':28, 'Day47':47, 'Day50':50, 'Day70':70}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c genotype")
    atypes = ['WT', 'Mut']
    ahash = {'DKC1 A386':0, 'DKC1 A386T':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == 70
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = ['Early', 'Late']
        ahash = {28:0, 70:1}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFernandez2021iAT2  = getFernandez2021iAT2 


def getRamirez2021drugsMm(self, tn=1, tb=1):
    self.prepareData('COV362')
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("([74]) .*", "\\1", str(k)) for k in atype]
    ahash = {'Control None 7':1, 'Control None 14':1, 'Bleomycin None 7':2,
            'Bleomycin None 14':3, 'Bleomycin Nintedanib 7':4,
            'Bleomycin Nintedanib 14':5}
    dval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c days after treatment")
    ahash = {'(-)':0, '(day 0-6)':1, '(day 7-13)':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c days after treatment")
    atypes = ['7', '14']
    ahash = {'7':0, '14':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = dval
        atypes = ['UnRx', '7dRx']
        ahash = {2:0, 4:1}
    if (tn == 3):
        atype = dval
        atypes = ['UnRx', '14dRx']
        ahash = {3:0, 5:1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRamirez2021drugsMm = getRamirez2021drugsMm


def getBorok2020at2mm(self, tn=1, tb=1):
    self.prepareData('COV389')
    atype = self.h.getSurvName('c genotype')
    atypes = ['WT', 'KO']
    ahash = {'Sftpc+/creERT2; Grp78f/f without Tmx':0,
            'Sftpc+/creERT2; Grp78f/f with Tmx':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBorok2020at2mm = getBorok2020at2mm


def getSimonis2021CoV2(self, tn=1, tb=1):
    self.prepareData('COV390')
    atype = self.h.getSurvName('c disease state')
    ahash = {'SARS-CoV-2 convalescent patient':1,
            'SARS-CoV-2 naïve individual':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c stimulation')
    atypes = ['NC', 'LPS', 'SP']
    ahash = {'stimulated with LPS':1, 'unstimulated':0,
            'stimulated with S-protein':2}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn >= 2):
        atypes = ['N', 'C']
        atype = tval
        ahash = {2:0, 1:1}
        atype = [atype[i] if aval[i] == (tn - 2)
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSimonis2021CoV2 = getSimonis2021CoV2


def getAng2021CoV2(self, tn=1, tb=1):
    self.prepareData('COV391')
    atype = self.h.getSurvName('c pcr status')
    ahash = {'positive':2, 'negative':1, 'healthy_control':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['H', 'CoV2']
    ahash = {'COVID-19':1, 'Healthy':0}
    if (tn == 2):
        atype = self.h.getSurvName('c pcr status')
        atypes = ['H', 'PCR-', 'PCR+']
        ahash = {'positive':2, 'negative':1, 'healthy_control':0}
    if (tn == 3):
        atype = self.h.getSurvName('c gene deletion status')
        atypes = ['H', 'ndel', 'del', 'NA']
        ahash = {'no_deletion_hash':1, 'deletion_star':2,
                'Not Applicable':3, 'healthy_control':0}
        atype = [atype[i] if tval[i] == 0 or tval[i] == 2
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAng2021CoV2 = getAng2021CoV2


def getRagan2021CoV2scblkHam(self, tn=1, tb=1):
    self.prepareData('COV392')
    atype = self.h.getSurvName('c treatment type')
    atypes = ['NV', 'V', 'V+C', 'V+O']
    ahash = {'vaccinated with SolaVAX and CpG1018 adjuvant':2,
            'no vaccination':0,
            'vaccinated with SolaVAX':1,
            'vaccinated with SolaVAX and ODN1668 adjuvant':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRagan2021CoV2scblkHam = getRagan2021CoV2scblkHam


def getShannon2020HBV(self, tn=1, tb=1):
    self.prepareData('COV394')
    atype = self.h.getSurvName('c visit')
    atypes = ['V3', 'V4', 'V5', 'V6', 'V7']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShannon2020HBV = getShannon2020HBV


def getLee2011Autoimmune(self, tn=1, tb=1):
    self.prepareData('COV400')
    atype = self.h.getSurvName('c disease')
    atypes = ['H', 'HC', 'RA', 'SLE', 'SOJIA', 'PTJIA']
    ahash = {'rheumatoid arthritis':2, 'healthy individual':0,
            'systemic lupus erythematosus':3, 'healthy child':1,
            'systemic-onset juvenile idiopathic arthritis':4,
            'polyarticular type juvenile idiopathic arthritis':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLee2011Autoimmune = getLee2011Autoimmune


def getTurnier2021SLE(self, tn=1, tb=1):
    self.prepareData('COV401')
    atype = self.h.getSurvName('c disease')
    atypes = ['N', 'SLE', 'JMnL', 'JML']
    ahash = {'cSLE skin lesion':1, 'JM Non-lesional skin':2,
            'Normal skin':0, 'JM Lesional skin':3}
    if (tn == 2):
        atypes = ['N', 'SLE']
        ahash = {'cSLE skin lesion':1, 'Normal skin':0}
    if (tn == 3):
        atypes = ['NL', 'L']
        ahash = {'JM Non-lesional skin':0, 'JM Lesional skin':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTurnier2021SLE = getTurnier2021SLE


def getPanwar2021SLE(self, tn=1, tb=1):
    self.prepareData('COV402')
    atype = self.h.getSurvName('c cell type')
    ahash = {'T cells':0, 'B cells':1, 'PMN':2, 'cDC':3, 'pDC':4, 'cMo':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state')
    atypes = ['H', 'SLE']
    ahash = {'healthy control':0, 'systemic lupus erythematosus (SLE)':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPanwar2021SLE = getPanwar2021SLE


def getBrooks2021SLE(self, tn=1, tb=1):
    self.prepareData('COV404')
    atype = self.h.getSurvName('c timepoint')
    ahash = {'1':1, '56':56, '84':84}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['P', 'T']
    ahash = {'Placebo':0, 'Tofacitinib':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBrooks2021SLE = getBrooks2021SLE


def getZhang2021SLE(self, tn=1, tb=1):
    self.prepareData('COV405')
    atype = self.h.getSurvName('c disease state')
    atypes = ['H', 'SLE']
    ahash = {'Healthy Control':0,
            'Patients of Systematic Lupus Erythematosus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2021SLE = getZhang2021SLE


def getJiang2021SLE(self, tn=1, tb=1):
    self.prepareData('COV406')
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['H', 'SLE']
    ahash = {'SLE patient':1, 'normal control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJiang2021SLE = getJiang2021SLE


def getBarnes2004jra(self, tn=1, tb=1):
    self.prepareData('COV407.1')
    atype = self.h.getSurvName('c Course Type')
    atypes = ['C', 'Pa', 'Po', 'JS']
    ahash = {'Ctrl':0, 'Pauci':1, 'Poly':2, 'JSpA':3}
    if (tn == 2):
        atypes = ['C', 'JRA']
        ahash = {'Ctrl':0, 'Pauci':1, 'Poly':1, 'JSpA':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBarnes2004jra = getBarnes2004jra


def getCharney2021MISC(self, tn=1, tb=1):
    self.prepareData('COV408')
    atype = self.h.getSurvName('c age range (years)')
    ahash = {'6-11':0, '18-23':1, '0-5':0, '12-17':0, '24-29':1, '35-40':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['H', 'CoV', 'M']
    ahash = {'covid':1, 'MISC':2, 'healthy':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCharney2021MISC = getCharney2021MISC


def getdeCevins2021MISC(self, tn=1, tb=1):
    self.prepareData('COV409')
    atype = self.h.getSurvName('c age group')
    ahash = {'pediatric':0, 'adult':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease group')
    atypes = ['H', 'K', 'M', 'MY']
    ahash = {'MISC (CoV2+)':2, 'MISC_MYO (CoV2+)':3, 'KD (CoV2+)':1, 'CTL':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['M', 'MYO+']
        ahash = {'MISC (CoV2+)':0, 'MISC_MYO (CoV2+)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getdeCevins2021MISC = getdeCevins2021MISC


def getdeCevins2021MISCscblk(self, tn=1, tb=1):
    self.prepareData('COV409.3')
    atype = self.h.getSurvName('c MajorCellTypes')
    ahash = {'Tcells':0, 'Bcells':1, 'Myeloid cells':2, 'HSC':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Group')
    atypes = ['H', 'Inf', 'CoV', 'K', 'M', 'MY']
    ahash = {'CTL':0, 'Acute-Inf(CoV2-)':1, 'Acute-Inf(CoV2+)':2,
            'MIS-C(CoV2+)':4, 'MIS-C_MYO(CoV2+)':5, 'KD(CoV2-)':3}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getdeCevins2021MISCscblk = getdeCevins2021MISCscblk


def getdeCevins2021MISCscblkII(self, tn=1, tb=1):
    self.prepareData('COV409.4')
    atype = self.h.getSurvName('c Group')
    atypes = ['H', 'Inf', 'CoV', 'K', 'M', 'MY']
    ahash = {'CTL':0, 'Acute-Inf(CoV2-)':1, 'Acute-Inf(CoV2+)':2,
            'MIS-C(CoV2+)':4, 'MIS-C_MYO(CoV2+)':5, 'KD(CoV2-)':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getdeCevins2021MISCscblkII = getdeCevins2021MISCscblkII


def getSchulert2020SJIA(self, tn=1, tb=1):
    self.prepareData('COV410')
    atype = self.h.getSurvName('c sample type')
    atypes = ['C', 'NOS', 'A', 'I', 'MAS']
    ahash = {'Control':0, 'NOS':1, 'Active':2, 'Inactive':3, 'NOS/MAS':4}
    if (tn == 2):
        atypes = ['C', 'D']
        ahash = {'Control':0, 'Active':1, 'NOS/MAS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSchulert2020SJIA = getSchulert2020SJIA


def getBrown2018SJIA(self, tn=1, tb=1):
    self.prepareData('COV411')
    atype = self.h.getSurvName('c disease state')
    atypes = ['C', 'A', 'I']
    ahash = {'Active SJIA':1, 'Inactive SJIA':2, 'Control patient':0}
    if (tn == 2):
        atypes = ['C', 'D']
        ahash = {'Active SJIA':1, 'Control patient':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBrown2018SJIA = getBrown2018SJIA


def getGorelik2013SJIA(self, tn=1, tb=1):
    self.prepareData('COV412')
    atype = self.h.getSurvName('c age_at_onset')
    atypes = ['<6', '>=6']
    ahash = {'LT6':0, 'GTE6':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGorelik2013SJIA = getGorelik2013SJIA


def getFall2007SJIA(self, tn=1, tb=1):
    self.prepareData('COV413')
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(".*:", "", str(k)) for k in atype]
    atypes = ['N', 'sJIA']
    ahash = {' normal control':0, ' new onset sJIA':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFall2007SJIA = getFall2007SJIA


def getGharib2016sarcoidosis(self, tn=1, tb=1):
    self.prepareData('COV414')
    atype = self.h.getSurvName('c individual')
    atypes = ['N', 'S']
    ahash = {'Sarcoidosis Patient':1, 'Normal Control Subject':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGharib2016sarcoidosis = getGharib2016sarcoidosis


def getAubert2012nomid(self, tn=1, tb=1):
    self.prepareData('COV415')
    atype = self.h.getSurvName('c disease development')
    atypes = ['N', 'L', 'Pre', 'Post']
    ahash = {'post-treatment non-lesional':3, 'lesional':1,
            'pre-treatment non-lesional':2, 'normal':0}
    if (tn == 2):
        atypes = ['C', 'L']
        ahash = { 'lesional':1, 'pre-treatment non-lesional':0, 'normal':0}
    if (tn == 3):
        atypes = ['C', 'L']
        ahash = { 'lesional':1, 'normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAubert2012nomid = getAubert2012nomid


def getAlmeida2011nomid(self, tn=1, tb=1):
    self.prepareData('COV416')
    atype = self.h.getSurvName('c src1')
    atypes = ['N', 'T']
    ahash = {'bone tumor':1, 'cartilage':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAlmeida2011nomid = getAlmeida2011nomid


def getCanna2014mas(self, tn=1, tb=1):
    self.prepareData('COV417')
    atype = self.h.getSurvName('c src1')
    atypes = ['H', 'M', 'Pre', 'Post']
    ahash = {'patient with NLRC4-MAS':1,
            'NOMID patients with active disease prior to anakinra treatment':2,
            'NOMID patients with inactive disease after anakinra treatment':3,
            'healthy pediatric controls':0}
    if (tn == 2):
        atypes = ['H', 'D']
        ahash = {'patient with NLRC4-MAS':1,
                'healthy pediatric controls':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCanna2014mas = getCanna2014mas


def getFrank2009jra(self, tn=1, tb=1):
    self.prepareData('COV418')
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" [0-9].*", "", str(k)) for k in atype]
    atypes = ['C', 'JIA', 'JDM']
    ahash = {'Neutrophil JIA':1, 'Neutrophil Control':0,
            'Neutrophil JDM':2, 'PBMC JIA':1,
            'PBMC JDM':2, 'Neutrophil control':0, 'PBMC Control':0}
    if (tn == 2):
        atypes = ['C', 'JIA', 'JDM']
        ahash = {'PBMC JIA':1, 'PBMC JDM':2, 'PBMC Control':0}
    if (tn == 3):
        atypes = ['C', 'JIA', 'JDM']
        ahash = {'Neutrophil JIA':1, 'Neutrophil Control':0,
                'Neutrophil JDM':2, 'Neutrophil control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFrank2009jra = getFrank2009jra


def getWong2016jia(self, tn=1, tb=1):
    self.prepareData('COV419')
    atype = self.h.getSurvName('c disease stage')
    atypes = ['C', 'U', 'T', 'R']
    ahash = {'healthy control':0,
            'JIA patient with clinical remission on medication':3,
            'JIA patient with active, untreated disease':1,
            'JIA patient with active disease on treatment':2}
    if (tn == 2):
        atypes = ['C', 'D']
        ahash = {'healthy control':0,
                'JIA patient with active, untreated disease':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWong2016jia = getWong2016jia


def getZhang2021CoV2pbmc(self, tn=1, tb=1):
    self.prepareData('COV420')
    atype = self.h.getSurvName('c immune response')
    atypes = ['U', 'A', 'S', 'R', 'P']
    ahash = {'uninfected':0, 'asymptomatic':1, 'symptomatic':2,
            'recovering':3, 're-detectable positive patients':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2021CoV2pbmc = getZhang2021CoV2pbmc


def getVanBooven2021GWI(self, tn=1, tb=1):
    self.prepareData('COV421')
    atype = self.h.getSurvName('c condition')
    ahash = {'Healthy Control':0, 'GWI':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c time point')
    atypes = ['T0', 'T1', 'T2']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c src1')
        atypes = ['HT0', 'HT1', 'HT2', 'T0', 'T1', 'T2']
        ahash = {'HCGWI_T1_PBMC':1, 'HCGWI_T0_PBMC':0, 'HCGWI_T2_PBMC':2,
                'GWI_T2_PBMC':5, 'GWI_T1_PBMC':4, 'GWI_T0_PBMC':3}
    if (tn == 4):
        atype = self.h.getSurvName('c src1')
        atypes = ['T2', 'T0']
        ahash = {'GWI_T2_PBMC':0, 'GWI_T0_PBMC':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanBooven2021GWI = getVanBooven2021GWI


def getMontoya2014pbmc(self, tn=1, tb=1):
    self.prepareData('MACV279')
    atype = self.h.getSurvName('c src1')
    atypes = ['M', '4', '10', '15']
    ahash = {'Media 6h':0, 'IL-10 24h':2, 'Media 24h':0, 'IL-10 6h':2,
            'IL-15 6h':3, 'IL-4 6h':1, 'IL-15 24h':3, 'IL-4 24h':1}
    if (tn == 2):
        atypes = ['C', 'IL15']
        ahash = {'Media 6h':0, 'Media 24h':0, 'IL-15 6h':1, 'IL-15 24h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMontoya2014pbmc = getMontoya2014pbmc


def getLey2020mm(self, tn=1, tb=0):
    self.prepareData('MACV280')
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("_.*", "", str(k)) for k in atype]
    atype = [re.sub("[0-9]+$", "", str(k)) for k in atype]
    ahash = {'unstim':0, '0.5hr':1, '1hr':2, '2hr':3, 'PD':4,
            'stim':5, 'p105AA':6}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c genotype')
    atypes = ['WT', 'Tpl2', 'Elk', 'Nfkb']
    ahash = {'WT':0, 'Map3k8[D270A]':1, 'Elk1-/-Elk4-/-':2,
            'Nfkb1[SS/AA]':3}
    if (tn == 2):
        atypes = ['WT', 'Tpl2']
        ahash = {'WT':0, 'Map3k8[D270A]':1}
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLey2020mm = getLey2020mm


def getRowley2014mm(self, tn=1, tb=0):
    self.prepareData('MACV281')
    atype = self.h.getSurvName('c treatment')
    ahash = {'unstimulated':0, 'stimulated with LPS':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c genotype/variation')
    atypes = ['WT', 'Tpl2']
    ahash = {'Tpl2 deficient':1, 'Tpl2 wild type':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRowley2014mm = getRowley2014mm


def getMacParland2018sam(self, tn=1, tb=0):
    self.prepareData('MACV286')
    atype = self.h.getSurvName('c CellType')
    atypes = ['NI', 'I']
    ahash = { 'Non-inflammatory_Macrophage':0, 'Inflammatory_Macrophage':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMacParland2018sam = getMacParland2018sam


def getRamachandran2019sam(self, tn=1, tb=0):
    self.prepareData('MACV288')
    atype = self.h.getSurvName('c src1')
    ahash = {'Liver':0, 'PBMC':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease status')
    atypes = ['H', 'F']
    ahash = {'fibrotic':1, 'healthy':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRamachandran2019sam = getRamachandran2019sam


def getDobie2019sam(self, tn=1, tb=0):
    self.prepareData('MACV290')
    atype = self.h.getSurvName('c disease status')
    atypes = ['H', '6', '72']
    ahash = {'CCl4 6wk':1, 'Healthy':0, 'CCl4 72h':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDobie2019sam = getDobie2019sam


def getDobie2019samII(self, tn=1, tb=0):
    self.prepareData('MACV291')
    atype = self.h.getSurvName('c Sample')
    atype = [re.sub(" [0-9]$", "", str(k)) for k in atype]
    atypes = ['C', 'U', 'BDL', 'CCl']
    ahash = {'Control':0, 'Peak BDL':2, 'Peak CCl4':3, 'Uninjured':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDobie2019samII = getDobie2019samII


def getJaitin2019lam (self, tn=1, tb=0):
    self.prepareData('MACV283')
    atype = self.h.getSurvName('c replicate id')
    atype = [re.sub(" _[0-9].*", "", str(k)) for k in atype]
    atypes = ['C', 'O']
    ahash = {'OAT_CD45_obese':1, 'OAT_CD45_CTL':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJaitin2019lam  = getJaitin2019lam 


def getKerenShaul2017dam(self, tn=1, tb=0):
    self.prepareData('MACV282')
    atype = self.h.getSurvName('c organ')
    ahash = {'Spinal cord':0, 'Brain':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c selection marker')
    ahash = {'CD45+':0, 'CD45+CD11b+Gr1-':0,
            'CD45-intCD11b-intCD11c+Gr1-':2,
            'CD45-intCD11b-intGr1-':2,
            'CD45-intCD11b-intCD11c+Gr1-;  CD45-intCD11b-intGr1-':2,
            'CD45-intCD11bintCD11c+Gr1-':2,
            'CD45-intCD11bintGr1-':2}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c treatment')
    atypes = ['Age', 'AD', 'ALS']
    ahash = {'ALS':2, "Alzheimer's disease":1, 'Aging':0}
    if (tn == 2):
        atype = [atype[i] if sval[i] == tb
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKerenShaul2017dam = getKerenShaul2017dam


def getHe2021(self, tn=1, ta=0, tb=0):
    self.prepareData('MACV293')
    atype = self.h.getSurvName('c treatment')
    ahash = {'TRAM':1, 'DMSO':0, 'PANO':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c cell line')
    ahash = {'THP-1':1, '':0}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c polarization')
    atypes = ['C', 'IFNLPS', 'IL4']
    ahash = {'IL4':2, 'IFNLPS':1, 'control':0, 'Control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta and sval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = ['DMSO', 'TRAM', 'PANO']
        ahash = {0:0, 1:1, 2:2}
        atype = [atype[i] if aval[i] == ta and sval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'IL4', 'tIL4', 'pIL4']
        ahash = {'0-0-0':0, '2-0-0':1, '2-0-1':2, '2-0-2':3}
    if (tn == 5):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'LPS', 'tLPS', 'pLPS']
        ahash = {'0-0-0':0, '1-0-0':1, '1-0-1':2, '1-0-2':3}
    if (tn == 6):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'IL4', 'tIL4', 'pIL4']
        ahash = {'0-1-0':0, '2-1-0':1, '2-1-1':2, '2-1-2':3}
    if (tn == 7):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'LPS', 'tLPS', 'pLPS']
        ahash = {'0-1-0':0, '1-1-0':1, '1-1-1':2, '1-1-2':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHe2021 = getHe2021


def getPG2021oxLDL(self, tn=1, ta=0, tb=0):
    self.prepareData('PG21', "/Users/dtv004/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c Sample ID')
    atypes = ['WT0', 'WT24', 'GIV0', 'GIV24']
    ahash = { 'TGPM_WT_0h_S16':0, 'TGPM_GIVKO_24h_S21':3,
            'TGPM_WT_24h_S7':1, 'TGPM_GIVKO_0h_S3':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPG2021oxLDL = getPG2021oxLDL


def getArijs2009uc(self, tn=1):
    self.prepareData("PLP142")
    atype = self.h.getSurvName("c WK8RSPHM")
    atypes = ["R", "NR"]
    ahash = {"Yes": 0, "No": 1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArijs2009uc = getArijs2009uc


def getBohne2014hcv(self, tn=1):
    self.prepareData("LIV69")
    atype = self.h.getSurvName("c src1")
    atypes = ["H", "Tb", "T12", 'NTb', 'NTr']
    ahash = {'tolerant_12 months after IS discontinuation':2,
            'non-tolerant_before immunosuppression (IS)':3,
            'non-tolerant_at the time of rejection':4,
            'tolerant_before immunosuppression (IS)':1,
            'healthy living liver donor':0}
    if (tn == 2):
        atype = self.h.getSurvName('c patient status')
        ahash = {'stable HCV-infected liver recipients':1,
                'no stable HCV-infected liver recipients':2,
                'non-transplanted controls':0}
        atypes = ['H', 'S', 'nS']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBohne2014hcv = getBohne2014hcv


def getNingappa2022(self, tn=1):
    self.prepareData("MACV295")
    atype = self.h.getSurvName("c groups")
    atypes = ['late_NR', 'pre_NR', 'early_NR', 'pre_R', 'late_R', 'early_R']
    ahash = {}
    if (tn == 2):
        atypes = ['pre_R', 'pre_NR']
    if (tn == 3):
        atypes = ['early_R', 'early_NR']
    if (tn == 4):
        atypes = ['late_R', 'late_NR']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNingappa2022 = getNingappa2022


def getChruscinski2022(self, tn=1):
    self.prepareData("LIV70")
    atype = self.h.getSurvName("c group")
    atypes = ['T', 'NT']
    ahash = {'Operatiionally Tolerant':0, 'Non-Tolerant':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChruscinski2022 = getChruscinski2022


def getGoossens2016(self, tn=1):
    self.prepareData("LIV71")
    atype = self.h.getSurvName("c 32 gene prediction")
    atypes = ['Good', 'Poor']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGoossens2016 = getGoossens2016


def getTrepo2018III(self, tn=1, tb=0):
    self.prepareData("LIV72")
    atype = self.h.getSurvName('c disease status')
    atypes = ['AH', 'AS', 'AC']
    ahash = {'Alcoholic cirrhosis':2,
            'Mild acute alcoholic hepatitis':0,
            'Alcoholic steatosis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTrepo2018III = getTrepo2018III


def getPinyol2021Liver(self, tn=1, ta=0, tb=0):
    self.prepareData("LIV53")
    atype = self.h.getSurvName('c tissue')
    atypes = ['H', 'A', 'N', 'C', 'T']
    ahash = {'NASH-HCC tumor':4, 'NASH liver':2,
             'Non-tumoral NASH liver adjacent to HCC':1,
             'Cirrhotic liver':3,'Healthy liver':0}
    if (tn == 2):
        atypes = ['H', 'T']
        ahash = {'NASH-HCC tumor':1, 'Healthy liver':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPinyol2021Liver = getPinyol2021Liver


def getLim2013HCC(self, tn=1, ta=0, tb=0):
    self.prepareData("LIV45")
    atype = self.h.getSurvName('c ajcc stage')
    ahash = {'1':1, '2':2, '3a':3, '3b':3, '3c':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c tissue')
    atypes = ['NT', 'T']
    ahash = {'liver tumor':1, 'adjacent non-tumor liver':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLim2013HCC = getLim2013HCC


def getVillanueva2015CirHCC(self, tn=1, ta=0, tb=0):
    self.prepareData("LIV41")
    atype = self.h.getSurvName('c gender')
    ahash = {'female':0, 'male':1, 'NA':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c src1')
    atypes = ['Cir', 'HCC']
    ahash = {'hepatocellular carcinoma':1, 'cirrhosis':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVillanueva2015CirHCC = getVillanueva2015CirHCC


def getHonda2012(self, tn=1, ta=0):
    self.prepareData("LIV73")
    atype = self.h.getSurvName("c time")
    ahash = {'control':0, '8 weeks':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c clincal outcome")
    atypes = ['NR', 'R', 'NA']
    ahash = {'':2, 'Recurrence of HCC':1, 'nonRecurrence of HCC':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c patient")
        atypes = ['NR', 'R']
        ahash = {}
        for i in self.h.aRange():
            if aval[i] == 0 or aval[i] == 1:
                ahash[atype[i]] = aval[i]
        atype = [atype[i] if tval[i] == 0
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHonda2012 = getHonda2012

def getGlobal(self, tn=1):
    self.prepareData("GL1")
    atype = self.h.getSurvName('c TissueState')
    atypes = ['nan', 'solid', 'cellline', 'liquid', 'bonemarrow',
            'mix', 'unknown', 'lymphoma', 'lymph',
            'cellline_p3', 'saliva']
    ahash = {'':0}
    if tn == 2:
        atypes = ['solid', 'liquid', 'cellline']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobal = getGlobal

def getGlobalMouse(self, tn=1):
    self.prepareData("GL2")
    atype = self.h.getSurvName('c TissueState')
    atypes = ['solid', 'nan', 'cellline', 'bonemarrow', 'lymph',
            'liquid', 'solid&cellline', 'mix', 'unknown', 'lymphoma']
    ahash = {'':1}
    if tn == 2:
        atypes = ['solid', 'liquid', 'cellline']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobalMouse = getGlobalMouse

def getGlobalRat(self, tn=1):
    self.prepareData("GL46")
    atype = self.h.getSurvName('c TissueState')
    atypes = ['nan', 'solid', 'cellline', 'liquid', 'lymph', 'bonemarrow']
    ahash = {'':0}
    if tn == 2:
        atypes = ['solid', 'liquid', 'cellline']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobalRat = getGlobalRat

def getGlobalDog(self, tn=1):
    self.prepareData("MACV336.3")
    atype = self.h.getSurvName('c tissue')
    atypes = ['solid', 'liquid']
    ahash = {'Blood':1}
    for k in atype[2:]:
        if k not in ahash:
            ahash[k] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobalDog = getGlobalDog

def getGlobalMonkey(self, tn=1):
    self.prepareData("MACV336.4")
    atype = self.h.getSurvName('c tissue')
    atypes = ['solid']
    ahash = {}
    for k in atype[2:]:
        if k not in ahash:
            ahash[k] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobalMonkey = getGlobalMonkey

def getGlobalBaboon(self, tn=1):
    self.prepareData("GL21")
    atype = self.h.getSurvName('c tissue')
    atypes = ['solid']
    ahash = {}
    for k in atype[2:]:
        if k not in ahash:
            ahash[k] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobalBaboon = getGlobalBaboon

def getGlobalNeu(self, tn=1):
    self.prepareData("GL1")
    cfile = "/booleanfs2/sahoo/Data/BooleanLab/RxCovea/human-gpl570-tissue.txt"
    df = bone.pd.read_csv(cfile, sep="\t", index_col=0)
    atype = ["", ""] + [str(df['c Tissue'][k]) if k in df['c Tissue']
                        else "" for k in self.h.headers[2:]]
    atypes = ["Other", 'Neu']
    ahash = {}
    for title in atype[2:]:
        if "neutrophil" in title.lower():
            ahash[title] = 1
        else:
            ahash[title] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGlobalNeu = getGlobalNeu

def getRojo2019(self, tn=1):
    self.prepareData("MACV333")
    atype = self.h.getSurvName('c genotype')
    atypes = ['WT', 'KO']
    ahash = {'FIRE -/-':1, 'FIRE +/+':0}
    if tn == 2:
        tissue = self.h.getSurvName('c tissue')
        bhash = {"Peyer's Patches/Ileum":0, 'hippocampus':1}
        ttype = [bhash[k] if k in bhash else None for k in tissue]
        atype = [atype[i] if ttype[i] == 1 else None for i in range(len(atype))]
    if tn == 3:
        tissue = self.h.getSurvName('c tissue')
        bhash = {"Peyer's Patches/Ileum":0, 'hippocampus':1}
        ttype = [bhash[k] if k in bhash else None for k in tissue]
        atype = [atype[i] if ttype[i] == 0 else None for i in range(len(atype))]
    if tn == 4:
        tissue = self.h.getSurvName('c tissue')
        btypes = ['PP', 'HIP']
        bhash = {"Peyer's Patches/Ileum":0, 'hippocampus':1}
        atype = [atypes[ahash[atype[i]]]+"-"+btypes[bhash[tissue[i]]]
                 if tissue[i] in bhash and atype[i] in ahash
                 else None for i in range(len(atype))]
        atypes = ['WT-HIP', 'KO-HIP', 'WT-PP', 'KO-PP']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRojo2019 = getRojo2019

def getGenentechCellLines(self, tn=1):
    self.prepareData("GL9")
    atype = self.h.getSurvName("c Cell_line")
    atypes = ['HCT 116', 'SW 480', 'DLD-1', 'Caco-2']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGenentechCellLines = getGenentechCellLines

def getXue2014(self, mt=1):
    self.prepareData("MAC2")
    atype = self.h.getSurvName("c Type")
    atypes = ['M0', 'M1', 'M2']
    ahash = {\
            'M_GMCSF_IL4_72h':2,
            'M_GMCSF_IFNg_72h':1,
            'M0_GMCSF_72h':0}
    if mt == 2:
        ahash = {\
                'M_GMCSF_IL4_24h':2,
                'M_GMCSF_IFNg_24h':1,
                'M0_GMCSF_24h':0}
    if mt == 3:
        ahash = {\
                'M_GMCSF_IL4_12h':2,
                'M_GMCSF_IFNg_12h':1,
                'M0_GMCSF_12h':0}
    if mt == 4:
        ahash = {\
                'M0_GMCSF_0h':0,
                'M0_GMCSF_12h':0,
                'M0_GMCSF_24h':0,
                'M0_GMCSF_48h':0,
                'M0_GMCSF_6h':0,
                'M0_GMCSF_72h':0,
                'M0_MCSF_0h':0,
                'M1/2_GMCSF_24h':0,
                'M_GMCSF_IFNg_30min':0,
                'M_GMCSF_IFNg_1h':0,
                'M_GMCSF_IFNg_2h':0,
                'M_GMCSF_IFNg_4h':1,
                'M_GMCSF_IFNg_6h':1,
                'M_GMCSF_IFNg_12h':1,
                'M_GMCSF_IFNg_24h':1,
                'M_GMCSF_IFNg_72h':1,
                'M_GMCSF_IL4_30min':2,
                'M_GMCSF_IL4_1h':2,
                'M_GMCSF_IL4_2h':2,
                'M_GMCSF_IL4_4h':2,
                'M_GMCSF_IL4_6h':2,
                'M_GMCSF_IL4_12h':2,
                'M_GMCSF_IL4_24h':2,
                'M_GMCSF_IL4_72h':2,
                'M_MCSF_IL4_72h':2}
    if mt == 5:
        atype = [re.sub("M[01_].*", "M", str(k)) for k in atype]
        ahash = {'Monocyte_CD14+':0, 'M':1, 'B cell':2, 'T_CD3+':3,
                 'DC_upLPS_10_24h':5, 'Tconv':3, 'DC_reg':5, 'iTreg_IL-2':3,
                 'DC_imm':5, 'Tstim':3,  'iTreg':3, 'Tresting':3, 'DC_mat':5,
                 'Tnaive':3, 'NK_cells':4, 'Treg':3, 'Tmemory':3}
        atypes = ['Mono', 'Mac', 'B', 'T', 'NK', 'DC']
    if mt == 6:
        atype = [re.sub("M[01_].*", "M", str(k)) for k in atype]
        ahash = {'DC_upLPS_10_24h':1, 'DC_reg':0,
                 'DC_imm':2, 'DC_mat':3}
        atypes = ['DC', 'LPS', 'iDC', 'maDC']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getXue2014 = getXue2014

def getMorell2018 (self, tn=1):
    self.prepareData("MAC8")
    atype = self.h.getSurvName('c cell type')
    atypes = ['Mo', 'Mac']
    ahash = {'peripheral blood monocyte':0, 'alveolar macrophage':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMorell2018 = getMorell2018
def getMisharin2017Mm (self, tn=1):
    self.prepareData("MAC9")
    atype = self.h.getSurvName('c cell type')
    atypes = ['Mo', 'Mac', 'iMac']
    ahash = {'Alveolar Macrophages':1, 'Lung Monocytes':0,
             'Interstitial Macrophages':2}
    if (tn == 2):
        atypes = ['Mo', 'Mac']
        ahash = {'Alveolar Macrophages':1, 'Lung Monocytes':0}
    if (tn == 3):
        atypes = ['Mo', 'iMac']
        ahash = {'Interstitial Macrophages':1, 'Lung Monocytes':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMisharin2017Mm = getMisharin2017Mm
def getMartinez2006 (self, tn=1):
    self.prepareData("MAC19.1")
    atype = self.h.getSurvName('c Type')
    atypes = ['Mo', 'Mac']
    ahash = {'Monocyte at 3 days':1,
             'classical or M1 activated macrophages':1,
             'Macrophage at 7 days':1,
             'Alternative or M2 activated macrophages':1,
             'Monocyte at T0':0, 'Monocyteat 3 days':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinez2006 = getMartinez2006
def getMaouche2008 (self, tn=1):
    self.prepareData("MAC21")
    atype = self.h.getSurvName('c Cell type')
    atypes = ['Mo', 'Mac']
    ahash = {'monocyte':0, 'macrophage':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaouche2008 = getMaouche2008
def getCassetta2019 (self, tn=1):
    self.prepareData("MAC29")
    atype = self.h.getSurvName('c cell type')
    atypes = ['Mo', 'TAM']
    ahash = {'monocytes':0, 'tumour associated macrophages':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCassetta2019 = getCassetta2019
def getWoetzel2014(self, tn=1):
    self.prepareData("MAC31")
    atype = self.h.getSurvName("c disease state")
    atype2 = self.h.getSurvName("c clinical status")
    atype = [atype[i] + atype2[i] for i in range(len(atype))]
    atypes = ['HC', 'RA', 'OA']
    ahash = {'healthy control':0,
            'rheumatoid arthritis':1,
            'synovial tissue isolated from osteoarthritic joint':2,
            'osteoarthritis':2,
            'normal control':0}
    if (tn == 2):
        atypes = ['HC', 'OA']
        ahash = {'healthy control':0,
                'synovial tissue isolated from osteoarthritic joint':1,
                'osteoarthritis':1, 'normal control':0}
    if (tn == 3):
        atypes = ['HC', 'RA']
        ahash = {'healthy control':0, 'normal control':0,
                'rheumatoid arthritis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWoetzel2014 = getWoetzel2014
def getSutherland2011 (self, tn=1):
    self.prepareData("MACV128")
    atype = self.h.getSurvName('c health status')
    atypes = ['H', 'S', 'PS']
    ahash = {'HEALTHY':0, 'POST_SURGICAL':2, 'SEPSIS':1}
    if (tn == 2):
        atypes = ['H', 'S']
        ahash = {'HEALTHY':0, 'SEPSIS':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSutherland2011 = getSutherland2011
def getParnell2013(self, tn=1, tb=0):
    self.prepareData("MACV272")
    atype = self.h.getSurvName('c disease status')
    atypes = ['H', 'S', 'NS']
    ahash = {'healthy':0, 'sepsis survivor':1, 'sepsis nonsurvivor':2}
    if (tn == 2):
        atypes = ['H', 'S']
        ahash = {'healthy':0, 'sepsis survivor':1, 'sepsis nonsurvivor':1}
    if (tn == 3):
        atypes = ['S', 'NS']
        ahash = {'sepsis survivor':0, 'sepsis nonsurvivor':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getParnell2013 = getParnell2013
def getParnell2012(self, tn=1, tb=0):
    self.prepareData("MACV129")
    atype = self.h.getSurvName('c sample type')
    atypes = ['HC', 'B', 'V', 'B+V', 'I']
    ahash = { 'influenza A pneumonia':2, 'healthy control':0,
             'bacterial pneumonia':1,
             'systemic inflammatory response':4,
             'mixed bacterial and influenza A pneumonia':3}
    if (tn == 2):
        atypes = ['HC', 'B']
        ahash = {'healthy control':0, 'bacterial pneumonia':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getParnell2012 = getParnell2012
def getSweeney2015(self, tn=1):
    self.prepareData("MACV130")
    atype = self.h.getSurvName("c disease")
    atypes = ['C', 'S', 'SS', 'SIRS']
    ahash = {'SepticShock':2, 'SIRS':3, 'Sepsis':1, 'Control':0}
    if (tn == 2):
        atypes = ['C', 'SIRS']
        ahash = {'SIRS':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSweeney2015 = getSweeney2015
def getMcHugh2015(self, tn=1):
    self.prepareData("MACV131.2")
    atype = self.h.getSurvName("c group")
    atypes = ['C', 'S']
    ahash = {'post-surgical':0, 'Sepsis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMcHugh2015 = getMcHugh2015
def getPresnell2015(self, tn=1):
    self.prepareData("MACV132")
    atype = self.h.getSurvName("c study group")
    atypes = ['H', 'Sm', 'So', 'T2D']
    ahash = {'Uninfected type 2 diabetes mellitus':3,
             'Other sepsis':2, 'Septicemic melioidosis':1,
             'Uninfected healthy':0}
    if (tn == 2):
        atypes = ['H', 'S']
        ahash = {'Other sepsis':1, 'Septicemic melioidosis':1,
                 'Uninfected healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPresnell2015 = getPresnell2015
def getSuarez2015(self, tn=1):
    self.prepareData("MACV133")
    atype = self.h.getSurvName("c condition")
    atypes = ['C', 'S']
    ahash = {'COINFECTION':1, 'Healthy Control':0, 'VIRUS':1, 'BACTERIA':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'B', 'V', 'CO']
        ahash = {'COINFECTION':3, 'Healthy Control':0, 'VIRUS':2, 'BACTERIA':1}
    if (tn == 3):
        atype = self.h.getSurvName("c race")
        atypes = ['W', 'B', 'A']
        ahash = {'White':0, 'Black':1, 'Asian':2}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSuarez2015 = getSuarez2015
def getEmonts2008(self, tn=1, tb='Blood'):
    self.prepareData("MACV134")
    atype = self.h.getSurvName("c src1")
    atype = [re.sub("-", "", str(k)) for k in atype]
    disease = [re.sub(".\/.*", "", str(k)) for k in atype]
    atype = [re.sub("^[^/]*.\/", "", str(k)) for k in atype]
    celltype = [str(k).split("/")[0] if len(str(k).split("/")) > 0
             else None for k in atype]
    ahash = {'Monocyte', 'Lymphocyte', 'Blood'}
    timepoint = [str(k).split("/")[1] if len(str(k).split("/")) > 1
             else '' for k in atype]
    ahash = {'', 'T24', 'T8', 'T0', 'T72'}
    atype = disease
    atypes = ['H', 'S']
    ahash = {'Control':0, 'Patient':1}
    if (tn == 2):
        atype = [atype[i] if celltype[i] == tb
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i]+"-"+timepoint[i] if celltype[i] == tb
                 else None for i in range(len(atype))]
        ahash = {'Control-':0, 'Patient-T24':3, 'Patient-T8':2,
                 'Patient-T0':1, 'Patient-T72':4}
        atypes = ['C', 'T0', 'T8', 'T24', 'T72']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEmonts2008 = getEmonts2008
def getPankla2009I(self, tn=1):
    self.prepareData("MACV135")
    atype = self.h.getSurvName("c Illness")
    atypes = ['H', 'Sm', 'So', 'Sr', 'T2D']
    ahash = {'Sepsis/Other infections':2, 'Control/Recovery':3,
             'Sepsis/Melioidosis':1, 'Control/type 2 diabetes':4,
             'Control/healthy':0}
    if (tn == 2):
        atypes = ['H', 'Sm']
        ahash = {'Sepsis/Melioidosis':1,'Control/healthy':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPankla2009I = getPankla2009I
def getPankla2009II(self, tn=1):
    self.prepareData("MACV136")
    atype = self.h.getSurvName("c Illness")
    atypes = ['H', 'Sm', 'So', 'T2D']
    ahash = {'Control/healthy':0, 'Sepsis/Melioidosis':1,
             'Sepsis/Other infections':2, 'Control/type 2 diabetes':3}
    if (tn == 2):
        atypes = ['H', 'Sm']
        ahash = {'Control/healthy':0, 'Sepsis/Melioidosis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPankla2009II = getPankla2009II
def getZaas2009(self, tn=1):
    self.prepareData("MACV137")
    atype = self.h.getSurvName('c symptom group')
    atypes = ['A', 'S']
    ahash = {'Symptomatic':1, 'Asymptomatic':0}
    if (tn == 2):
        atype = self.h.getSurvName('c timepoint')
        atypes = ['B', 'P']
        ahash = {'baseline':0, 'T - peak symptoms':1}
    if (tn == 3):
        sym = self.h.getSurvName('c symptom group')
        time = self.h.getSurvName('c timepoint')
        atype = [str(sym[i])+"-"+str(time[i]) for i in range(len(sym))]
        ahash = {'Symptomatic-baseline':1, 'Symptomatic-T - peak symptoms':3,
                 'Asymptomatic-baseline':0, 'Asymptomatic-T - peak symptoms':2}
        atypes = ['AB', 'SB', 'AP', 'SP']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZaas2009 = getZaas2009
def getParnell2011 (self, tn=1):
    self.prepareData("MACV138")
    atype = self.h.getSurvName("c src1")
    group = [str(k).split("_")[0] if len(str(k).split("_")) > 0
             else '' for k in atype]
    timepoint = [str(k).split("_")[1] if len(str(k).split("_")) > 1
             else '' for k in atype]
    atype = group
    atypes = ['V', 'B', 'I']
    ahash = {'Severe Influenza':2, 'bacterial pneumonia':1, 'Vaccine':0}
    if (tn == 2):
        atype = [str(group[i])+"-"+str(timepoint[i]) for i in range(len(atype))]
        ahash = {'Vaccine-baseline':0,'Vaccine-day 7':1,
                 'bacterial pneumonia-day 1':2, 'bacterial pneumonia-day 2':3,
                 'bacterial pneumonia-day 3':4, 'bacterial pneumonia-day 4':5,
                 'bacterial pneumonia-day 5':6,
                 'Severe Influenza-day 1':7, 'Severe Influenza-day 2':8, 
                 'Severe Influenza-day 3':9, 'Severe Influenza-day 4':10,
                 'Severe Influenza-day 5':11}
        atypes = ['V0', 'V7', 'B1', 'B2', 'B3', 'B4', 'B5', 'I1', 'I2', 'I3', 'I4', 'I5']
    if (tn == 3):
        atype = [str(group[i])+"-"+str(timepoint[i]) for i in range(len(atype))]
        ahash = {'Vaccine-baseline':0, 'bacterial pneumonia-day 2':1,
                 'Severe Influenza-day 1':2}
        atypes = ['V0', 'B2', 'I1']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getParnell2011 = getParnell2011
def getSmith2014I(self, tn=1):
    self.prepareData("MACV139")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("...$", "", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'Con':0, 'Inf':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2014I = getSmith2014I
def getSmith2014II(self, tn=1):
    self.prepareData("MACV140")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("...$", "", str(k)) for k in atype]
    atypes = ['C', 'I']
    ahash = {'Con':0, 'Inf':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2014II = getSmith2014II
def getSmith2014III(self, tn=1):
    self.prepareData("MACV141")
    atype = self.h.getSurvName("c Title")
    atype = [re.sub("...$", "", str(k)) for k in atype]
    atypes = ['Con', 'Inf', 'Vir', 'NEC']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2014III = getSmith2014III
def getAhn2015(self, tn=1):
    self.prepareData("MACV142")
    atype = self.h.getSurvName("c pathogen")
    atypes = ['C', 'S']
    ahash = {'-':0,
            'Staphylococcus aureus':1,
            'Escherichia coli':1,
            'Staphylococcus aureus and Streptococcus pneumoniae':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'Ec', 'Sa', 'S']
        ahash = {'-':0,
                'Staphylococcus aureus':2,
                'Escherichia coli':1,
                'Staphylococcus aureus and Streptococcus pneumoniae':3}
    if (tn == 3):
        atype = self.h.getSurvName("c ethnicity")
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'White':0, 'Unknown':3, 'Black':1, 'unknown':3, 'Asian':2,
                'black':1, 'white':0}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAhn2015 = getAhn2015
def getAhn2015Mm(self, tn=1):
    self.prepareData("MACV143")
    atype = self.h.getSurvName('c infection status')
    atypes = ['C', 'I']
    ahash = {'infected':1, '-':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAhn2015Mm = getAhn2015Mm
def getMejias2013I(self, tn=1):
    self.prepareData("MACV144")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'RSV', 'RSVf', 'HRV', 'Inf']
    ahash = {'Dallas  Acute RSV':1, 'Dallas healthy control':0,
             'Finnish Acute RSV':1, 'Dalls RSV 1-2 months follow up':2,
             'Dallas Acute HRV':3, 'Dallas Influenza A':4,
             'Finnish healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMejias2013I = getMejias2013I
def getMejias2013II(self, tn=1):
    self.prepareData("MACV145")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'RSV']
    ahash = {'Columbus Acute RSV':1, 'Columbus healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMejias2013II = getMejias2013II
def getHu2013(self, tn=1):
    self.prepareData("MACV146")
    atype = self.h.getSurvName('c group')
    atypes = ['Ctrl', 'RNA-virus','DNA-virus']
    ahash = {'Adenovirus':2, 'Virus-negative Control':0,'HHV6':2,
            'Enterovirus':1,'Rhinovirus':1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atypes = ['HC', 'I']
        ahash = {'Adenovirus':1, 'Virus-negative Control':0,'HHV6':1,
                'Enterovirus':1,'Rhinovirus':1}
    if (tn == 3):
        atype = self.h.getSurvName("c ethnicity")
        atypes = ['W', 'B', 'A', 'O']
        ahash = {'White':0, 'Black':1, 'Other':3}
        atype = [atype[i] if aval[i] == 0
                else None for i in range(len(atype))]
    if (tn == 4):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'Rhinovirus':1}
    if (tn == 5):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'Enterovirus':1}
    if (tn == 6):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'Adenovirus':1}
    if (tn == 7):
        atypes = ['HC', 'I']
        ahash = {'Virus-negative Control':0, 'HHV6':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHu2013 = getHu2013
def getHerberg2013(self, tn=1):
    self.prepareData("MACV147")
    atype = self.h.getSurvName('c infecting pathogen')
    atypes = ['C', 'B', 'RSV', 'Inf']
    ahash = {'gram positive bacterial infection':1, 'RSV':2, 
             'Influenza A H1N1/09':3, 'none':0}
    if (tn == 2):
        atypes = ['C', 'B']
        ahash = {'gram positive bacterial infection':1,'none':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHerberg2013 = getHerberg2013
def getKwissa2014(self, tn=1):
    self.prepareData("MACV148")
    atype = self.h.getSurvName('c status')
    atypes = ['C', 'DF', 'DHF', 'CV']
    ahash = {'control':0, 'convalescent':3, 'DHF':2, 'DF':1}
    if (tn == 2):
        atypes = ['C', 'DF']
        ahash = {'control':0, 'DF':1}
    if (tn == 3):
        atypes = ['C', 'DHF']
        ahash = {'control':0, 'DHF':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKwissa2014 = getKwissa2014
def getTabone2018(self, tn=1):
    self.prepareData("MACV149")
    atype = self.h.getSurvName('c sapsii')
    atypes = ['C', 'SH', 'SL']
    ahash = {'SAPSII-High':1, 'NA':0, 'SAPSII-Low':2}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'SAPSII-High':1, 'NA':0, 'SAPSII-Low':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTabone2018 = getTabone2018
def getvandeWeg2015(self, tn=1):
    self.prepareData("MACV152")
    atype = self.h.getSurvName('c Characteristics[disease]')
    atypes = ['C', 'DHF']
    ahash = {'Dengue Hemorrhagic Fever':1, 'Healthy Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getvandeWeg2015 = getvandeWeg2015
def getTsalik2016(self, tn=1):
    self.prepareData("MACV154")
    atype = self.h.getSurvName('c infection_status')
    atypes = ['C', 'B', 'V']
    ahash = {'non-infectious illness':0, 'viral':2, 'bacterial':1}
    if (tn == 2):
        atypes = ['C', 'B']
        ahash = {'non-infectious illness':0, 'bacterial':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTsalik2016 = getTsalik2016
def getLill2013(self, tn=1):
    self.prepareData("MACV156")
    atype = self.h.getSurvName('c sample group')
    atypes = ['C', 'B']
    ahash = {'healthy control':0, 'bacterial meningitis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLill2013 = getLill2013
def getPayen2009(self, tn=1):
    self.prepareData("MACV179")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" *[VT].*$", "", str(k)) for k in atype]
    atypes = ['C', 'S', 'GN', 'GP', 'M']
    ahash = {'Sepsis gram-negative':2, 'Sepsis gram-positive':3,
             'Sepsis':1, 'Sepsis mixed infection':4, 'Control':0}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'Sepsis':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPayen2009 = getPayen2009
def getShalova2015(self, tn=1):
    self.prepareData("MACV196")
    atype = self.h.getSurvName('c src1')
    atypes = ['C', 'S', 'SR']
    ahash = {'healthy donor':0, 'patient recovering from sepsis':2,
             'patient with sepsis':1}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'healthy donor':0,'patient with sepsis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShalova2015 = getShalova2015
def getMartinezPaz2020(self, tn=1):
    self.prepareData("MACV255")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['C', 'S', 'NS']
    ahash = {'septic shock':1, 'control patient':0, 'non-septic shock':2}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'septic shock':1, 'control patient':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMartinezPaz2020 = getMartinezPaz2020
def getBandyopadhyay2020(self, tn=1): #URINE sample
    self.prepareData("MACV257")
    tissue = self.h.getSurvName('c tissue')
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['C', 'S']
    ahash = {'SEPTIC':1, 'VASCULAR':0}
    if (tn == 2):
        atype = [atype[i] if tissue[i] == 'URINE'
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tissue[i] == ''
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBandyopadhyay2020 = getBandyopadhyay2020
def getNg2020(self, tn=1):
    self.prepareData("MACV258")
    atype = self.h.getSurvName('c clinical group')
    atypes = ['C', 'CS', 'PS']
    ahash = {'Confirmed late-onset sepsis':1, 'No late-onset sepsis':0,
             'Possible late-onset sepsis':2}
    if (tn == 2):
        atypes = ['C', 'S']
        ahash = {'Confirmed late-onset sepsis':1, 'No late-onset sepsis':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNg2020 = getNg2020
def getNovakovic2023(self, tn=1):
    self.prepareData("MACV346")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(".$", "", str(k)) for k in atype]
    atypes = ['M', 'NM', '3', '6']
    ahash = {'MMA 3 days ':2, 'monocyte ':0, 'Naive macrophage ':1, 'MMA 6 days ':3}
    if (tn == 2):
        atypes = ['Mo', 'Mac']
        ahash = {'MMA 3 days ':1, 'monocyte ':0, 'Naive macrophage ':1, 'MMA 6 days ':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNovakovic2023 = getNovakovic2023
def getPannaraj2022(self, tn=1):
    self.prepareData("MACV303")
    atype = self.h.getSurvName('c src1')
    atypes = ['ML0', 'ML1', 'ML7', 'MI0', 'MI1', 'MI7',
             'BL0', 'BL1', 'BL7', 'BI0', 'BI1', 'BI7']
    ahash = { 'breast milk, LAIV vaccine, 0D':0, 'breast milk, LAIV vaccine, 1D':1,
             'breast milk, LAIV vaccine, 7D':2, 'breast milk, IIV vaccine, 0D':3,
             'breast milk, IIV vaccine, 1D':4, 'breast milk, IIV vaccine, 7D':5,
             'whole blood, LAIV vaccine, 0D':6, 'whole blood, LAIV vaccine, 1D':7,
             'whole blood, LAIV vaccine, 7D':8, 'whole blood, IIV vaccine, 0D':9,
             'whole blood, IIV vaccine, 1D':10, 'whole blood, IIV vaccine, 7D':11}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPannaraj2022 = getPannaraj2022
def getOckenhouse(self, tn=1): # GSE5418
    self.prepareData("GL14")
    atype = self.h.getSurvName('c Tissue')
    atype = [bone.re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['B', 'A', 'E', 'T']
    ahash = {'experimental':2, 'treated':3, 'baseline':0, 'acute':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOckenhouse = getOckenhouse
def getMuehlenbachs(self, tn=1): # GSE7586
    self.prepareData("GL14")
    atype = self.h.getSurvName('c title')
    atype = [bone.re.sub("placenta-....-", "", str(k)) for k in atype]
    atype = [bone.re.sub("-rep.*", "", str(k)) for k in atype]
    atypes = ['UU', 'UI', 'UP', 'II', 'IP']
    ahash = {'uninflamed-infected':1, 'inflamed-infected':3, 'inflamed-past':4,
       'uninflamed-past':2, 'uninflamed-uninfected':0}
    if (tn == 2):
        atypes = ['UU', 'II']
        ahash = {'inflamed-infected':1, 'uninflamed-uninfected':0}
    if (tn == 3):
        atypes = ['U', 'I']
        ahash = {'inflamed-infected':1, 'uninflamed-uninfected':0, 'uninflamed-infected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMuehlenbachs = getMuehlenbachs
def getVentoTormo201810x(self, tn=1):
    self.prepareData("MACV361")
    atype = self.h.getSurvName('c annotation')
    atypes = sorted(bone.hu.uniq(atype[2:]))
    ahash = {}
    if (tn == 2):
        sample = self.h.getSurvName('c SampleID')
        atype = ["-".join([str(atype[i]), str(sample[i])])
                 for i in range(len(atype))]
        atype = sample
        atypes = sorted(bone.hu.uniq(atype[2:]))      
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVentoTormo201810x = getVentoTormo201810x
def getVentoTormo201810xII(self, tn=1):
    self.prepareData("MACV361.2")
    atype = self.h.getSurvName('c Fetus')
    btype = self.h.getSurvName('c location')
    atype = ["-".join([str(atype[i]), str(btype[i])])
             for i in range(len(atype))]
    atypes = sorted(bone.hu.uniq(atype[2:]))
    ahash = {}
    if (tn == 2):
        atypes = ['D6', 'D7', 'D8', 'D9', 'D10', 'D12']
        ahash = {'D6-Decidua':0, 'D7-Decidua':1, 'D8-Decidua':2,
                 'D9-Decidua':3, 'D10-Decidua':4, 'D12-Decidua':5}
    if (tn == 3):
        atypes = ['D6', 'D7', 'D8', 'D9']
        ahash = {'D6-Blood':0, 'D7-Blood':1, 'D9-Blood':3, 'D8-Blood':2}
    if (tn == 4):
        atypes = ['D8', 'D9', 'D10', 'D11', 'D12']
        ahash = {'D8-Placenta':0, 'D9-Placenta':1, 'D10-Placenta':2, 
                 'D11-Placenta':3, 'D12-Placenta':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVentoTormo201810xII = getVentoTormo201810xII

def getVentoTormo2018SS2II(self, tn=1):
    self.prepareData("MACV362.2")
    atype = self.h.getSurvName('c Fetus')
    btype = self.h.getSurvName('c location')
    atype = ["-".join([str(atype[i]), str(btype[i])])
             for i in range(len(atype))]
    atypes = sorted(bone.hu.uniq(atype[2:]))
    ahash = {}
    if (tn == 2):
        atypes = ['D1', 'D2', 'D3', 'D4', 'D5']
        ahash = {'D1-Decidua':0, 'D2-Decidua':1, 'D3-Decidua':2,
                 'D4-Decidua':3, 'D5-Decidua':4}
    if (tn == 3):
        atypes = ['D4', 'D5']
        ahash = { 'D4-Blood':0, 'D5-Blood':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVentoTormo2018SS2II = getVentoTormo2018SS2II
def getVentoTormo2018III(self, tn=1):
    self.prepareData("MACV361.3")
    atype = self.h.getSurvName('c Fetus')
    btype = self.h.getSurvName('c location')
    atype = ["-".join([str(atype[i]), str(btype[i])])
             for i in range(len(atype))]
    atypes = sorted(bone.hu.uniq(atype[2:]))
    ahash = {}
    if (tn == 2):
        atypes = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D12']
        ahash = {'D1-Decidua':0, 'D2-Decidua':1, 'D3-Decidua':2,
                 'D4-Decidua':3, 'D5-Decidua':4,
                 'D6-Decidua':5, 'D7-Decidua':6, 'D8-Decidua':7,
                 'D9-Decidua':8, 'D10-Decidua':9, 'D12-Decidua':10}
    if (tn == 3):
        atypes = ['D4', 'D5', 'D6', 'D7', 'D8', 'D9']
        ahash = { 'D4-Blood':0, 'D5-Blood':1, 'D6-Blood':2, 
                 'D7-Blood':3, 'D9-Blood':4, 'D8-Blood':5}
    if (tn == 4):
        atypes = ['D8', 'D9', 'D10', 'D11', 'D12']
        ahash = {'D8-Placenta':0, 'D9-Placenta':1, 'D10-Placenta':2, 
                 'D11-Placenta':3, 'D12-Placenta':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVentoTormo2018III = getVentoTormo2018III
def getAltman2020Blood(self, tn=1, tb=0):
    self.prepareData("MACV256")
    atype = self.h.getSurvName('c src1')
    disease = [str(k).split("-")[1] if len(str(k).split("-")) > 1
                     else None for k in atype]
    stype = [str(k).split("-")[2] if len(str(k).split("-")) > 2
                     else None for k in atype]
    atypes = ['C', 'COPD']
    ahash = {'whole blood-COPD-Control':0,
            'whole blood-COPD-COPD':1}
    if (tn == 2):
        atypes = ['C', 'Staph']
        ahash = {'whole blood-Staph-Control':0,
                'whole blood-Staph-Staph':1}
    if (tn == 3):
        atypes = ['C', 'Sepsis']
        ahash = {'whole blood-Sepsis-Control':0,
                'whole blood-Sepsis-melioidosis':1}
    if (tn == 4):
        atypes = ['C', 'TB']
        ahash = {'whole blood-TB-Control':0,
                'whole blood-TB-PTB':1}
    if (tn == 5):
        atypes = ['C', 'Melanoma']
        ahash = {'whole blood-Melanoma-Control':0,
                'whole blood-Melanoma-Melanoma':1}
    if (tn == 6):
        atypes = ['C', 'Bcell def']
        ahash = {'whole blood-B-cell deficiency-Bcell':1,
                'whole blood-B-cell deficiency-Control':0}
    if (tn == 7):
        atypes = ['C', 'Flu']
        ahash = {'whole blood-Flu-Control':0,
                'whole blood-Flu-FLU':1}
    if (tn == 8):
        atypes = ['C', 'HIV']
        ahash = {'whole blood-HIV-Control':0,
                'whole blood-HIV-HIV':1}
    if (tn == 9):
        atypes = ['C', 'JDM']
        ahash = {'whole blood-Juvenile Dermatomyositis-Control':0,
                'whole blood-Juvenile Dermatomyositis-JDM':1}
    if (tn == 10):
        atypes = ['C', 'KD']
        ahash = {'whole blood-Kawasaki-Control':0,
                'whole blood-Kawasaki-Kawasaki':1}
    if (tn == 11):
        atypes = ['C', 'Liver Transplant']
        ahash = {'whole blood-Liver Transplant-Control':0,
                'whole blood-Liver Transplant-Transplant':1}
    if (tn == 12):
        atypes = ['C', 'MS']
        ahash = {'whole blood-MS-Control':0,
                'whole blood-MS-MS Patient':1}
    if (tn == 13):
        atypes = ['C', 'Pregnancy']
        ahash = {'whole blood-Pregnancy-Control':0,
                'whole blood-Pregnancy-Pregnancy':1}
    if (tn == 14):
        atypes = ['C', 'RSV']
        ahash = {'whole blood-RSV-Control':0,
                'whole blood-RSV-RSV':1}
    if (tn == 15):
        atypes = ['C', 'SLE']
        ahash = {'whole blood-SLE-Control':0,
                'whole blood-SLE-SLE':1}
    if (tn == 16):
        atypes = ['C', 'SoJIA']
        ahash = {'whole blood-SoJIA-Control':0,
                'whole blood-SoJIA-SoJIA':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAltman2020Blood = getAltman2020Blood

def getSharma2022Mac(self, tn=1, ta=None):
    self.prepareData("MACV365")
    pdom = self.h.getSurvName('c Placental_Domain')
    infl = self.h.getSurvName('c Inflam_vasc_none')
    cell = self.h.getSurvName('c cell population')
    labor = self.h.getSurvName('c Preterm Labor')
    pterm = self.h.getSurvName('c Preterm_0no1yes')
    vpterm = self.h.getSurvName('c Very Preterm0no1yes')
    bpd = self.h.getSurvName('c BPD Classify')
    atype = cell
    atypes = ['cMNC', 'iMNC', 'ncMNC']
    ahash = {}
    if (tn == 2):
        atype = infl
        atypes = ['I', 'V']
        if ta == None:
            atype = [atype[i] if labor[i] == 'Yes'
                     else None for i in range(len(atype))]
        else:
            atype = [atype[i] if labor[i] == 'Yes' and cell[i] == ta
                     else None for i in range(len(atype))]
    if (tn == 3):
        atype = infl
        atypes = ['I', 'V']
        if ta == None:
            atype = [atype[i] if labor[i] == 'No'
                     else None for i in range(len(atype))]
        else:
            atype = [atype[i] if labor[i] == 'No' and cell[i] == ta
                     else None for i in range(len(atype))]
    if (tn == 4):
        atype = labor
        atypes = ['No', 'Yes']
    if (tn == 5):
        atype = bpd
        atypes = ['no', 'yes']
        atype = [atype[i] if pterm[i] == '1' and cell[i] == 'iMNC' and
                 labor[i] == 'No'
                 else None for i in range(len(atype))]
    if (tn == 6):
        atype = bpd
        atypes = ['no', 'yes']
        atype = [atype[i] if pterm[i] == '1' and cell[i] == 'cMNC' and
                 labor[i] == 'No'
                 else None for i in range(len(atype))]
    if (tn == 7):
        atype = bpd
        atypes = ['no', 'yes']
        atype = [atype[i] if pterm[i] == '1' and cell[i] == 'ncMNC' and
                 labor[i] == 'No'
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSharma2022Mac = getSharma2022Mac
def getJohnson2021(self, tn=1):
    self.prepareData("MACV387")
    atype = self.h.getSurvName('c listeria infection')
    atypes = ['NI', 'I']
    ahash = {'No':0, 'Yes':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJohnson2021 = getJohnson2021
def getAzari2021(self, tn=1):
    self.prepareData("MACV388")
    atype = self.h.getSurvName('c infection')
    atypes = ['NI', 'I']
    ahash = {'Non-infected (NI)':0, 'WT L. monocytogenes infected (WT)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAzari2021 = getAzari2021
def getWeisblum2017(self, tn=1):
    self.prepareData("MACV389")
    tissue = self.h.getSurvName('c tissue')
    atype = self.h.getSurvName('c infection')
    atypes = ['NI', 'I']
    ahash = {'HCMV':1, 'Mock':0, 'ZIKV':1}
    if (tn == 2):
        atype = [atype[i] if tissue[i] == 'Decidua'
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tissue[i] == 'Placenta'
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWeisblum2017 = getWeisblum2017
def getCappelletti2021(self, tn=1):
    self.prepareData("MACV390.2")
    atype = self.h.getSurvName('c treatment')
    atype = [re.sub("intra.*of ", "", str(k)) for k in atype]
    atypes = ['NI', 'I']
    ahash = {'saline':0, 'lipopolysaccharide (LPS)':1, 'live E coli':1,
             'live E coli + intra-amniotic/intra-muscular antibiotics (Abx)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCappelletti2021 = getCappelletti2021
def getYang2023Mm(self, tn=1):
    self.prepareData("MACV375.2")
    atype = self.h.getSurvName('c tissue')
    atypes = ['U', 'P']
    ahash = {'Uterus, Yolk sac and Placenta':1, 'Uterus':0,
             'Yolk sac and Placenta':1, 'Uterus (CBA/J disease model)':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2023Mm = getYang2023Mm
def getHe2021(self, tn=1, ta=0, tb=0):
    self.prepareData('MACV293')
    atype = self.h.getSurvName('c treatment')
    ahash = {'TRAM':1, 'DMSO':0, 'PANO':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c cell line')
    ahash = {'THP-1':1, '':0}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c polarization')
    atypes = ['C', 'IFNLPS', 'IL4']
    ahash = {'IL4':2, 'IFNLPS':1, 'control':0, 'Control':0}
    aval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta and sval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 3):
        atype = tval
        atypes = ['DMSO', 'TRAM', 'PANO']
        ahash = {0:0, 1:1, 2:2}
        atype = [atype[i] if aval[i] == ta and sval[i] == tb
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'IL4', 'tIL4', 'pIL4']
        ahash = {'0-0-0':0, '2-0-0':1, '2-0-1':2, '2-0-2':3}
    if (tn == 5):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'LPS', 'tLPS', 'pLPS']
        ahash = {'0-0-0':0, '1-0-0':1, '1-0-1':2, '1-0-2':3}
    if (tn == 6):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'IL4', 'tIL4', 'pIL4']
        ahash = {'0-1-0':0, '2-1-0':1, '2-1-1':2, '2-1-2':3}
    if (tn == 7):
        atype = ["-".join([str(i) for i in [aval[i], sval[i], tval[i]]])
                         for i in range(len(aval))]
        atypes = ['C', 'LPS', 'tLPS', 'pLPS']
        ahash = {'0-1-0':0, '1-1-0':1, '1-1-1':2, '1-1-2':3}
    if (tn == 8):
        atype = self.h.getSurvName('c src1')
        atypes = ['THP-1', 'PBMC']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHe2021 = getHe2021
def getJensen2011Mm(self, tn=1, ta=2):
    self.prepareData("MACV343")
    atype = self.h.getSurvName('c cell line')
    ahash = {'J774 Macrophage cell':0, 'RAW264.7 Macrophage cell':1,
             'DC2.4 Dendritic Cell':2}
    cell = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName('c infection')
    atypes = ['Un', 'Me49', 'CEP']
    ahash = {'uninfected':0, 'Me49':1, 'CEP':2}
    if (tn == 2):
        atypes = ['Un', 'I']
        ahash = {'uninfected':0, 'Me49':1, 'CEP':1}
        atype = [atype[i] if cell[i] == ta
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJensen2011Mm = getJensen2011Mm
def getMantani2023Rt(self, tn=1, ta=2):
    self.prepareData("MACV334")
    atype = self.h.getSurvName('c tissue')
    atypes = ['IL', 'PP']
    ahash = {'ileum':0, "ileal Peyer's patch":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMantani2023Rt = getMantani2023Rt
def getTorow2023Mm(self, tn=1, ta=2):
    self.prepareData("MACV335")
    atype = self.h.getSurvName('c treatment')
    atypes = ['U', 'T', 'P']
    ahash = {'untreated':0, '8h post PBS gavage':2,
             '8h post R848 gavage (0.4ug/g body weight)':1}
    if (tn == 2):
        atypes = ['U', 'T']
        ahash = {'untreated':0,
                 '8h post R848 gavage (0.4ug/g body weight)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTorow2023Mm = getTorow2023Mm
def getShulman2022Mm(self, tn=1, ta=2):
    self.prepareData("MACV337")
    atype = self.h.getSurvName('c Type')
    atypes = ['GCB', 'Mo', 'PP']
    ahash = {'GCB1_STm':0, 'GCB1_Control':0, 'GCB1_LPS':0, 'GCB2_STm':0,
             'GCB2_STm_IFNg':0, 'GCB2_Control':0, 'Monocytes_STm':1,
             'Monocytes_STm_IFNg':1, 'Monocytes_Control':1, 'PP_Control':2,
             'PP_STm':2}
    if (tn == 2):
        atypes = ['U', 'I']
        ahash = {'PP_Control':0, 'PP_STm':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getShulman2022Mm = getShulman2022Mm

def getMulder2021Tissue(self, tn=1, ta=2):
    self.prepareData("MACV399.7")
    atype = self.h.getSurvName('c Tissue')
    atypes = ['Blood', 'Spleen', 'Tonsil', 'Kidney','Breast', 'Stomach', 'Colon',
              'Liver', 'Lung', 'Pancreas', 'Skin', 'Ascites']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMulder2021Tissue = getMulder2021Tissue

def getScicluna2015(self, tn=1):
    self.prepareData("MAC104")
    atype = self.h.getSurvName('c diabetes_mellitus')
    atypes = ['No_DM', 'NA', 'DM']
    ahash = {}
    if (tn == 2):
        atypes = ['No_DM', 'DM']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getScicluna2015 = getScicluna2015

def getTakahama2024Sepsis(self, tn=1, ta=None):
    self.prepareData("MACV400")
    tissue = self.h.getSurvName('c tissue')
    treatment = self.h.getSurvName('c treatment')
    atype = tissue
    atypes = ['PBMC', 'BM', 'Spl', 'Thy', 'Ing', 'H', 'Kidney', 'SI', 'Colon',
              'Liver', 'Lung', 'Brain', 'Skin']
    ahash = {'Bone marrow':1, 'Spleen':2, 'Lung':10, 'Kidney':6, 'Heart':5,
             'Inguinal lymph node':4, 'Thymus':3, 'Liver':9, 'Colon':8,
             'Small intestine':7, 'Brain':11, 'Skin':12, 'PBMC':0}
    if tn == 2:
        atype = treatment
        atypes = ['0', '0.25', '1', '2', '3', '5']
        ahash = {'LPS_d2':3, 'LPS_d0':0, 'LPS_d3':4,
                 'LPS_d0.25':1, 'LPS_d1':2, 'LPS_d5':5}
        if ta != None:
            atype = [atype[i] if tissue[i] == ta
                    else None for i in range(len(atype))]
    if tn == 3:
        atype = treatment
        atypes = ['c0.25', 'c1', 'mi0.25', 'mi1', 'mo0.25', 'mo1', 's0.25', 's1']
        ahash = {'CLP_sham_d0.25':0, 'CLP_mild_d1':3, 'CLP_sham_d1':1,
                 'CLP_moderate_d0.25':4, 'CLP_mild_d0.25':2, 'CLP_severe_d1':7,
                 'CLP_moderate_d1':5, 'CLP_severe_d0.25':6}
        if ta != None:
            atype = [atype[i] if tissue[i] == ta
                    else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTakahama2024Sepsis = getTakahama2024Sepsis

def getGaw2023(self, tn=1):
    self.prepareData("MACV403")
    atype = self.h.getSurvName('c cell type')
    atypes = ['MIM', 'HBC']
    ahash = {'fetal Hofbauer cell':1, 'maternal intervillous monocyte':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGaw2023 = getGaw2023

def getSun2020(self, tn=1):
    self.prepareData("MACV404")
    atype = self.h.getSurvName('c tissue')
    atypes = ['D', 'P']
    ahash = {'Primary placental tissue':1, 'Primary decidua tissue':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2020 = getSun2020

def getCampbell2023scblk(self, tn=1):
    self.prepareData("MACV405.2")
    atype = self.h.getSurvName('c cell.type')
    atypes = ['fM', 'HBC', 'mM', 'm16M']
    ahash = {'Fetal CD14+ Monocytes':0, 'Fetal Hofbauer Cells':1,
             'Maternal CD14+ Monocytes':2, 'Maternal FCGR3A+ Monocytes':3}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCampbell2023scblk = getCampbell2023scblk

def getZoran2022(self, tn=1):
    self.prepareData("MACV406")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("Sp[0-9]*_", "", str(k)) for k in atype]
    atype = [re.sub("_", "", str(k)) for k in atype]
    atypes = ['Af3h', 'K6h', 'Nm3h', 'Sa3h', 'Af6h', 'K3h', 'Nm6h', 'Sa6h']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZoran2022 = getZoran2022
def getDozmorov2009(self, tn=1):
    self.prepareData("MACV408")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['mock', 'spore']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDozmorov2009 = getDozmorov2009

def getLe2019(self, tn=1, ta=None):
    self.prepareData("MACV409.3")
    cell = self.h.getSurvName('c cell_type')
    atype = self.h.getSurvName('c source_name')
    atypes = ['RPMI', 'Pa', 'Ca', 'Af', 'Mt', 'Sp', 'LPS', 'TNF', 'IL1b']
    ahash = {'Pseudomonas aeruginosa_24h':1, 'RPMI_sup':0, 'RPMI_4h':0,
             'Pseudomonas aeruginosa_4h':1, 'Candida albicans_24h':2,
             'Streptococcus pneumoniae_24h':5, 'Streptococcus pneumoniae_4h':5,
             'LPS_sup+ IL1RA+ TNF_Ab':6, 'RPMI_24h':0, 'Candida albicans_4h':2,
             'Aspergillus fumigatus_4h':3, 'RPMI':0, 'Spneu':5, 'Can_sup':2,
             'LPS_sup':6, 'Aspergillus fumigatus_24h':3, 'Can_sup+ IL1RA+ TNF_Ab':2,
             'LPS':6, 'RPMI_sup+IL1RA':0, 'Mycobacterium tuberculosis_24h':4,
             'TNF':7, 'IL-1b':8, 'Mycobacterium tuberculosis_4h':4,
             'Strep_sup':5, 'Strep_sup+ IL1RA+ TNF_Ab':5, 'Candida':2}
    if (tn == 2):
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['C', 'I']
        ahash = {'RPMI_24h':0, 'Candida albicans_24h':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLe2019 = getLe2019

def getNorvershtern2011(self, tn=1):
    self.prepareData("G19")
    atype = self.h.getSurvName('c cell type')
    atypes = ['Other', 'Neu', 'NeuM', 'Eos', 'Bas', 'NK']
    ahash = {'Granulocyte (Neutrophil)':1,
             'Granulocyte (Neutrophilic Metamyelocyte)':2,
             'Eosinophill':3, 'Basophils':4,
            'Mature NK cell_CD56+ CD16+ CD3-':5}
    for title in atype[2:]:
        if title not in ahash:
            ahash[title] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNorvershtern2011 = getNorvershtern2011

def getAllantaz2011(self, tn=1):
    self.prepareData("NEU6", cfile = "/Users/rohan/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c cell type (ch1)")
    atypes = ['B', 'T', 'Neu', "Mo", 'Eos', 'NK']
    ahash = {'CD19+ B cells':0, 'CD14+ monocytes':3, 'CD4+ T cells':1,
             'CD8+ T cells':1, 'Eosinophils':4, 'NK cells':5, 'Neutrophils':2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAllantaz2011 = getAllantaz2011

def getMonaco2017(self, tn=1):
    self.prepareData("NEU12", cfile = "/Users/rohan/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c cell type (ch1)")
    atypes = ['Other', 'Neu', 'Mo', 'Bas', 'NK']
    ahash = {'Classical monocytes':2, 'Intermediate monocytes':2,
             'Non classical monocytes':2, 'Natural killer cells':4,
             'Low-density neutrophils':1, 'Low-density basophils':3}
    for title in atype[2:]:
        if title not in ahash:
            ahash[title] = 0
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMonaco2017 = getMonaco2017

def getTharp2023(self, tn=1):
    self.prepareData("MACV325")
    atype = self.h.getSurvName("c cell type")
    atypes = ['WL', 'Neu', 'Mac']
    ahash = {'Whole lung':0, 'Interstitial Macrophage':2,
             'Neutrophil':1, 'Alveolar Macrophage':2}
    if (tn == 2):
        atypes = ['iM', 'AM']
        ahash = {'Interstitial Macrophage':0,
                 'Alveolar Macrophage':1}
    if (tn == 3):
        atypes = ['Neu', 'AM']
        ahash = {'Neutrophil':0,
                 'Alveolar Macrophage':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTharp2023 = getTharp2023
def getBacac2023(self, tn=1):
    self.prepareData("MACV330")
    atype = self.h.getSurvName("c timepoint")
    atypes = ['4h', '20h']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBacac2023 = getBacac2023

def getSKY(self, tn=1):
    self.prepareData("SKY1")
    atype = self.h.getSurvName("c Type")
    atypes = ['Pre', 'Post']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSKY = getSKY

def getChandran2021(self, tn=1):
    self.prepareData("MACV466")
    atype = self.h.getSurvName("c treatment")
    atypes = ['T1', 'T2', 'T3', 'T4']
    ahash = {'T1 (baseline)':0, 'T4 (three months after meditation)':3,
            'T3 (after meditation)':2, 'T2 (before meditation)':1}
    if (tn == 2):
        atypes = ['T2', 'T3']
        ahash = {'T3 (after meditation)':1, 'T2 (before meditation)':0}
    if (tn == 3):
        atypes = ['T2', 'T4']
        ahash = {'T4 (three months after meditation)':1,
                'T2 (before meditation)':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChandran2021 = getChandran2021

def getAggarwal2022(self, tn=1):
    self.prepareData("MACV467")
    atype = self.h.getSurvName("c prakriti")
    atypes = ['V', 'K', 'PK', 'VK', 'VP', 'P']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAggarwal2022 = getAggarwal2022

def getQu2013(self, tn=1):
    self.prepareData("MACV468")
    atype = self.h.getSurvName("c treatment")
    atypes = ['CB', 'CA', 'YB', 'YA']
    ahash = {'before Yoga intervention':2,
            'after Control intervention':1,
            'after Yoga intervention':3,
            'before Control intervention':0}
    if (tn == 2):
        atypes = ['YB', 'YA']
        ahash = {'before Yoga intervention':0,
                'after Yoga intervention':1}
    if (tn == 3):
        atypes = ['CB', 'CA']
        ahash = {'after Control intervention':1,
                'before Control intervention':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQu2013 = getQu2013

def getHitreck2023(self, tn=1):
    self.prepareData("MUSCLE19")
    atype = self.h.getSurvName("c Type")
    atypes = ['R', 'nR', 'KR', 'KnR']
    ahash = {'MUT SVE':2, 'MUT non-SVE':3, 'WT SVE':0, 'WT non-SVE':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHitreck2023 = getHitreck2023

def getArjaans2020(self, tn=1):
    self.prepareData("MACV469")
    atype = self.h.getSurvName("c bpd")
    atypes = ['N', 'Mild', 'Mod', 'S']
    ahash = {'Severe':3, 'Mild':1, 'Moderate':2, 'None':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getArjaans2020 = getArjaans2020

def getMestan2024(self, tn=1):
    self.prepareData("MACV470")
    atype = self.h.getSurvName("c SampleType")
    atypes = ['S', 'QC', 'Cal', 'B']
    ahash = {'Sample':0, 'QC':1, 'Calibrator':2, 'Buffer':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    if (tn == 2):
        atype = self.h.getSurvName("c BPD_NIH _DICHOT")
        atypes = ['N', 'BPD']
        ahash = {'BPD':1, 'No BPD':0}
        atype = [atype[i] if tval[i] == 0
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMestan2024 = getMestan2024

def getMestan2024II(self, tn=1):
    self.prepareData("MACV470.2")
    bpd = self.h.getSurvName("c BPD_NIH _DICHOT")
    atype = bpd
    atypes = ['N', 'BPD']
    ahash = {'BPD':1, 'No BPD':0}
    if (tn == 2):
        atype = self.h.getSurvName("c BPD_NIH_5CAT")
        atypes = ['M', 'S']
        ahash = {'sev':1, 'mild':0}
    if (tn == 3):
        atype = self.h.getSurvName("c km_ai")
        atypes = ['BPDnoAI', 'BPDai']
        ahash = {'0.0':1, '1.0':0}
        atype = [atype[i] if bpd[i] == 'BPD'
                else None for i in range(len(atype))]
    if (tn == 4):
        atype = self.h.getSurvName("c km_mvm")
        atypes = ['BPDnomvm', 'BPDmvm']
        ahash = {'0.0':1, '1.0':0}
        atype = [atype[i] if bpd[i] == 'BPD'
                else None for i in range(len(atype))]
    if (tn == 5):
        atype = self.h.getSurvName("c primarypathdomain")
        atypes = ['MVM', 'NONE']
        ahash = {}
    if (tn == 6):
        atype = self.h.getSurvName("c km_ai")
        atypes = ['noAI', 'ai']
        ahash = {'0.0':1, '1.0':0}
    if (tn == 7):
        atype = self.h.getSurvName("c primarypathdomain")
        atypes = ['NONE', 'CI', 'AI', 'MVM', 'MIXED']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMestan2024II = getMestan2024II

def getJackson2023kd(self, tn=1):
    self.prepareData("COV432")
    atype = self.h.getSurvName('c Characteristics[disease]')
    atypes = ['V', 'B', 'KD', 'U']
    ahash = {'adenovirus infection':0,
            'escherichia coli infection':1,
            'group A streptococcal infection':1,
            'influenza':0,
            'juvenile idiopathic arthritis':3,
            'Kawasaki disease':2,
            'malaria':1,
            'meningococcal infection':3,
            'pneumococcal infection':3,
            'rhinovirus infection':0,
            'Respiratory Syncytial Virus Infection':0,
            'Staphylococcus aureus infection':1,
            'Tuberculosis':1}
    if (tn == 2):
        atypes = ['B', 'KD']
        ahash = {'Kawasaki disease':1, 'escherichia coli infection':0,
                'group A streptococcal infection':0,
                'Staphylococcus aureus infection':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJackson2023kd = getJackson2023kd

def getJackson2023kdII(self, tn=1):
    self.prepareData("COV432.2")
    atype = self.h.getSurvName('c Characteristics[disease]')
    atypes = ['N', 'V', 'B', 'KD', 'MISC']
    ahash = {'Kawasaki disease':3,
            'viral disease':1,
             'normal':0,
             'bacterial diseasel':2,
             'multisystem inflammatory syndrome in children':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJackson2023kdII = getJackson2023kdII

def getLance2007(self, tn=1, tb=0):
    self.prepareData("LP3")
    treatment = self.h.getSurvName('c Treatment')
    place = self.h.getSurvName('c Collection            Place ')
    bpd = self.h.getSurvName('c BPD                  Yes/No')
    time = self.h.getSurvName('c Time           (Letter)')
    sid = self.h.getSurvName('c Original               Sample Name')
    sid = [re.sub("[A-Z]\s*$", "", str(k)) for k in sid]
    sid = [re.sub("\s", "", str(k)) for k in sid]
    shash = {'M01':'I', 'M07':'I', 'M10': 'V', 'M11':'I',
             'M13': 'I', 'M14': 'I', 'M17': 'I', 'M18': 'I',
             'GM02': 'I', 'GM03': 'I', 'GM08': 'V', 'GM14': 'V', 'GM16': 'V',
             'GM17': 'N', 'GM19': 'I', 'GM21': 'I', 'GM22': 'I', 'GM26': 'I'}
    atype = treatment
    atypes = ['CTRL', 'LPS']
    ahash = {}
    if (tn == 2):
        atype = place
        atypes = ['VAND', 'GM', 'GMR', 'UCSD']
    if (tn == 3):
        atype = [atype[i] if place[i] == 'GMR'
                 else None for i in range(len(atype))]
    if (tn == 4):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if time[i] == 'A'
                 else None for i in range(len(atype))]
        #atype = [atype[i] if time[i] != 'A'
        #         else None for i in range(len(atype))]
    if (tn == 5):
        atype = [treatment[i] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['CTRL', 'LPS']
    if (tn == 6):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V']
    if (tn == 7):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V']
        atype = [atype[i] if time[i] != 'A' and time[i] != 'B'
                 else None for i in range(len(atype))]
    if (tn == 8):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V']
        atype = [atype[i] if time[i] == 'A' or time[i] == 'B'
                 else None for i in range(len(atype))]
    if (tn == 9):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V', 'N']
        atype = [atype[i] if time[i] != 'A'
                 else None for i in range(len(atype))]
    if (tn == 10):
        atype = ['AB' if time[i] == 'A' or time[i] == 'B'
                 else 'C' for i in range(len(sid))]
        atypes = ['AB', 'C']
    if (tn == 11):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if time[i] == 'A' or time[i] == 'B'
                 else None for i in range(len(atype))]
    if (tn == 12):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if time[i] != 'A' and time[i] != 'B'
                 else None for i in range(len(atype))]
    if (tn == 13):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if (time[i] == 'A' or time[i] == 'B') and
                treatment[i] == tb
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLance2007 = getLance2007

def getBell2023(self, tn=1, ta = 'CordBlood'):
    self.prepareData("MACV367")
    atype = self.h.getSurvName("c variable-disease")
    atypes = ['nonBPD', 'BPD']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c variable-time")
        atypes = ['CB', 'D14', 'D28']
        ahash = {'Day14':1, 'CordBlood':0, 'Day28':2}
    if (tn == 3):
        timepoint = self.h.getSurvName("c variable-time")
        atype = [atype[i] if timepoint[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 4):
        o2 = self.h.getSurvName("c variable-o2 treatment")
        atype = [str(atype[i])+'-'+str(o2[i])for i in range(len(atype))]
        timepoint = self.h.getSurvName("c variable-time")
        atype = [atype[i] if timepoint[i] == ta
                 else None for i in range(len(atype))]
        atypes = ['nonBPD-no', 'BPD-yes', 'nonBPD-yes']
        ahash = {}
    if (tn == 5):
        o2 = self.h.getSurvName("c variable-o2 treatment")
        atype = [str(atype[i])+'-'+str(o2[i])for i in range(len(atype))]
        timepoint = self.h.getSurvName("c variable-time")
        atype = [atype[i] if timepoint[i] == ta
                 else None for i in range(len(atype))]
        atypes = ['nonBPD-no', 'BPD-yes']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBell2023 = getBell2023
def getWindhorst2023(self, tn=1):
    self.prepareData("MACV368")
    atype = self.h.getSurvName("c disease state")
    atypes = ['No', 'Mild', 'M/S']
    ahash = {'No BPD':0, 'Mild BPD':1, 'Moderate/severe BPD':2}
    if (tn == 2):
        weight = [0,0] + ana.h.getSurvName("c birth weight")[2:]
        atypes = ['No', 'M/S']
        ahash = {'No BPD':0, 'Moderate/severe BPD':1}
        atype = [atype[i] if ((atype[i] == 'No BPD') and (float(weight[i]) < 1300)) or
                 ((atype[i] == 'Moderate/severe BPD') and (float(weight[i]) < 1000))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWindhorst2023 = getWindhorst2023
def getWang2022(self, tn=1):
    self.prepareData("MACV369")
    atype = self.h.getSurvName("c disease status")
    atypes = ['nonBPD', 'BPD']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2022 = getWang2022
def getRyan2019(self, tn=1):
    self.prepareData("MACV370")
    atype = self.h.getSurvName("c src1")
    atypes = ['NC', 'BC', 'ND', 'BD']
    ahash = {'No_BPD_DHA_peripheral blood':2,
             'BPD_DHA_peripheral blood':3,
             'No_BPD_Control_peripheral blood':0,
             'BPD_Control_peripheral blood':1}
    if (tn == 2):
        atypes = ['noBPD', 'BPD']
        ahash = {'No_BPD_Control_peripheral blood':0,
                 'BPD_Control_peripheral blood':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRyan2019 = getRyan2019
def getCheng2018(self, tn=1):
    self.prepareData("MACV372")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['C', 'BPD']
    ahash = {'normal':0, 'Bronchopulmonary dysplasia (BPD)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCheng2018 = getCheng2018

def getEldredge2024(self, tn=1, ta=0):
    self.prepareData("MACV506")
    atype = self.h.getSurvName("c cell type")
    atypes = ['E', 'M']
    ahash = {'cd14+cd16- monocytes':1, 'cd14+cd16+ monocytes':1,
              'epithelial':0}
    if tn == 2:
        cell = self.h.getSurvName("c cell type")
        ahash = {'cd14+cd16- monocytes':1, 'cd14+cd16+ monocytes':2,
                  'epithelial':0}
        ctype = [ahash[i] if i in ahash else None for i in cell]
        atype = self.h.getSurvName("c severity")
        atypes = ['M', 'S']
        ahash = {'severe I':1, 'severe II':1, 'moderate':0, 'mild':0}
        atype = [atype[i] if ctype[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atypes = ['E', 'M']
        ahash = {'cd14+cd16- monocytes':1,
                  'epithelial':0}
    if tn == 4:
        atypes = ['E', 'M']
        ahash = {'cd14+cd16+ monocytes':1,
                  'epithelial':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEldredge2024 = getEldredge2024

def getZec2023mm(self, tn=1, ta=0):
    self.prepareData("MACV515")
    atype = self.h.getSurvName("c src1")
    atypes = ['M0', 'M1']
    ahash = {'wildtype':1, 'wildtype_mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZec2023mm = getZec2023mm

def getLi2024mm(self, tn=1, ta=0):
    self.prepareData("MACV516")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1']
    ahash = {'no':0, 'IL-4':1, 'IL-4/ARV825':1,
            'IL-4/DMSO':1, 'IL-4/JQ1':1, 'IL-4/ZL0420':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mm = getLi2024mm

def getLi2024mmII(self, tn=1, ta=0):
    self.prepareData("MACV516.2")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M2']
    ahash = {'none':0, 'IL-4/ARV825':1, 'IL-4/ZL0420':1,
            'IL-4':1, 'IL-4/DMSO':1, 'IL-4/JQ1':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mmII = getLi2024mmII

def getLi2024mmIII(self, tn=1, ta=0):
    self.prepareData("MACV516.3")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub("[-_].*", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {}
    if tn == 2:
        atype = self.h.getSurvName("c Title")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1-1':1, 'M1-2':1, 'M1-3':1,
                'M2-1':2, 'M2-2':2, 'M2-3':2,
                'M0-a':0, 'M0-b':0, 'M0-c':0}
    if tn == 3:
        atype = self.h.getSurvName("c Title")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1-1':1, 'M1-2':1, 'M1-3':1,
                'M2-1':2, 'M2-2':2, 'M2-3':2}
    if tn == 4:
        atype = self.h.getSurvName("c Title")
        atypes = ['M1', 'M2']
        ahash = {'M1-1':0, 'M1-2':0, 'M1-3':0,
                'M2-1':1, 'M2-2':1, 'M2-3':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mmIII = getLi2024mmIII

def getDusek2008(self, tn=1, ta=0):
    self.prepareData("BL11")
    atype = self.h.getSurvName("c src1")
    atypes = ['N1', 'N2', 'M']
    ahash = {'Healthy Subject with No RR Practice':0,
            'Healthy Subject with 8 weeks of  RR Practice':1,
            'Long-term practitioner of daily RR practice':2}
    if tn == 2:
        atype = self.h.getSurvName("c Title")
        n1 = [34, 35, 32, 46, 39, 11, 12, 7, 47, 30, 33, 23, 3, 2, 1,
                41, 40, 6, 14]
        n2 = [6, 14, 3, 7, 35, 40, 41, 46, 47, 30, 34, 11, 12, 39, 23, 1,
                2, 32, 36]
        m = [24, 10, 25, 27, 13, 38, 29, 16, 18, 19, 15, 20, 21, 5, 22, 28,
                42, 43, 44]
        ahash = {}
        for k in n1:
            ahash[f"N1-{k:02d}"] = 0
        for k in n2:
            ahash[f"N2-{k:02d}"] = 1
        for k in m:
            ahash[f"M-{k:02d}"] = 2
    if tn == 3:
        atypes = ['N1', 'N2']
        ahash = {'Healthy Subject with No RR Practice':0,
                'Healthy Subject with 8 weeks of  RR Practice':1}
    if tn == 4:
        atypes = ['N1', 'M']
        ahash = {'Healthy Subject with No RR Practice':0,
                'Long-term practitioner of daily RR practice':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDusek2008 = getDusek2008

def getBuric2020(self, tn=1, ta=0):
    self.prepareData("BL12")
    atype = self.h.getSurvName("c sample group")
    atypes = ['C-Pre', 'C-Post', 'M-Pre', 'M-Post', 'Y-Pre', 'Y-Post']
    ahash = {'Control post-intervention':1,
            'Control pre-intervention':0,
            'mindfulness post-intervention':3,
            'mindfulness pre-intervention':2,
            'yoga post-intervention':5,
            'yoga pre-intervention':4}
    if tn == 2:
        atypes = ['M-Pre', 'M-Post']
        ahash = {'mindfulness post-intervention':1,
                'mindfulness pre-intervention':0}
    if tn == 3:
        atypes = ['Y-Pre', 'Y-Post']
        ahash = {'yoga post-intervention':1,
                'yoga pre-intervention':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBuric2020 = getBuric2020

def getKuo2015(self, tn=1, ta=0):
    self.prepareData("BL14")
    atype = self.h.getSurvName("c disease state")
    ahash = {'IBS':0, 'IBD':1}
    disease = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName("c time point")
    atypes = ['Pre', 'Post']
    ahash = {'after 9 weeks of intervention':1,
            'before intervention':0,
            'baseline before intervention':0}
    if tn == 2:
        atype = [atype[i] if disease[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKuo2015 = getKuo2015

def getSun2025BPD(self, tn=1, ta=0):
    self.prepareData("LP9")
    disease = self.h.getSurvName("c disease")
    atype = disease
    atypes = ['control', 'hBPD', 'aeBPD', 'eBPD']
    ahash = {}
    if tn == 2:
        atypes = ['C', 'BPD']
        ahash = {'control':0, 'hBPD':1, 'aeBPD':1, 'eBPD':1}
    if tn == 3:
        atype = self.h.getSurvName('c age')
        btype = self.h.getSurvName('c age-unit')
        atype = [str(atype[i])+"-"+str(btype[i]) for i in range(len(atype))]
        atypes = ['I', 'T', 'C']
        ahash = {'4-month':0, '5-month':0, '6-month':0, '7-month':0,
                '8-month':0, '316-day':0, '11-month':0, '12-month':1,
                '13-month':1, '14-month':1, '15-month':1, '18-month':1,
                '19-month':1, '20-month':1, '21-month':1, '3-year':2}
        atype = [atype[i] if disease[i] == 'control'
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = self.h.getSurvName('c sex')
        atypes = ['M', 'F']
        ahash = {'female':1, 'male':0}
        atype = [atype[i] if disease[i] == 'control'
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2025BPD = getSun2025BPD

def getSun2025BPDII(self, tn=1, ta=0):
    self.prepareData("LP9.2")
    cell = self.h.getSurvName("c predicted.cellref_celltype")
    disease = self.h.getSurvName("c disease")
    atype = disease
    atypes = ['control', 'hBPD', 'aeBPD', 'eBPD']
    ahash = {}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
        atype = [atype[i] if disease[i] == 'control'
                 else None for i in range(len(atype))]
    if tn == 6:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC2'))
                 else None for i in range(len(atype))]
    if tn == 7:
        atype = [atype[i] if ((cell[i] == 'cDC1') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2025BPDII = getSun2025BPDII

def getLiu2024bpd(self, tn=1, ta=0):
    self.prepareData("LP10")
    atype = self.h.getSurvName('c disease state')
    atype = [hu.re.sub(",.*", "", str(k)) for k in atype]
    atypes = ['H', 'BPD']
    ahash = {'healthy':0, 'Bronchopulmonary dysplasia (BPD)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiu2024bpd = getLiu2024bpd


def getMaag2017(self, tn=1):
    self.prepareData("BE35", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c Characteristics[disease]')
    atypes = ['N', 'BE', 'EAC']
    ahash = {'normal':0, "Barrett's esophagus":1,
             'esophageal adenocarcinoma ':2}
    if (tn == 2):
        atypes = ['BE', 'EAC']
        ahash = {"Barrett's esophagus":0,
                 'esophageal adenocarcinoma ':1}
    if (tn == 3):
        atypes = ['N', 'BE']
        ahash = {'normal':0, "Barrett's esophagus":1}
    if (tn == 4):
        atypes = ['N', 'EAC']
        ahash = {'normal':0,
                 'esophageal adenocarcinoma ':1}
    if tn == 5:
        atype = self.h.getSurvName('c Disease')
        ahash = {'normal':0, "Barrett's esophagus low-grade dysplasia":2,
                 "Barrett's esophagus non-dysplastic":1,
                 'esophageal adenocarcinoma':3}
        atypes = ['N', 'NDBE', 'BE-D', 'EAC']
    if tn == 6:
        atype = self.h.getSurvName('c Disease')
        ahash = {'normal':0, "Barrett's esophagus non-dysplastic":1}
        atypes = ['N', 'NDBE']
    if tn == 7:
        atype = self.h.getSurvName('c Disease')
        ahash = {"Barrett's esophagus low-grade dysplasia":1,
                 "Barrett's esophagus non-dysplastic":0}
        atypes = ['NDBE', 'BE-D']
    if tn == 8:
        atype = self.h.getSurvName('c Disease')
        ahash = {"Barrett's esophagus low-grade dysplasia":0,
                 'esophageal adenocarcinoma':1}
        atypes = ['BE-D', 'EAC']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMaag2017 = getMaag2017

def getGuo2022cellrefHs(self, tn=1, ta=0):
    self.prepareData("LP17")
    cell = self.h.getSurvName("c celltype_level3")
    age = self.h.getSurvName("c Age")
    atype = age
    atypes = ['Y', 'M', 'O']
    ahash = {}
    for i in self.h.aRange():
        if atype[i] == '':
            continue
        if float(atype[i]) < 45:
            ahash[atype[i]] = 0
        elif float(atype[i]) < 60:
            ahash[atype[i]] = 1
        else:
            ahash[atype[i]] = 2
    aval = [ahash[k] if k in ahash else None for k in atype]
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
        atype = [atype[i] if aval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuo2022cellrefHs = getGuo2022cellrefHs

def getGuo2022cellrefMm(self, tn=1, ta=0):
    self.prepareData("LP18")
    cell = self.h.getSurvName("c celltype_level3")
    age = self.h.getSurvName("c Time")
    atype = age
    atypes = ['E16.5', 'E18.5', 'PND1', 'PND3', 'PND7', 'PND10',
            'PND14', 'PND28']
    ahash = {}
    aval = [ahash[k] if k in ahash else None for k in atype]
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuo2022cellrefMm = getGuo2022cellrefMm

def getVanneste2023lungAMmm(self, tn=1, ta=0):
    self.prepareData("BL34")
    atype = self.h.getSurvName("c src1")
    atypes = ['Mo', 'AM', 'IM']
    ahash = {'Ly6C+ Mo':0, 'AM':1,
            'CD206+ control IM':2, 'CD206+ refilled IM':2,
            'CD206- refilled IM':2, 'CD206- control IM':2}
    if tn == 2:
        atypes = ['IM', 'AM']
        ahash = {'AM':1,
                'CD206+ control IM':0, 'CD206+ refilled IM':0,
                'CD206- refilled IM':0, 'CD206- control IM':0}
    if tn == 3:
        atypes = ['AM', 'CD206-IM', 'CD206+IM']
        ahash = {'AM':0,
                'CD206+ control IM':2, 'CD206+ refilled IM':2,
                'CD206- refilled IM':1, 'CD206- control IM':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getVanneste2023lungAMmm = getVanneste2023lungAMmm

def getLi2024mBALMm(self, tn=1, ta=0):
    self.prepareData("MACV523.2")
    atype = self.h.getSurvName("c Cell.Types")
    atypes = ['AM', 'IM', 'Lym', 'DC', 'Epi']
    ahash = {'AMs':0, 'IMs':1, 'Cyc.AMs':0, 'Hb.AMs':0, 'Lym':2,
            'Res.DCs':3, 'Mig.DCs':3, 'Epi':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mBALMm = getLi2024mBALMm

def getZahalka2022Mm(self, tn=1, ta=0):
    self.prepareData("BL63")
    atype = self.h.getSurvName("c Type")
    atypes = ['C', 'LPS', 'HISP', 'HISP-LPS']
    ahash = {'lps_hisp':3, 'lps_medium':1, 'ctrl_hisp':2, 'ctrl_medium':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZahalka2022Mm = getZahalka2022Mm

def getRothchild2024Mm(self, tn=1, ta=0):
    self.prepareData("BL64")
    atype = self.h.getSurvName("c Type")
    atypes = ['AM-P-N', 'AM-P-P', 'AM-L-N', 'AM-L-P',
            'BM', 'BM-L', 'BM-P-P', 'BM-L-N', 'BM-L-P']
    ahash = {
            'AM, PBS-bead bystander':0,
            'AM, PBS-bead positive':1,
            'AM, LPS-bead bystander':2,
            'AM, LPS-bead positive':3,
            'BMDM, PBS-bead positive':6,
            'BMDM, LPS-bead bystander':7,
            'BMDM, LPS-bead positive':8,
            'BMDM, untreated':4,
            'BMDM, soluble LPS':5}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRothchild2024Mm = getRothchild2024Mm

def getSanku2024MoDC(self, tn=1, ta=0):
    self.prepareData("BL65")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'Mf', 'EV', 'EV-', 'E/S']
    ahash = {'unconditioned':0, 'Live Mf':1,
            'Mf-derived extracellular vesicles (EVs)':2,
            'EV depleted E/S':3,
            'Mf-derived excretory/secretory (E/S) products':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSanku2024MoDC = getSanku2024MoDC

def getMorrow2019copd(self, tn=1, ta=0):
    self.prepareData("MACV305")
    atype = self.h.getSurvName('c cell type')
    ahash = {'alveolar macrophage':0, 'bronchial epithelium':1,
            'whole blood':2}
    cell = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName("c copd")
    atypes = ['cont', 'case']
    ahash = {}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = cell
        atypes = ['AM', 'BE', 'WB']
        ahash = {0:0, 1:1, 2:2}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMorrow2019copd = getMorrow2019copd

def getBlackburn2022copd(self, tn=1, ta=0):
    self.prepareData("LU22")
    atype = self.h.getSurvName("c Sex")
    atypes = ['Male', 'Female']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlackburn2022copd = getBlackburn2022copd

def getBlackburn2022copdII(self, tn=1, ta=0):
    self.prepareData("LU22.2")
    cell = self.h.getSurvName("c predicted.cellref_celltype.pruned")
    atype = self.h.getSurvName("c Sex")
    atypes = ['Male', 'Female']
    ahash = {}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
    if tn == 6:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC2'))
                 else None for i in range(len(atype))]
    if tn == 7:
        atype = [atype[i] if ((cell[i] == 'cDC1') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlackburn2022copdII = getBlackburn2022copdII

def getLv2024copd(self, tn=1, ta=0):
    self.prepareData("LU23")
    atype = self.h.getSurvName("c Title")
    disease = [hu.re.sub("[0-9].*", "", str(k)) for k in atype]
    atype = disease
    atypes = ['CTL', 'COPD']
    ahash = {}
    if tn == 2:
        atype = self.h.headers
        ahash = {'Male':0, 'GSM8576150':1, 'GSM8576152':1, 'GSM8576149':1,
                'GSM8576143':1,'GSM8576148':1,'GSM8576155':1,'GSM8576145':1}
        atype = [atype[i] if atype[i] in ahash
                 else 'Male' for i in range(len(atype))]
        atypes = ['M', 'F']
        #atype = [atype[i] if disease[i] == 'COPD'
        #         else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLv2024copd = getLv2024copd

def getLv2024copdII(self, tn=1, ta=0):
    self.prepareData("LU23.2")
    cell = self.h.getSurvName("c predicted.cellref_celltype.pruned")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub("[0-9].*", "", str(k)) for k in atype]
    atypes = ['CTL', 'COPD']
    ahash = {}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
    if tn == 6:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC2'))
                 else None for i in range(len(atype))]
    if tn == 7:
        atype = [atype[i] if ((cell[i] == 'cDC1') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLv2024copdII = getLv2024copdII

def getHeo2025copdBlood(self, tn=1, ta=0):
    self.prepareData("LU24")
    atype = self.h.getSurvName("c group")
    atypes = ['Control', 'COPD']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHeo2025copdBlood = getHeo2025copdBlood

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

def getMoon2024copdII(self, tn=1, ta=0):
    self.prepareData("LU25.2")
    cell = self.h.getSurvName("c predicted.cellref_celltype.pruned")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub(".* C", "C", str(k)) for k in atype]
    atype = [hu.re.sub(", scRNA", "", str(k)) for k in atype]
    atypes = ['Ct', 'COPDt', 'Co', 'COPDo']
    ahash = {'Control, Tissue':0, 'COPD, Tissue':1, 'Control, Organoid':2,
            'COPD, Organoid':3, 'COPD, Organid':3}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
    if tn == 6:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC2'))
                 else None for i in range(len(atype))]
    if tn == 7:
        atype = [atype[i] if ((cell[i] == 'cDC1') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMoon2024copdII = getMoon2024copdII

def getWang2024copd(self, tn=1, ta=0):
    self.prepareData("LU26")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub(".$", "", str(k)) for k in atype]
    atypes = ['C', 'COPD', 'Ct', 'COPDt']
    ahash = {'CTumor':3, 'Control':0, 'COPD':1, 'Tumor':2}
    if tn == 2:
        atypes = ['C', 'COPD']
        ahash = {'Control':0, 'COPD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2024copd = getWang2024copd

def getWang2024copdII(self, tn=1, ta=0):
    self.prepareData("LU26.2")
    cell = self.h.getSurvName("c predicted.cellref_celltype.pruned")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub(".$", "", str(k)) for k in atype]
    atypes = ['C', 'COPD', 'Ct', 'COPDt']
    ahash = {'CTumor':3, 'Control':0, 'COPD':1, 'Tumor':2}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if ((cell[i] == 'AM') or (cell[i] == 'IM'))
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC1') or
            (cell[i] == 'cDC2') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    if tn == 5:
        atype = cell
        atypes = ['pMON', 'iMON', 'AM', 'IM', 'pDC', 'cDC1', 'cDC2', 'maDC']
    if tn == 6:
        atype = [atype[i] if ((cell[i] == 'pDC') or (cell[i] == 'cDC2'))
                 else None for i in range(len(atype))]
    if tn == 7:
        atype = [atype[i] if ((cell[i] == 'cDC1') or (cell[i] == 'maDC'))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2024copdII = getWang2024copdII

def getAdams2020copd(self, tn=1, ta=0):
    self.prepareData("LU27")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'COPD', 'IPF']
    ahash = {}
    if tn == 2:
        atypes = ['Control', 'COPD']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAdams2020copd = getAdams2020copd

def getAdams2020copdII(self, tn=1, ta=0):
    self.prepareData("LU27.2")
    cell = self.h.getSurvName("c Subclass_Cell_Identity")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'COPD', 'IPF']
    ahash = {}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAdams2020copdII = getAdams2020copdII

def getAdams2020copdIII(self, tn=1, ta=0):
    self.prepareData("LU27.3")
    cell = self.h.getSurvName("c Manuscript_Identity")
    atype = self.h.getSurvName("c disease")
    atypes = ['Control', 'COPD', 'IPF']
    ahash = {}
    if tn == 2:
        atype = [atype[i] if cell[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAdams2020copdIII = getAdams2020copdIII

def getCao2023lungFetal(self, tn=1, ta=0):
    self.prepareData("LU28")
    atype = self.h.getSurvName("c time")
    atypes = ['4w', '5w', '6w', '7w', '7w+', '8w']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCao2023lungFetal = getCao2023lungFetal

def getHe2022lungFetal(self, tn=1, ta=0):
    self.prepareData("LU29")
    atype = self.h.getSurvName("c Characteristics[age]")
    atypes = ['T1', 'T2']
    ahash = {'17':1, '13':1, '15':1, '14':1, '6':0, '21':1, '11':0, '20':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHe2022lungFetal = getHe2022lungFetal

def getBarnes2023lungFetal(self, tn=1, ta=0):
    self.prepareData("LU30")
    atype = self.h.getSurvName("c Characteristics[age]")
    atypes = ['T1', 'T2', 'NA']
    ahash = {'12 pcw':0, 'CS23+8pcw':0, '20 pcw':1, '9pcw':0, '':2, '9 pcw':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBarnes2023lungFetal = getBarnes2023lungFetal

def getGuo2019lungMm(self, tn=1, ta=0):
    self.prepareData("LU31")
    atype = self.h.getSurvName("c age")
    atypes = ['Pre', 'Post']
    ahash = {'E16.5':0, 'E18.5':0, 'PND1':1, 'PND3':1,
            'PND7':1, 'PND14':1, 'PND28':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGuo2019lungMm = getGuo2019lungMm

def getTCGAPANCAN(self, tn=1):
    self.prepareData("GL44.4")
    atype = self.h.getSurvName("c sample_type")
    atypes = ['N', 'T']
    ahash = {'Primary Tumor':1, 'Solid Tissue Normal':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getTCGAPANCAN = getTCGAPANCAN

def getGTEX(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("GTEx1", "/Users/aglina/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c SEX');
    ahash = {'1':1, '2':0}
    atypes = ['F','M']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGTEX = getGTEX
 
def getSohail2024(ana, tn=1):
    ana.prepareData("LU41")
    atype = ana.h.getSurvName("c infection group")
    ahash = {'Uninfected':0,
            'Influenza A virus (IAV)':2,
            'Mycobacterium bovis (BCG)':3,
            'Pseudomonas aeruginosa (PA)':1}
    atypes = ['U', 'PA', 'IAV', 'BCG']
    if tn == 2:
        ahash = {'Uninfected':0,
                'Pseudomonas aeruginosa (PA)':1}
        atypes = ['U', 'PA']
    ana.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSohail2024 = getSohail2024

def getMingTing2024(ana, tn=1):
    ana.prepareData("PLP249")
    atype = ana.h.getSurvName("c treatment")
    ahash = {'untreated':0, 'non-AIEC T75':1, 'AIEC 541-15':2}
    atypes = ['U', 'nAIEC', 'AIEC']
    if tn == 2:
        ahash = {'untreated':0, 'AIEC 541-15':1}
        atypes = ['U', 'AIEC']
    ana.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMingTing2024 = getMingTing2024

def getKabiljo2024(ana, tn=1):
    ana.prepareData("BL37")
    atype = ana.h.getSurvName("c src1")
    ahash = {'Macrophages, patient-derived colorectal cancer-associated fibroblasts':1,
            'Macrophages, patient-derived colorectal cancer-associated fibroblasts, patient-derived colorectal cancer organoids':2,
            'Monocytes, patient-derived colorectal cancer-associated fibroblasts':0}
    atypes = ['Mo', 'MoC', 'MoCF']
    ana.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKabiljo2024 = getKabiljo2024

def getWilhelm2025paAM(ana, tn=1):
    ana.prepareData("LU42")
    atype = ana.h.getSurvName("c time")
    ahash = {}
    atypes = ['0h', '8d', '24h', '96h']
    atypes = ['0h', '8d']
    ana.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWilhelm2025paAM = getWilhelm2025paAM

def getGEO2025THP1(ana, tn=1):
    ana.prepareData("BL121")
    atype = ana.h.getSurvName("c PMA")
    cell_line = ana.h.getSurvName("c cell line")
    ahash = {'Y':1, 'N':0}
    atypes = ['Mo', 'Mac']
    atype = [atype[i] if cell_line[i] == 'THP1'
            else None for i in range(len(atype))] 
    ana.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGEO2025THP1 = getGEO2025THP1

