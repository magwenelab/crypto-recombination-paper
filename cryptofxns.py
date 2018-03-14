import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from itertools import cycle

## File paths
FILEPATHgf = '../FILES/GenotoFiltGeno2016-04-07.csv'
FILEPATHfgf = '../FILES/FiltGenotoHapFiltGeno2017-02-22.csv'
FILEPATHsi = '../FILES/Strain_list_Cryptococcus.xlsx'

## Colors used for ploting
colors = ["brown", "forestgreen", "olive", "orange", "maroon", "lightcoral", "navy" , "red"]

## Needed variables. 
thecol = 4; # the column in the genotype df where seg genotypes start. 
myinfo = ["Strain name","Cross"]

## The current (Jan 22, 2016) genotype file. 
geno_df = pd.read_csv(FILEPATHgf,skiprows=0,header=1);

## remove following samples that aren't in the genotype data 
## or are near identical segregants
geno_remove = [#"SS-C045",     # not yet genotyped
               #"SS-C040",     #  " "  
               #"SS-C030_cor", #  " "
                "SS-B377",     # nearly identical to SS-B312
                "SS-B385",     # nearly identical to SS-B317
                "SS-B410",     # nearly identical to SS-B314
                "SS-B411",     # nearly identical to SS-B395 and diplody Chr01
                "SS-B560",     # nearly identical to SS-B308
                "SS-B395",     # Chr01    
                "SS-B594",     # Chr01
                "SS-C273",     # Chr07
                "SS-B879",     # Chr10 
                "SS-B885",     # Chr10
                #"XL280",       # We dont need the folks for this analysis Chr10
                #"XL280a",
                #"SS-A837",
                "SS-B593"] # Non-recombinant All chrom from XL280 except Chr07"

axalp = [['SS-B869_Correction', 'SS-B872_Correction', 'SS-B873_Correction', 'SS-B874_Correction'], # 1
         ['SS-B876', 'SS-B880','SS-B879', 'SS-B885'], #2## SS-B885 and B879 are identical, except chrom 4 and disomic on chrom 10
         ['SS-B886_Correction', 'SS-B887_Correction'],#3
         ['SS-B890_Correction', 'SS-B896_Correction', 'SS-B898_Correction', 'SS-B892_Correction',#4 ## SS-B892_cor 58.5% like SS-B930
         'SS-B901_Correction', 'SS-B905', 'SS-B906_Correction', 'SS-B908_Correction'],#4
        ['SS-B926', 'SS-B927', 'SS-B929', 'SS-B930'],#5
         ['SS-B952', 'SS-B956'],#6
         ['SS-B960_Correction', 'SS-B961_Correction'],#7 ## SS-B960_cor
        ['SS-B996', 'SS-B997', 'SS-B998'],#8
        ['SS-C001'],#9
        ['SS-C026_Correction', 'SS-C029_Correction'],#10
        ['SS-C031_Correction'],#11
        ['SS-C033', 'SS-C036', 'SS-C039'],#12
        ['SS-C242'],#13
        ['SS-C271', 'SS-C272'],#14
        ['SS-C273'],#15
        ['SS-C290', 'SS-C291']]#16

## Length of each chromosome. 
chromosome_lengths = [2300525, 1632286, 2105722, 1783051, 1507536, 1438940, 1347778, 1194303, 1178642, 1085711, 1019830, 906713, 787979, 762695];

## Location of ssk1 on chromosome 2
ssk1_start = 930980 ## (-)
ssk1_end = 935461

## Location of erg6 on chromosome 2
erg6_start = 935614 ## (+)
erg6_end = 937336

## Location of the MAT locus on chromosome 4
mats = 1529240 ## 5' 
matst = 1661861 ## 3'

## Location of ura5 on chromosome 7
ura5_start = 1047561 ## (-)
ura5_end = 1048380

## Our scaling factors 
scales = [2,1,2.0/3,.5]

## Set Directories and Paths as variables 
folders = ['IMI/','BPS/','HAP/','FGT/','FIG/'];

CHRPATHS = ['../FILES/CHR' + str(i) +'/' for i in range(1,len(chromosome_lengths)+1)];

FILECHRPATHS = [ [chp + f for f in folders] for chp in CHRPATHS];

## Functions used in analysis of recombination frequencies. 
## =================================================================== ##


def load_centrom_locs(filepath='../FILES/JEC21_cens.tsv',
    cent_col=['Name','Start','Stop']):
    """
    Loads dataframe with locations of centromiers for each chromosome
    """
    temp_cent = pd.read_csv(filepath,sep='\t',header=None)
    temp_cent.columns = cent_col
    centros = temp_cent.ix[:,1:]
    centros = centros.astype(int)
    centros.columns = [0,1]
    centros.index = temp_cent.Name.tolist()
    return centros

def IMI(w1,v2):
   """Intermarker: For consecutive markers [v1,w1] and [v2,w2],
   the intermarker interval (IMI) is defined as: 
   [x,y] = [w1 + 1, v2 - 1 ] in Z. The nominal size is then defined as
   y - x + 1. """
   return ( int(v2) - 1) - ( int(w1) + 1 )  + 1 

def haplotypes(imi,clen,method='min'):
    """
    Calculates the size of the Crossover or non-crosover gene conversion track. 
    Aka the size of the haplotype. 
    Inputs: IMI, a DataFrame containing intermarkers itervals with 4 columns:
            v: starting marker position of IMI
            w: ending marker position of IMI
            X: Number of detected genotype changes
            IMI: The inter marker interval size. 

            CLEN, the length of the chromosome.

            METHOD, the method used for calculating the 5' and 3' positions of each haplotype.
            default is 'min'

    Output: HAPDF, a DataFrame containing haplotype information with 3 columns:
            v: 5' starting position of the haplotype. 
            w: 3' ending position of the haplotype. 
            nb: The size of the haplotype.  
    """
    ## Find breakpoints
    bp = imi[imi.X != 0].dropna(axis=0,how='any');
    ## 1st and last marker positions
    pos1 = int(imi.v[imi.index.tolist()[0]]);
    pose = int(imi.w[imi.index.tolist()[-1]]);
    ## Starting and ending markers in haplotypes
    start = [pos1] + bp.w.tolist()
    stops = bp.v.tolist() + [pose]
    ## List all markers involved adding 0 and chromosome length. 
    markers = sorted([0] + start + stops + [clen])
    ## Find intervals and size via minimal method:
    hapdf = pd.DataFrame()
    if method is 'min':
        hapdf['v'] = [markers[i] for i in range(1,len(markers)-1,2)];
        hapdf['w'] = [markers[i+1] for i in range(1,len(markers)-1,2)];
    elif method is 'max':
        hapdf['v'] = [markers[i-1]+1 for i in range(1,len(markers)-1,2)];
        hapdf['w'] = [markers[i+2]-1 for i in range(1,len(markers)-1,2)];
    else: ## Method is midpoint
        hapdf['v'] = [(markers[i-1]+markers[i])*.5 for i in range(1,len(markers)-1,2)];
        hapdf['w'] = [(markers[i+1]+markers[i+2])*.5 for i in range(1,len(markers)-1,2)];
    hapdf['nb'] = [int(hapdf.w[k] - hapdf.v[k] + 1) for k in hapdf.index.tolist()];
    return hapdf

def haplotype_filter(hapdf,k1):
    """
    Filters haplotype dataframe by size.
    
    Inputs: HAPDF, a dataframe with a chromosomal haplotype information. 
            This includes columns named,
            "v" beginning of haplotypes (left to right) chromosomal position in bases. 
            "w" ending of haplotypes in chromosomal positons. 
            "nb" The size of the haplotype.
            "genotype" the genotype that acures most often for the markers within a haplotype.

            K1, the cutoff size for the haplotypes. Haplotypes < K1 are filtered out. 
            
    Output: HAP, a dataframe with filtered chromosomal haplotypes and the same columns, 
            start, stop, size, and genotype. 
    """
    hap = hapdf.copy() 
    while min(hap.nb.tolist()) < k1: ## If the smallest haplotype is less than our filter
        ri = hap[hap.nb == min(hap.nb.tolist())].index[0]
        if ri == 0: ## First Haplotype
            temp = hap.ix[ri:ri+1,:];
            assert temp.genotype[ri] != temp.genotype[ri+1]
            nv = temp.v[ri];
            nw = temp.w[ri+1];
            nsize = nw - nv + 1;
            hap.ix[ri] = [nv,nw,nsize,abs(hap.genotype[ri] - 1)]
            hap.drop([ri+1],inplace=True)
        elif ri == len(hap) - 1: ## Last Haplotype (index zero based)
            temp = hap.ix[ri-1:ri,:];
            assert temp.genotype[ri] != temp.genotype[ri-1]
            nv = temp.v[ri-1];
            nw = temp.w[ri];
            nsize = nw - nv + 1;
            hap.ix[ri] = [nv,nw,nsize,abs(hap.genotype[ri] - 1)]
            hap.drop([ri-1],inplace=True)
        else: ## Any Haplotype between 1st and last haplotype
            temp = hap.ix[ri-1:ri+1,:];
            nv = temp.v[ri-1];
            nw = temp.w[ri+1];
            nsize = nw - nv + 1;
            hap.ix[ri] = [nv,nw,nsize,abs(hap.genotype[ri] - 1)]
            hap.drop([ri-1,ri+1],inplace=True)
        ## Reset index Calculate the smallest haplotype and run again. 
        hap.reset_index(drop = True,inplace=True)
    return hap

def strain_info(sample_names,info,path = None):
    """Pairs strain names with information from 
    the excel file"""
    ## Define variables
    if path is None:
        path = FILEPATHsi
    ## Load in the strain list excel file
    xl = pd.ExcelFile(path);
    ## Parse the file
    straindf = xl.parse("Sheet1");
    ## put the sample names into a dataframe
    samplenames = pd.DataFrame(sample_names,
        columns = ['Strain_name']);
    ## take the common indicies. 
    strainsdf = straindf[straindf["Strain name"].isin(samplenames.Strain_name)].reset_index(drop = True);
    ## take only the needed info from the dataframe
    #mystrainsdf[info];
    ## return 
    return strainsdf

def progeny_list(geno_df,tc=thecol,mi=myinfo):
    """
    Lists and seperate the progeny by the cross they came from.

    Input: Genotype dataframe with segregatns names.

    Output: Two lists with smaple names. One for the progeny 
    from the uni-sexual cross (alpha X alpha) and bi-sexual 
    cross (a X alpha). 
    """
    straininfo = strain_info(geno_df.columns[tc:].tolist(),mi) 
    ## Alpha X Alpha Progeny
    alpxalp = straininfo['Strain name'][
        straininfo.Cross == str(
        straininfo['Cross'].unique()[0:2][1])].tolist()
    ## A X Alpha Progeny 
    axalp = straininfo['Strain name'][
    straininfo.Cross == str(
        straininfo['Cross'].unique()[0:2][0])].tolist()
    return(alpxalp,axalp)

def allelic_ANOVA(site, pheno):
    """This regression is equivalent to one-way ANOVA with 2 groups. Return F-statistic.

    Assumes site is coded as 0, 1
    """
    coding = np.array(site, np.float)
    pheno = np.array(pheno, np.float)
    
    meany = np.mean(pheno)
    meandummy = np.mean(coding)
    ctry = pheno - meany
    ctrdummy = coding - meandummy
    
    # regression coefficient and intercept
    b = np.dot(ctry, ctrdummy)/np.dot(ctrdummy, ctrdummy)
    intercept = meany - b * meandummy
    yhat = b * ctrdummy
    len_yhat = np.sqrt(np.dot(yhat,yhat))
    len_y = np.sqrt(np.dot(ctry,ctry))
    df_yhat = 1
    
    error = ctry  - yhat
    len_error = np.sqrt(np.dot(error,error))
    if abs(len_error**2) < 1e-5:
        #print(Exception("Zero length error in ANOVA"))#raise Exception("Zero length error in ANOVA")
        F = np.nan
    else:
        df_error = len(pheno) - 2
        # F-statistic
        F = (len_yhat**2/df_yhat) / (len_error**2/df_error)
    return F

def association_logPval(site, pheno):
    #pheno = pheno.ix[site.dropna().index];assert len(pheno) > 0;
    #site = site.dropna();assert len(site) > 0;
    pheno = pheno.dropna();assert len(pheno) > 0;
    site = site[pheno.dropna().index.tolist()];assert len(site) == len(pheno)
    F = allelic_ANOVA(site, pheno)
    if np.isnan(F):
        return F
    else:
        return -np.log10(ss.f.sf(F, 1, len(pheno)-2))

## FTNs for drawing SNPs
def chrom_offset(chromlens):
    return [sum(chromlens[:i]) for i in range(len(chromlens))]

def jitter(x,mult=0.025):
    "Adds noise"
    sdx = np.std(x);
    return np.array(x) + mult * np.random.uniform(-sdx, sdx, len(x))

def draw_QTLs(df, phenolabel, chrlist, clens, ylabel = "$-\log_{10}\, p$",
              colors=colors, alpha=0.3,mk = '.', ls = ' ',
              mylines=None,fprefix=None,lw =.3,yscale = None, Close = True, DPI = 600):

    clrcycle = cycle(colors)
    offset_dict = dict(zip(chrlist, chrom_offset(clens)))
    combcoord = df.Pos + np.array([offset_dict[i] for i in df.Chrom.tolist()])
        
    fig, ax = plt.subplots(4,4,figsize=(16,16), dpi=150)
    ax1 = plt.subplot2grid((4,4), (0,0), colspan=2)
    
    for i, chrom in enumerate(chrlist):
        sub = df.query("Chrom == '%s'" % chrom)
        clr = clrcycle.next()
        ax1.plot(sub.Pos + offset_dict[chrom], sub[phenolabel], 
            marker = mk, linestyle=ls, color=clr, alpha=alpha)
        plt.sca(ax.ravel()[i+2])
        plt.plot(sub.Pos, sub[phenolabel], 
            marker=mk, linestyle=ls, color=clr, alpha=alpha)
        if mylines is not None:
            for they in mylines:
                plt.hlines(they,min(sub.Pos),max(sub.Pos),linestyle='--',linewidth=lw)
        plt.ylim(0,1.05*np.max(df[phenolabel]))
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.text(0.9, 0.95,chrom, horizontalalignment='center',
                 verticalalignment='center',
                 transform = plt.gca().transAxes)
        plt.ylabel(ylabel, fontsize=14)
        
    chrom_midpts = 0.5 *  np.sum(np.array(zip([0]+list(np.cumsum(clens)),
                                      np.cumsum(clens))),axis=1)
    ax1.set_xticks(chrom_midpts)
    ax1.set_xticklabels([str(i) for i in range(1, len(clens)+1)])
    ax1.set_xlabel("Chromosomes")
    ax1.set_ylabel(ylabel)
    if yscale is not None:
        ax1.set_ylim(0,yscale)
    plt.suptitle(phenolabel, fontsize=18);
    plt.tight_layout();
    plt.subplots_adjust(top=0.95)
    if fprefix is not None:
        fig.savefig(fprefix + '-{}.pdf'.format(phenolabel), dpi=DPI)
    if Close is True:
        plt.close()

## FTNs for collapsing genotype information by blocks of recombination
def marker_blocks(geno,chrlist,tc=4,ref=0):
    chromblocks = []
    ## Sub sample geno by chromosome
    for chrom in chrlist:
        test = geno[geno.Chrom==chrom].reset_index(drop=True);
        endin = test.index.tolist()[-1];
        subtractdf = abs(test.ix[:endin-1,tc:].reset_index(drop=True) - test.ix[1:,tc:].reset_index(drop=True));
        subtract = subtractdf.sum(axis=1);
        ## Set blocks
        blocks = []
        block = []
        blockin = []
        for i,a in enumerate(subtract.tolist()):
            block.append(a)
            blockin.append(subtract.index.tolist()[i])
            if sum(block) > 0:
                blocks.append(blockin)
                block = []
                blockin = []
        ## Check work 
        for block in blocks:
            for a in test.ix[block,tc:].sum(axis=0).tolist():
                assert a == ref or a == len(block)
        chromblocks.append(blocks)
    return chromblocks

def collapse_geno(geno,chromblocks,chrlist,tc=4,ref=0):
    segs = geno.columns.tolist()[tc:]
    temp= [];
    for ch,chrom in enumerate(chrlist):
        sub = geno[geno.Chrom == chrom].reset_index(drop=True)
        blocks = chromblocks[ch]
        for block in blocks:
            v = int(sub.ix[block[0],'Pos'])
            w = int(sub.ix[block[-1],'Pos'])
            gts = sub.ix[block[0],segs].tolist()
            assert len(gts) == len(segs)
            temp.append([chrom,v,w]+gts)
    return pd.DataFrame(temp,columns=['Chrom','v','w']+segs);

def hotspots(x,mu):
    return ss.poisson.pmf(x,mu)#-np.log10(ss.poisson.pmf(x,mu))
## Approach for IMI analysis as seen in Stienmetz et al. (This approach is some what biased towards the bicross)
def IMIandCON(df,chrlist=None,chrlab='Chrom',poslab='Pos',tc=4,unit_scale = 1e-6):
    dfs = []
    if chrlist is None:
    	chrlist = sorted(df.Chrom.unique().tolist())
    for chrom in chrlist:
        sub = df[df[chrlab] == chrom].reset_index(drop=True)
        N = float(len(sub.columns.tolist()[tc:]))
        sub_in = sub.index.tolist()
        IMIs = [IMI(int(sub[poslab][i]),int(sub[poslab][i+1])) for i in range(len(sub_in[:-1]))]
        temp = sub.ix[:sub_in[-1],tc:].reset_index(drop=True) - sub.ix[1:,tc:].reset_index(drop=True)
        CONs = temp.abs().sum(axis=1);
        cdf = pd.DataFrame(columns=['Chrom','W','V','IMI','XO'],index=sub_in[:-1])
        cdf['Chrom'] = chrom
        cdf['W'] = sub[poslab].tolist()[:-1]
        cdf['V'] = sub[poslab].tolist()[1:]
        cdf['IMI'] = IMIs
        cdf['XO'] = CONs
        dfs.append(cdf)
    return pd.concat(dfs).reset_index(drop=True)

def collapse_IMI(df,Threshold = 500):
    chrlist = sorted(df.Chrom.unique().tolist())
    assert len(chrlist) == 14
    newdfs = []
    for chrom in chrlist:
        newsub = []
        sub = df[df.Chrom == chrom].reset_index(drop=True).copy()
        subin = sub.index.tolist()
        combinein = []
        for i in range(len(subin[:-1])):
            if (int(sub.W.tolist()[i+1]) - int(sub.V.tolist()[i])) <= Threshold:
                combinein.append(i+1)
        for i in combinein[::-1]:
            sub.ix[i-1,'V'] = sub.ix[i,'V']
            sub.ix[i-1,'IMI'] = IMI(sub.ix[i-1,'W'],sub.ix[i,'V'])
            sub.ix[i-1,'XO'] = sub.ix[i-1,'XO'] + sub.ix[i,'XO']
            sub.ix[i-1,'bXO'] = sub.ix[i-1,'bXO'] + sub.ix[i,'bXO']
            sub.ix[i-1,'uXO'] = sub.ix[i-1,'uXO'] + sub.ix[i,'uXO']
            sub.ix[i-1,'Pval'] = max(sub.ix[i-1,'Pval'],sub.ix[i,'Pval'])
            sub.ix[i-1,'bPval'] = max(sub.ix[i-1,'bPval'],sub.ix[i,'bPval'])
            sub.ix[i-1,'uPval'] = max(sub.ix[i-1,'uPval'],sub.ix[i,'uPval'])
            sub.drop(i,inplace=True)
        newdfs.append(sub)
    return pd.concat(newdfs).reset_index(drop=True)