#! /Users/croth/anaconda/bin/python
## Import the first 2 needed mods
import argparse
import datetime
today = str(datetime.date.today())
## Pares the arguments.   
parser = argparse.ArgumentParser()
parser.add_argument("inputf",
                    type=str,
                    help="Name of (or path) to VCF file.");
parser.add_argument("-outputf",
                    type=str,
                    help="The output file name that will contain data from this analysis.",
                    default = 'FromVCFtoGeno'+today+'.csv');
parser.add_argument("chl",
                    type=str,
                    help="Name of file (or path) containing the lengths of the chromosomes.");
parser.add_argument("-badc",
                    type=str,
                    help="Name of file (or path) containing the list of segregatns with anaumolus chromosomes. Default: '../Crypto-Bad-Chrom-Apr2016.csv'");
parser.add_argument("-QD",
                    type=int,
                    help="The normilized qulaity threshold value for each site. Default: 20",
                    default = 20);
parser.add_argument("-X",
                    type=int,
                    help="Read coverage threshold per variant site per segregant. Default: 15X",
                    default = 15);
parser.add_argument("-arf",
                    type=float,
                    help="The minumum threhold for ratio between reads suggesting a variant and total read depth. Default: 0.75",
                    default = 0.75);
parser.add_argument("-af",
                    type=float,
                    help="The minimum allelic frequency threshold. Used to filter out rare variants. Default: 0.1",
                    default = .1);
parser.add_argument("-howm",
                    type=float,
                    help="The needed prcentage of mapping populaiton with a given variant. Default: 0.9",
                    default =.9);
args = parser.parse_args();
## Import other needed mods
import os.path
import math
import numpy as np
import pandas as pd
import cryptofxns
import vcf
#import scipy.stats
#from matplotlib import pyplot as plt
#from itertools import chain
## Patch values. 
## These are needed files and values.
## In the future their nessecity may change.
inputfile = args.inputf #'../QTL_MERGE_AUGUST15.vcf'
outputfile = args.outputf #'../Crypto-Genotypes-Apr2016.csv'
chromlenpath = args.chl #'../XL280-chromlen.txt'
QDthres = args.QD
covX = args.X
threshold = args.arf
p = args.af
howm = args.howm
badchrpath = args.badc # "../Crypto-Bad-Chrom-Apr2016.csv"
## List the path to the file name. 
## VCF_FILE = '../QTL_MERGE_AUGUST15.vcf'
## Parse the VCF File
rdr = vcf.Reader(open(os.path.expanduser(inputfile), 'r'))
## For the records in the VCF file
recs = [rec for rec in rdr]
## How many records are in the vcf file?
print "Total Number of variant calls:", len(recs)
## What are the chromosome names and lengths? 
## Load chrom lengths
## chromlenpath = '../XL280-chromlen.txt'
clens = np.loadtxt(chromlenpath, dtype=np.int).tolist()
print "Lengths of Chromosomes: ", clens
## Get unique chromosome names
chrlist = list(np.unique([str(rec.CHROM) for rec in recs]));
print("Chromosomes: %s " % chrlist)
## What are the sample names?
## Get sample names
sample_names = [seg.sample for seg in recs[0]]
## For each record add the strain name as an info field. 
## In addition to some of the site paramaters, we will edit this field instead of the records in the VCF file. 
for rec in recs:
    rec.add_info('segregants',[seg.sample for seg in rec][:]) # <- this makes a "slice" that we can edit for each rec.
## Basic Filtering: Take only sites that are SNPs 
## Take only SNPs from our data (removing indels)
## A "global" filter
snps = [rec for rec in recs if rec.is_snp]
print "Sites that are SNPs:",len(snps)
## there shouldn't be any sites where all samples are fixed for the ref allele
snps = [rec for rec in snps if not rec.is_monomorphic]
print "Sites taht are non-monomorphic SNPs: ", len(snps)
## remove any non-heterozygotic sites, since this is a (primarily) haploid organism
snps = [rec for rec in snps if not (rec.num_het > 0)]
print "Sites that are non-heterozygotic SNPs:", len(snps)
## We may want to return to these indels for analysis. 
## Calculate the number of indels
#indels = [indel for indel in recs if indel.is_indel]
#print "Number of indels:",len(indels)
## Check to see that number of SNPs and Indels add to original size
#print len(recs) - (len(indels) + len(snps))
## Consider only bialleiac SNPs: only one alternate allele for segreagant. 
## Here we want to test only a site who has one altenate allele. 
## Another global filter
## Primarily b/c that is what our QTL model is designed to test for. 
#bisites = [rec for rec in snps if (rec.call_rate == 1.0) and (len(rec.ALT) == 1)]
## I've changed this code to allow for sites that are not called for everyone. 
## Instead we should set a minimum number of individuals that we need to be called later on.
bisites = [rec for rec in snps if (len(rec.ALT) == 1)]
print "Sites that are biallelic SNPs:",len(bisites)
## Delete some memory hogs!
## get rid of memory hogs we dont need anymore
del recs
del snps
## Aneuploidy and Disomy.
## From our draw SNP program we saw several characters with either disomic or possessing aneuploidy. For now, we have taken note of these segregants and listed them in a file. 
## badchrpath = "../Crypto-Bad-Chrom-Apr2016.csv"
BADchr = pd.read_csv(badchrpath)
## On the level of the segregant, filter out anomalous sites at bad chromosomes. 
def bad_chrom_out(Bisites,badchr): 
    FiltSNPs = [];
    for rec in Bisites:
    ## If rec in bad chrom, take only the individuals who are "normal" 
        if rec.CHROM in badchr.columns: 
            for seg in rec: 
                if badchr[rec.CHROM].str.contains(seg.sample).any() == True and seg.sample in rec.INFO['segregants']:
                    rec.INFO['segregants'].remove(seg.sample)
            FiltSNPs.append(rec)
        else:
            FiltSNPs.append(rec)
    return FiltSNPs
## call the ftn above
filtSNPs = bad_chrom_out(bisites,BADchr);
## Check the number of records in the filtered data
#len(filtSNPs)
## Filter SNP sites based of normilized quality scores
#QDthres = 20;
FiltSNP = [rec for rec in filtSNPs if rec.INFO['QD'] > QDthres];
print "Sites with quality score > ", QDthres, ": ", len(FiltSNP)
## Extract Biallelic SNPs with high coverage
## Sites with <= 15X coverage for all individuals. This guarantees probability of ~80% correct calling. Don't know if that is true but sounds good. Will double check. 
## This is a filter on the individual domain. 
#covX = 15.0;
hcbiallelic = [];
for rec in FiltSNP:
    for seg in rec: 
        if seg.sample in rec.INFO['segregants'] and seg.data.DP < covX:
            rec.INFO['segregants'].remove(seg.sample)
        elif seg.sample in rec.INFO['segregants'] and seg.data.DP == None:
            rec.INFO['segregants'].remove(seg.sample)
    hcbiallelic.append(rec)
print "Sites with high coverage and biallelic: ", len(hcbiallelic), " with X >= ", covX
## Take sites with 75% of the reads are mapping a SNP
## The allelic read depth must be 80% of the reads required to map a site. 
#threshold = 0.75
realsnps = [];
for rec in hcbiallelic:
    for seg in rec:
        if seg.sample in rec.INFO['segregants'] and float(seg.data.AD[int(seg.data.GT)]) / int(seg.data.DP) < threshold:
            rec.INFO['segregants'].remove(seg.sample)
    realsnps.append(rec)
print "Sites with reliable allelic read depth SNPs: ",len(realsnps)
## Take SNPs where ~90% of the mapping population are called. 
#howm = .9
numpop = round(len(sample_names) * howm)
print numpop, "segregants per site to ~ 90% of the mapping population. "
goodsnps = [];
for rec in realsnps:
    #print len(rec.INFO['segregants'])
    if len(rec.INFO['segregants']) >= numpop:
        goodsnps.append(rec)
print "Sites with ~90% of the population called: ", len(goodsnps)
## Take only sites that vary. 
## 10% of segregants in the population should have an allele.
## This filter removes rare alleles
## Filter out any fixed or private alleles
## This is another example of a global filter scheme
#M = .1
#thres = round(len(straininfo) * M) # ~ 11 individuals. 
varsites = [];
for rec in goodsnps:
    indgenotype = [int(seg.data.GT) for seg in rec if seg.sample in rec.INFO['segregants']];
    #print indgenotype
    if float(sum(indgenotype))/len(rec.INFO['segregants']) >= p:
        varsites.append(rec)
print "Sites with allelic variation: ", len(varsites)
## Create and DataFrame with Genotype Info
## Here we will want to reset an entry to "NA" if they are disomic or anueploidic. 
varsnps = varsites
chroms = [i.CHROM for i in varsnps]
pos = [i.POS for i in varsnps]
refs = [str(i.REF) for i in varsnps]
alts = [str(i.ALT[0]) for i in varsnps]
## Make data frame
genoDF = pd.DataFrame()
genoDF['Chrom'] = chroms
genoDF['Pos'] = pos
genoDF['Ref'] = refs
genoDF['Alt'] = alts
## list samplenames
samplenames = [i.sample for i in varsnps[0]] ## DO NOT CHANGE THIS. ORDER IS IMPORTANT!
## Code genomes
allgenos = [cryptofxns.code_genotypes(rec, coding=(1,0)) for rec in varsnps]
allgenoDF = pd.DataFrame(allgenos)
allgenoDF.columns = samplenames
## Combine dataframes 
combinedDF = pd.concat([genoDF,allgenoDF], axis=1)
#combinedDF.shape
## Save Filtered genotype data
#outputfile = '../Crypto-Genotypes-Apr2016.csv';
## Write meta data from this analysis (includes thresholds used and dates ran).
# Write the function to do so here.
def appendsave(DF,file,args):
    f = open(file,'w');
    #s = '#test\n'
    s = '# '+ today + ', ' + 'VCF: ' + str(args.inputf) + ', ' + 'Chromfile: ' + str(args.chl) + ', ' +  'Norm QD: ' + str(args.QD) + ', ' +'RD: ' + str(args.X) + ', ' + 'Allelic RD Freq: ' + str(args.arf) + ', ' + 'Min Allelic Freq: ' + str(args.af) + ', ' + '% Map Pop: ' + str(args.howm) + ', ' + 'File of Anomalies: ' + str(args.badc) + '\n'; 
    f.write(s)
    f.close()
    with open(file,'a') as f:
        DF.to_csv(f,index=False)
    f.close()
    pass 
#combinedDF.to_csv(outputfile, index=False)
appendsave(combinedDF,outputfile,args)
## Write a new, filtered VCF file!
## The Writer class provides a way of writing a VCF file. 
## Currently, you must specify a template Reader which provides the metadata:
#vcf_reader = vcf.Reader(filename=VCF_FILE)
#vcf_writer = vcf.Writer(open('../QTL_FILTERED_JAN16.vcf', 'w'), rdr)
#for rec in varsnps:
#    vcf_writer.write_record(rec)