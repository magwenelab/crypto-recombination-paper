#! /Users/croth/anaconda/bin/python
## Import the first 2 needed mods
import argparse
import datetime
today = str(datetime.date.today())
## Pares the arguments.   
parser = argparse.ArgumentParser()
parser.add_argument("inputf",
                    type=str,
                    help="Name of (or path) to CSV file.");
parser.add_argument("-outputf",
                    type=str,
                    help="The output file name that will contain data from this analysis.",
                    default = 'GenotoFiltGeno'+today+'.csv');
parser.add_argument("strainfile",
                    type=str,
                    help="Name of file (or path) containing the list of segregatns with mating info. Default: '../Strain_list_Cryptococcus.xlsx'",
                    default = '../Strain_list_Cryptococcus.xlsx');
parser.add_argument("parents",
                    nargs='+',
                    help="The parents used in the cross. If there are 3 this implies both a unisexual and bisexual cross.");
parser.add_argument("-donor",
                    type=str,
                    help="Segregant to use as template genotype if parent is missing.",
                    default = 'SS-B960_cor');
parser.add_argument("-MatPos",
                    type=float,
                    help="The starting position for analysis and copying from donor of the Matting locus. Default: 1.4 * 10**6",
                    default = 1.4 * 10**6);
parser.add_argument("-straininfo",
                   nargs='+',
                    help="The information to be taken out of the strain file. Default: '['Strain name','Cross']'",
                    default = ["Strain name","Cross"]);
args = parser.parse_args();
import math
import os.path
import numpy as np
import pandas as pd
import cryptofxns
## Load Genotype Data
## The genotype data was created, filtered, and saved using the first file in the pipeline, "Crypto-VCFtoGeno-Mar2016.ipynb".
## Load genotype data
## args.CSVpath = "../Crypto-Genotypes-filter-Feb2016.csv"
geno = pd.read_csv(args.inputf,skiprows=0,header=1)
print "Genotype data dimensions:", geno.shape
#geno.head()
## Creation of synthetic parent: XL280a
## Add to the data a clone of XL280 only for chromosome 4 it will need to have synthetic "a" mating locus. This locus will be copied from SS-B960_cor. My analysis of the "wood floor" plots has shown a unique "a" mating locus signiture. I have concluded that the donor SS-B960_cor is not recombinant from chromosomal position 1.4* 10**6 to the end of chromosome IV. This interval includes the "a" mating locus. 
## List who the parents are.
#parentA = 'XL280'
#parenta = 'XL280a'
#parentB = 'SS-A837'
print args.parents
print type(args.parents)
if len(np.unique(args.parents)) == 2:
    parentA = args.parents[0];
    parentB = args.parents[1];
elif len(np.unique(args.parents)) == 3:
    parentA = args.parents[0];
    parenta = args.parents[1];    
    parentB = args.parents[2];
print parentA
print type(parentA)
## if we have a doner ... 
## Slice and copy the data fro XL280
list1 = ['Chrom','Pos','Ref','Alt'] + [parentA];
list1a = ['Chrom','Pos','Ref','Alt'] + [parenta];
XL280a = geno[['Chrom','Pos','Ref','Alt',parentA]][:]
## Change the name of the column
XL280a.columns = [list1a]
## check our slice
#XL280a.head()
## Take chr04 for the segregant named SS-B960_cor.
## This segregant was chosen from examining the recombination spots and a locus signiture
donor = args.donor#'SS-B960_cor'
## from diagnostic plots 
chr04donor = geno[['Chrom','Pos','Ref','Alt',donor]][geno.Chrom == 'Chr04'][:]
## check this slice
#chr04donor.head()
## This is the position in SS-B960_cor's chr04 that has the a locus 
## and no recombination from this point to the end of the chromosome. 
matlocus_position = args.MatPos#1.4 * 10**6
## Take the index of those SNP's
chr04alocus_index = chr04donor[chr04donor.Pos > matlocus_position].index;
## print the number of SNP's
#print len(chr04alocus_index)
## Updatea the values in the XL280a df with the new "a" synthetic locus SNPs
XL280a.ix[chr04alocus_index,'XL280a']  = chr04donor.ix[chr04alocus_index,donor];
## Add it to the larger combinedDF
genoDF = pd.concat([geno, XL280a.XL280a], axis=1, join_axes=[geno.index]);
## Check if the above worked
print genoDF.shape
## More Data filtering
## After examining the SNPs between our two parents, we need to drop any site that is invariant between our three parental lines SS-A837 and XL280, sepereatly according to mating type.
## Progeny According to Mating. 
## The progeny that are from XL280alpha x SS-A4837alpha (unisexual reproduction) cross are those between B307 and B621, and the progeny from XL280a x SS-A837alpha (bisexual reproduction) cross are those between B869 and C291.
myinfo = args.straininfo#["Strain name","Cross"]
straininfo = cryptofxns.strain_info(genoDF.columns[4:].tolist(),myinfo,args.strainfile)
alpxalp = straininfo['Strain name'][straininfo.Cross == str(straininfo['Cross'].unique()[0:2][1])].tolist()
axalp = straininfo['Strain name'][straininfo.Cross == str(straininfo['Cross'].unique()[0:2][0])].tolist()
print len(alpxalp)," progeny from alpha x alpha cross"
print len(axalp)," progeny from a x alpha cross"
## Seperate the progeny
info1DF = genoDF[['Chrom','Pos','Ref','Alt',parentA,parentB]][:];
info2DF = genoDF[['Chrom','Pos','Ref','Alt',parenta,parentB]][:];
## Make new data frames. 
unis = pd.concat([info1DF,genoDF[alpxalp]], axis=1, join_axes=[info1DF.index])
bis = pd.concat([info2DF,genoDF[axalp]], axis=1, join_axes=[info2DF.index])
## Remove any sites where two parents are not different
unisDF = unis[unis[parentA] != unis[parentB]]
bisDF = bis[bis[parenta] != bis[parentB]]
## How many SNPs are invariant between our two partens?
print "Number of sites with parent's calls invariant Uni: ",unis.shape[0]-unisDF.shape[0]
print "Number of SNPs after filter Uni:",unisDF.shape[0]
print
print "Number of sites with parent's calls invariant bi: ",bis.shape[0]-bisDF.shape[0]
print "Number of SNPs after filter bi:",bisDF.shape[0]
## Reset SNP sites where XL280a or XL280alpha is 1.
## We need to "swap" the genotype states where XL280 is 1, sperately for uni- and bi- sexual mating. 
## For uni-
## rewrite the genotype fields reference and alternate alleles for uni
unisDF.ix[unisDF.XL280==1,['Alt','Ref']] = unisDF.ix[unisDF.XL280==1,['Ref','Alt']] 
## For uni-
## Alter the genotype states.
unisDF.ix[unisDF.XL280 == 1,unisDF.columns[4:]] = abs(unisDF.ix[unisDF.XL280 == 1,unisDF.columns[4:]]-1)
## for bi-
## rewrite the genotype fields reference and alternate alleles for bi
bisDF.ix[bisDF.XL280a==1,['Alt','Ref']] = bisDF.ix[bisDF.XL280a==1,['Ref','Alt']] 
bisDF.ix[bisDF.XL280a == 1,bisDF.columns[4:]] = abs(bisDF.ix[bisDF.XL280a == 1,bisDF.columns[4:]]-1);
## We need to "swap" the genotype states where SS-A837 is 0, sperately for uni- and bi- sexual mating. 
## For uni-
## rewrite the genotype fields reference and alternate alleles for uni
unisDF.ix[unisDF[parentB]==0,['Alt','Ref']] = unisDF.ix[unisDF[parentB]==0,['Ref','Alt']];
## For uni-
## Alter the genotype states.
unisDF.ix[unisDF[parentB] == 0,unisDF.columns[4:]] = abs(unisDF.ix[unisDF[parentB] == 0,unisDF.columns[4:]]-1);
## for bi-
## rewrite the genotype fields reference and alternate alleles for bi
bisDF.ix[bisDF[parentB]==0,['Alt','Ref']] = bisDF.ix[bisDF[parentB]==0,['Ref','Alt']];
bisDF.ix[bisDF[parentB] == 0,bisDF.columns[4:]] = abs(bisDF.ix[bisDF[parentB] == 0,bisDF.columns[4:]]-1);
## See final shape of loci data
print "Shape of Uni DF: ", unisDF.shape
print "Shape of Bi DF", bisDF.shape
## Combine and Save Filtered genotype data
## Combine the genotype info
filtgeno = bisDF.combine_first(unisDF)
print "Final DataFrame shape:", filtgeno.shape 
## Save
def appendsave(DF,file,args):
    f = open(file,'w');
    #s = '#test\n'
    s = '# '+ today + ', ' + 'CSV: ' + str(args.inputf) + ', ' + 'Strainfile: ' + str(args.strainfile) + ', ' +  'Parents: ' + str(args.parents) + ', ' +'Donor: ' + str(args.donor) + ', ' + 'Mat Locus Pos: ' + str(args.MatPos) + ', '  + '\n'; 
    f.write(s)
    f.close()
    with open(file,'a') as f:
        DF.to_csv(f,index=False)
    f.close()
    pass 
appendsave(filtgeno,args.outputf,args)
#filtgeno.to_csv('../Crypto-FilterGenotypes-filter-Feb2016.csv', index=False)