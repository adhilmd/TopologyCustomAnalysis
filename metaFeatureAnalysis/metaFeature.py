#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:20:41 2019

@author: amohamme
"""

#python3 metaFeature.py -bdg "/lustre/home/amohamme/data/Chip/NegSupercoiling/R1323/Outputf/DataProcessing/bed_graphfiles/proc_R1323_Top1_S-28_IP-A.bdg" -bed "/lustre/home/amohamme/data/Chip/NegSupercoiling/R1323/Outputf/DataProcessing/bed_graphfiles/proc_R1323_Top1_S-28_IP-A.bed" -ft "/lustre/home/amohamme/data/AnnotationFiles/saccer2011/Sc_proteincoding_header.txt" -st "yes" -stn "strand" -bs "1000" -fl "yes" -flb "500" -sm "Top1Neg" -t "16" -o "/storage/home/amohamme/My_Tools/metaFeature/"

import argparse

def argparser():
    parser = argparse.ArgumentParser(prog = '''\nmetaFeature.py''', description='''\n----------Meta Feature analysis using bedgraph and bedfile--------- \n
    \n[Date: 7th March 2018], \n[help: python metaFeature.py -h]\n''', usage = 'metaFeature.py *args')
    parser.add_argument('-bdg','--bdg', type=str, dest='bdgfile', help="bedGraph file with no header (Mandatory)", action = 'store', required = True)
    parser.add_argument('-bed','--bed', type=str, dest='bedfile', help="bedfile with header (Mandatory)", action = 'store', required = True)
    parser.add_argument('-ft','--feature', type=str, dest='featurefile', help="Feature bed file with header, four columns are mandatory with the column names (chrom,start,end,name) (Mandatory)", action = 'store', required = True)
    parser.add_argument('-st','--strand', type=str, dest='strand', help="Value:[yes or no] If yes, normalize the strand orientaion for aggregation (Default = no)", action = 'store', default="no")
    parser.add_argument('-stn','--strandname', type=str, dest='strname', help="Column name for strand from feature file, column should contain (+ or 1, - or -1)  (Default=None)", action = 'store', default=None)
    parser.add_argument('-bs','--basesmerging', type=int, dest='basesmerging', help="Total base for feature aggragation (Mandatory)", action = 'store', required = True)
    parser.add_argument('-fl','--flanking', type=str, dest='flanking', help="Value:[yes or no] If yes, flanking region will be considered (Default = no)", action = 'store', default ="no")
    parser.add_argument('-flb','--flankingbases', type=int, dest='flankingbases', help="Total flanking bases. (Required only if -fl is yes)", action = 'store', default = 0)
    parser.add_argument('-sm','--samplename', type=str, dest='samplename', help="Sample name", action = 'store', required = True)
    parser.add_argument('-sb','--subtype', type=str, dest='subtype', help="Subtype column name from feature file, where subtype is used for aggregation (Default=None)", action = 'store', default=None)
    parser.add_argument('-t','--multithreads', type=int, dest='procs', help="Number of threads (Default=2)", action = 'store', default=2)
    parser.add_argument('-o','--outdir', type=str, dest='outdir', help="outdir (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

#required library
from pybedtools import BedTool  
import pandas as pd
import os
import time
from multiprocessing import Pool
from functools import partial
import warnings
import shutil
warnings.filterwarnings("ignore")       


def strandflip(mrt,basesmerging):
    strp = mrt[mrt['strand'] == 1]
    strn = mrt[mrt['strand'] == -1]
    if (strn.shape[0] >= 1):
        strn.is_copy = False
        strn['value'] = list(map(lambda x: (x-(basesmerging+1))*-1, strn['value'].tolist()))
    fn = pd.concat([strp,strn])
    return(fn)

def flankingflip(fn,flankingbases,basesmerging,strand):
    if (strand == 'yes'):
        fn1 = fn[fn['name_x'].str.contains("_1$")]
        fn2 = fn[fn['name_x'].str.contains("_2$")]
        fn1p = fn1[fn1['strand'] == 1]
        fn2n = fn2[fn2['strand'] == -1]
        fnup = pd.concat([fn1p,fn2n])
        fnup['value'] = fnup['value'] - (flankingbases)
        fn1n = fn1[fn1['strand'] == -1]
        fn2p = fn2[fn2['strand'] == 1]
        ffup = pd.concat([fn1n,fn2p])
        ffup['value'] = ffup['value']  + (basesmerging)
        fnb = pd.concat([fnup,ffup])
    else:
        fn1 = fn[fn['name_x'].str.contains("_1$")] 
        fn2 = fn[fn['name_x'].str.contains("_2$")] 
        fn1.is_copy = False  
        fn2.is_copy = False
        fn1['value'] = fn1['value'] - (flankingbases)
        fn2['value'] = fn2['value'] + (basesmerging)
        fnb = pd.concat([fn1,fn2])
    return(fnb)

def fcontinous(mdf,featurechunk,basesmerging):
    ge = mdf['name'].unique().tolist()
    lzero = len(str(basesmerging))
    gdict=dict()
    for i in ge:
        gdict[i] = list(range(1,basesmerging+1))
    gd = pd.melt(pd.DataFrame(gdict))
    fcnt = featurechunk[['name','cnt']]
    gd=gd.merge(fcnt,left_on='variable',right_on='name')
    gd['value']=gd['value'].map(str)
    gd['cnt']=gd['cnt'].map(str)
    gd['ID'] = gd['value'].apply(lambda x: x.zfill(lzero))
    gd['ID'] = gd['cnt'] +'.'+ gd['ID']
    gd['ID'] = gd['ID'].map(float)
    mdf['normb']=mdf['normb'].map(str)
    mdf['cnt']=mdf['cnt'].map(str)
    mdf['ID'] = mdf['normb'].apply(lambda x: x.zfill(lzero))
    mdf['ID'] = mdf['cnt'] +'.'+ mdf['ID']
    mdf['ID'] = mdf['ID'].map(float)
    gd = gd.sort_values(by=['ID'])
    mdf = mdf.sort_values(by=['ID'])
    ftwmerge = pd.merge_asof(gd, mdf, on='ID', direction='nearest')
    ftwmerge = ftwmerge[['name_x','ID','scorer','strand','value','Dens']]
    ftwmerge['value'] = ftwmerge['value'].map(int)
    return(ftwmerge)

def postprocess(mdf,featurechunk,basesmerging):
    mdf.columns = ['chromr','startr','endr','scorer'] + featurechunk.columns.tolist() + ['overlap'] + ['Dens']
    mdf.loc[mdf['Dens'] > 0, 'Dens'] = 1
    mdf['normb'] = ((mdf['startr']-mdf['start'])/(mdf['end']-mdf['start']))*basesmerging
    mdf.loc[mdf['normb'] < 1, 'normb'] = 1
    mdf = mdf[~mdf['normb'].isnull()]
    mdf['normb']=mdf['normb'].round().map(int)
    mdf=mdf[['name','scorer','strand','cnt','normb','Dens']]
    mdfnew = mdf.groupby(['cnt','normb'], as_index=False).agg({"scorer": "median", "Dens":"mean"})
    mdfnew['Dens'] = mdfnew['Dens'] + 0.01
    mdfnew['Dens'] = mdfnew['Dens'].round().map(int)
    mdf = mdf.drop(['scorer','Dens'], axis=1)
    mdf = mdf.drop_duplicates()
    mdf = mdf.merge(mdfnew,on=['cnt','normb'])
    mdf=mdf[['name','scorer','strand','cnt','normb','Dens']]
    mdf = mdf.fillna(0)
    mdf['cnt'] = mdf['cnt'].map(int)
    mdf['normb'] = mdf['normb'].map(int)
    mdf['Dens'] = mdf['Dens'].map(int)
    mdf = mdf.sort_values(by=['cnt'])
    return(mdf)

def featuremerge(featurechunk,bdgfile,beddf,basesmerging,strand,flanking,flankingbases):
    y = BedTool.from_dataframe(featurechunk)
    x = BedTool(bdgfile)
    mp = x.intersect(y,wo=True)
    z = BedTool.from_dataframe(beddf)
    mpz = mp.intersect(z,wao=True)
    mdf = pd.read_table(mpz.fn)
    mdf = mdf.ix[:,list(range(0,14))+[17]]
    if (flanking == 'yes'):
        fnb = postprocess(mdf,featurechunk,flankingbases)
        fnb = fcontinous(fnb,featurechunk,flankingbases)
        if (strand == 'yes'):
            fnb = strandflip(fnb,flankingbases)
        fnb = flankingflip(fnb,flankingbases,basesmerging,strand) 
    else:
        fnb = postprocess(mdf,featurechunk,basesmerging)
        fnb = fcontinous(fnb,featurechunk,basesmerging)
        if (strand == 'yes'):
            fnb = strandflip(fnb,basesmerging) 
    return(fnb)

def main(featuredf,procs,bdgfile,bedfile,strand,basesmerging,flanking,flankingbases,samplename,sitems):
    print (os.getcwd())
    progstarts = time.time()    
    featuredf = featuredf.sort_values(by=['chrom','start'])
    featuredf = featuredf.fillna("PT")
    c = featuredf.groupby(["name"]).cumcount()
    c = c.replace(0, '').astype(str)
    featuredf['name'] += c
    featuredf['cnt'] = list(range(featuredf.shape[0]))
    nsplit = int(len(featuredf)/procs)
    featuresplit = [featuredf.iloc[i*nsplit:(i+1)*nsplit].copy() if i != procs-1 else featuredf.iloc[i*nsplit:featuredf.shape[0]].copy() for i in range(procs)]
    beddf = pd.read_csv(bedfile, sep="\t", header=0)
    beddf = beddf.iloc[:,0:3]
    beddf.columns = ['chromb','startb','endb']
    beddf = beddf.sort_values(by=['chromb','startb'])
    pool = Pool(processes=procs)
    mrt=pool.map(partial(featuremerge,bdgfile=bdgfile,beddf=beddf,basesmerging=basesmerging,strand=strand,flanking='no',flankingbases=0), featuresplit)
    pool.close()
    pool.join()
    bfn = pd.concat(mrt) 
    del mrt,featuresplit
    if (flanking == 'yes'):
        featuredf['ustart'] = featuredf['start']-flankingbases
        featuredf['dend'] = featuredf['end']+flankingbases
        featuredf = featuredf[['chrom','ustart','dend']+featuredf.columns[1:len(featuredf.columns)-2].tolist()]
        featuredf.loc[featuredf['ustart'] < 0, 'ustart'] = 0
        up = featuredf[['chrom','ustart','start']+featuredf.columns[5:len(featuredf.columns)].tolist()]
        down = featuredf[['chrom','end','dend']+featuredf.columns[5:len(featuredf.columns)].tolist()]
        up.is_copy = False
        down.is_copy = False
        up['name'] = [x+'_1' for x in up['name'].tolist()]
        up = up.rename(index=str, columns={"ustart": "start", "start": "end"})
        down['name'] = [x+'_2' for x in down['name'].tolist()]
        down = down.rename(index=str, columns={"end": "start", "dend": "end"})
        featuredfflank = pd.concat([up,down])
        featuredfflank.is_copy = False
        featuredfflank['cnt'] = list(range(featuredfflank.shape[0]))
        nsplit = int(len(featuredfflank)/procs)
        featuresplitflank = [featuredfflank.iloc[i*nsplit:(i+1)*nsplit].copy() if i != procs-1 else featuredfflank.iloc[i*nsplit:featuredfflank.shape[0]].copy() for i in range(procs)]
        pool = Pool(processes=procs)
        mrtf=pool.map(partial(featuremerge,bdgfile=bdgfile,beddf=beddf,basesmerging=basesmerging,strand=strand,flanking=flanking,flankingbases=flankingbases), featuresplitflank)
        pool.close()
        pool.join()    
        fnb = pd.concat(mrtf) 
        findat = pd.concat([bfn,fnb])
        del bfn,fnb,mrtf
    else:
        findat = bfn
        del bfn
    aggfn = findat.groupby('value', as_index=False).agg({"scorer": "median", "Dens":"sum"})
    aggfn['samplename'] = samplename
    aggfn = aggfn[["value","scorer","Dens","samplename"]]
    nowtime = time.time()
    runtime = nowtime-progstarts
    return([aggfn,runtime])

if __name__ == "__main__":
    args = argparser()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    tempdir = args.outdir+"/"+args.samplename+"_temp/"
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)
    os.chdir(tempdir)
    if (args.flanking == 'yes' and args.flankingbases == 0):
        print ("Error: Please provide flanking bases (-flb) value or make -fl no and execute the program again")
        quit()
    if (args.strand == 'yes' and args.strname == None):
        print ("Error: Please provide a column name from feature file for strand information (-stn) or make -st no and execute the program again")
        quit()
    featuredf = pd.read_csv(args.featurefile, sep="\t", header=0)
    if (args.strand == "yes"):
        featuredf = featuredf.rename(index=str, columns={args.strname: "strand"})
        featuredf['strand'] = featuredf['strand'].fillna("1")
        featuredf['strand'] = featuredf['strand'].map(str)
        featuredf[['strand']] = featuredf[['strand']].replace(["+","-"], ["1","-1"])
        featuredf['strand'] = featuredf['strand'].map(int)
    if (args.subtype != None):
        subtypes = featuredf[args.subtype].map(str)
        subtypes = featuredf[args.subtype].unique().tolist()
        subtlist = []
        for items in subtypes:
            itypes = featuredf[featuredf[args.subtype] == items]
            result = main(itypes,args.procs,args.bdgfile,args.bedfile,args.strand,args.basesmerging,args.flanking,args.flankingbases,args.samplename,items)
            result[0]['subtype'] = items
            subtlist.append(result[0])
        subtlist = pd.concat(subtlist)
        subtlist.to_csv(args.outdir+"/"+args.samplename+"_metaFeature.tsv",sep="\t",index=False)
    else:
        items = "NA"
        result = main(featuredf,args.procs,args.bdgfile,args.bedfile,args.strand,args.basesmerging,args.flanking,args.flankingbases,args.samplename,items)
        result[0]['subtype'] = args.samplename
        result[0].to_csv(args.outdir+"/"+args.samplename+"_metaFeature.tsv",sep="\t",index=False)
    shutil.rmtree(tempdir)
    print ("Done")
