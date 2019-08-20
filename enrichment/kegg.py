#!/usr/bin/python3

"""
Author : Sitao Zhu 
Date   : 2019-08-14
Version: 1th
Desc   : for RNAseq KEGG enrichment analysis
input file :

foreground(fg):
Pn21.1406       K21971|1|0.0|658|nnu:104597730|methyltransferase NSUN6 [EC:2.1.1.-]
Pn2.310 K21843|1|5e-20|97.8|dct:110108134|tetratricopeptide repeat protein 7
Pn2.398 K21843|1|5e-20|97.8|dct:110108134|tetratricopeptide repeat protein 7
----------

background(bg):
Pn1.1001        K15803|1|7e-46|166|nta:107812685|(-)-germacrene D synthase [EC:4.2.3.75]
Pn1.1002        K15803|1|3e-25|95.1|nta:107793646|(-)-germacrene D synthase [EC:4.2.3.75]
-----------

kegg fasta:
>ath:AT1G30550  S-adenosyl-L-methionine-dependent methyltransferase superfamily protein; K14292 trimethylguanosine synthase [EC:2.1.1.-]
MKKEAESLIEKEHGTNPKISRYWIQRYDLFSKYDQGIEMDEEGWYSVTPEEIAIKQAERC
RGKVVIDCFSGVGGNTIQFAKVCSSVIAIDIDPMKIALAMNNAKVYGVANRIDFVTGDFM
QLAPSLKGDVLFLSPPWGGPTYSKVESYKLDMLLPRDGYSLFQTALSITPNIIMFLPKNI
DLAQLEELACLSSPPLTLEIEENSIGGEIKAITAYFSSNVV
-----------

komap:
K00001  00010 00350 01100 01110 00071
K00002  00010 00040 01100 01110 00561
K00006  01110 00564
K00008  00040 01100 00051
K00010  01100 00562
K00011  00040 01100 00052 00051 00790 00561
------------

maptitle:
00010   Glycolysis / Gluconeogenesis    Metabolism      Carbohydrate metabolism
00020   Citrate cycle (TCA cycle)       Metabolism      Carbohydrate metabolism
00030   Pentose phosphate pathway       Metabolism      Carbohydrate metabolism
00040   Pentose and glucuronate interconversions        Metabolism      Carbohydrate metabolism
00051   Fructose and mannose metabolism Metabolism      Carbohydrate metabolism
00052   Galactose metabolism    Metabolism      Carbohydrate metabolism
------------

"""

import sys,re,os
import argparse
from collections import defaultdict
from scipy.stats import hypergeom
import rpy2.robjects as R
import numpy as np

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector



def get_species_ko(kegg_fasta,species):
    '''
    exact given specie ko from kegg fasta
    :param kegg_fasta: kegg fasta file
    :param species: three small-caption letters for specie
    :return: the given species ko
    '''
    with open(kegg_fasta,'r') as F:
        for line in F:
            pass

def koMap(komap):
    """
    read komap file into a dict
    :param komap: komap file  # K00001  00010 00350 01100 01110 00071
    :return: a dict
    """
    kos = {}
    with open(komap) as F:
        for line in F:
            fields = line.strip().split('\t') #K00001  00010 00350 01100 01110 00071
            kos[fields[0]] = fields[1].split(' ')
    return kos

def mapTitle(maptitle):
    """
    save maptitle into memory
    :param maptitle: a ko annotation file
    :return:
    """
    titles = defaultdict(dict)
    with open(maptitle) as F:
        for line in F:
            fields = line.strip().split('\t') # 00010   Glycolysis / Gluconeogenesis    Metabolism      Carbohydrate metabolism
            titles[fields[0]]['level3'] = fields[1]
            titles[fields[0]]['level2'] = fields[2]
            titles[fields[0]]['level1'] = fields[3]
    return titles

def adjust_pvalues(pvalues, method='BH'):
    #pvalue_lst = [v.r['p.value'] for v in pvalues]
    #return R['p.adjust'](R.FloatVector(pvalues), method=method)
    stats = importr('stats')
    p_adjust = stats.p_adjust(FloatVector(pvalues), method='BH')

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def hyper_geom(fg,bg,kos,titles):
    """
    deal with foreground gene ko
    :param fg: foreground
    :return:
    """

    genes1 = {}
    genes2 = {}

    sum1 = 0
    sum2 = 0

    stats = defaultdict(dict)

    with open(fg) as F:
        for line in F:
            line = line.strip()
            if line.startswith("#") or line == "":continue
            fields = line.split("\t") # \t
            if len(fields) == 0 or len(fields) == 1: continue
            fields[1] = re.sub('\|[^!]*!*','!',fields[1])
            fields[1] = re.sub('!$','',fields[1])
            tmp = fields[1].split('!')
            if not genes1.get(fields[0],0):
                sum1 += 1 # the total number of foreground gene
                genes1[fields[0]] = 1
            for KO in tmp:
                if kos.get(KO,0):
                    for ko in kos[KO]:
                        if titles.get(ko,0):
                            if not stats.get(ko,0):
                                stats[ko]['pathway'] = titles[ko]['level3']
                                stats[ko]['genes'] = []
                                stats[ko]['kos'] = []
                                if bg:
                                    stats[ko]['genes2'] = []
                                    stats[ko]['kos2'] = []
                            if 
                            stats[ko]['genes'].append(fields[0]) # gene
                            stats[ko]['kos'].append(KO) # ko
    if not bg:
        num = 0
    else:
        with open(bg,'r') as F:
            for line in F:
                line = line.strip()
                if line.startswith("#") or line == "": continue
                fields = line.split('\t')  # geneid K15803|1|7e-46|
                if len(fields) == 0 or len(fields) == 1 : continue
                fields[1] = re.sub('\|[^!]*!*', '!', fields[1])
                fields[1] = re.sub('!$', '', fields[1])
                tmp = fields[1].split('!')
                if not genes2.get(fields[0],0):
                    sum2 += 1   # the total number of background gene
                    genes2[fields[0]] = 1 # gene

                for KO in tmp: # bg(gene set) annotation ko list
                    if kos.get(KO, 0): # kos KEGG database KO data
                        for ko in kos[KO]:
                            if titles.get(ko, 0): # titles KO and it function annotation
                                if not stats.get(ko):continue
                                stats[ko]['genes2'].append(fields[0])
                                stats[ko]['kos2'].append(KO)
        num1 = 0
        num2 = 0
        pValues = []
        for ko in sorted(stats.keys()):
            num1 = len(stats[ko]['genes'])  # the gene number of given ko in fg
            num2 = len(stats[ko]['genes2']) # the gene number of given ko in bg
            # num1-1 :选出来特定ko的基因数 num2:特定ko的基因总数 sum2:基因数总数 sum1:选出来的基因数
            # phyper(num1 -1, num2,     sum2-num2, sum1) R version
            #  选出来白球的数目  白球的总数，黑球的总数，  选出来球的总数
            print(num1, sum2, num2, sum1)
            pValues.append(hypergeom.cdf(num1, sum2, num2, sum1))
        print(pValues)
        print(len(pValues))
        qValues = p_adjust_bh(pValues)
        print(qValues)

        i = 0
        for ko in sorted(stats.keys()):
            stats[ko]['pvalues'] = pValues[i]
            stats[ko]['qvalues'] = qValues[i]
            i += 1

        for ko in sorted(stats,key=lambda a:stats[a]['pvalues']):
            num1 = len(stats[ko]['genes'])
            num2 = len(stats[ko]['genes2'])
            content = stats[ko]['pvalues']
            #content = [stats[ko]['pathway'],num1,num2,stats[ko]['pvalues'],stats[ko]['qvalues'],ko,titles[ko]['level1'],titles[ko]['level2'],stats[ko]['genes'],stats[ko]['kos']]
            print(content)




def main():
    """
    %prog [-options]
    The program for enrichment analysis
    :return:
    """
    parser = argparse.ArgumentParser(prog='kegg.py')
    parser.add_argument('-fg', '--foreground', help='the selected gene ko (may from different expression gene)')
    parser.add_argument('-bg', '--background', help='the genome kegg annotation file(all ko)')
    parser.add_argument('-kegg','--kegg', help='kegg fasta')
    parser.add_argument('-komap', '--komap', help='ko_map.tab file')
    parser.add_argument('-maptitle','--maptitle', help='kegg map title file')
    parser.add_argument('-o','--output', help='outfile')

    args = parser.parse_args()
    print(args)
    fg = args.foreground
    bg = args.background
    kegg = args.kegg
    komap = args.komap
    maptitle = args.maptitle
    output = args.output



    kos = koMap(komap)
    titles = mapTitle(maptitle)
    hyper_geom(fg,bg,kos,titles)

if __name__ == "__main__":
    main()

























