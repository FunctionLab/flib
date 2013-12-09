#!/usr/bin/env python
#from __future__ import division
import os
import sys
import argparse
import re
#import random
from collections import defaultdict

parser = argparse.ArgumentParser(description='INSERT DESCRIPTION')
parser.add_argument(
    '-g',
    dest='genemap',
    nargs='?',
    type=argparse.FileType('r'),
    help='Genemap File from OMIM')
parser.add_argument(
    '-m',
    dest='mim2gene',
    nargs='?',
    type=argparse.FileType('r'),
    help='Mim2Gene File from OMIM')
parser.add_argument(
    '-d',
    dest='dofile',
    nargs='?',
    type=argparse.FileType('r'),
    help='Doid OBO File')
parser.add_argument(
    '-u',
    dest='umls_mode',
    action='store_true',
    default=False,
    help='Output UML terms')
parser.add_argument(
    '-o',
    dest='out',
    nargs='?',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='Mapping output file')
arg = parser.parse_args()

"""
Idea: Use the mim2gene and genemap files from omim to extract a mapping of omim disease phentype to its annotated genes

Logic:
1) Parse the mim2gene file, identifying a mapping between mim number and gene id for type: gene and gene/phenotype, this excluded moved/removed and phenotype types
2) Specify limitations necessary for parsing the genemap file, so far I think: limit Phenotype mapping method to 3 and status codes for both confirmed and provisional
3) Parse the genemap file, line by line if the status [col 7] is c or p then proceed to step 4, store this info for the output
3a) Check if a geneid exists, else skip
4) Separate the disorders (col 14) by semicolon ;, identify their phenotype mapping, classification and omim id
4a) Ignore the disorders with [non disease] and ?unconfirmed?, note if {suscpetibility} for output, also ignore any entried with phenotype other than (3)
"""

class mim_disease:
    def __init__(self):
        self.mimid = ''
        self.is_susceptibility = 0 #Whether it has {}
        self.phe_mm = ''  #Phenotype mapping method
        self.genetuples = [] #(Gene ID, Gene Status)


limit_type = ['gene','gene/phenotype']
limit_pheno = '(3)'
limit_status = ['C','P']

#This should be standardized, but you never know
#Most disorders that have omimphenotypes fit this expression
find_mimid = re.compile('\, [0-9]* \([1-4]\)')

mim_gene_map = {}
mimdiseases = {}

mim_gene_list = arg.mim2gene.read().splitlines()
genemap_list = arg.genemap.read().splitlines()

doid_omim = {}
doid_umls = defaultdict(set)

obo_lines = arg.dofile.read().splitlines()[::-1]
while obo_lines:
     l = obo_lines.pop()
     if l == '[Term]':
         while l != '':
             l = obo_lines.pop()
             if l.startswith('id:'):
                 doid = re.search('DOID:[0-9]+',l).group(0)
             if l.startswith('xref: OMIM:'): 
                 omim = re.search('[0-9]+',l).group(0)
                 if not doid_omim.has_key(doid): doid_omim[doid] = []
                 if omim not in doid_omim[doid]: doid_omim[doid].append(omim)
             if l.startswith('xref: UMLS_CUI:'):
                 umls = re.search('C[0-9]+',l).group(0)
                 doid_umls[doid].add(umls)

#Get mim to gid mapping
for l in mim_gene_list:
    l_split = l.split('\t')
    mim = l_split[0]
    type = l_split[1]
    gid = l_split[2]

    if type in limit_type:
        if mim_gene_map.has_key(mim): print mim, 'key exists'
        mim_gene_map[mim] = gid

if arg.umls_mode:
    for doid,umlsids in doid_umls.iteritems(): 
        omim = set()
        if doid in doid_omim:
            omim = doid_omim[doid]

        for mimid in omim:
            if mimid not in mimdiseases:
                mimdiseases[mimid] = mim_disease()
                mimdiseases[mimid].mimid = mimid
            for umlsid in umlsids:
                mimdiseases[mimid].genetuples.append((umlsid,'NA'))


else:
    for l in genemap_list:
        #The choice of fields relies on info from the genemap.key file from omim
        l_split = l.split('|')
        status = l_split[6].strip()
        mim_gene = l_split[9].strip()
        disorders = l_split[13].strip()

        #continuation of disorder field
        d2 = l_split[14].strip()
        d3 = l_split[15].strip()
        if d2 != '': disorders= disorders + ' ' + d2
        if d3 != '': disorders= disorders + ' ' + d3
        
        if disorders != '' and status in limit_status and mim_gene_map.has_key(mim_gene):
            #print 'Status ok, not blank and genemap has key'
            
            geneid = mim_gene_map[mim_gene]
            tuple_gid_status = (geneid, status)
            
            disorders_list = disorders.split(';')
            for d in disorders_list:
                if '[' not in d and '?' not in d:
                    mim_info = re.search(find_mimid,d)
                    if mim_info:
                        #print 'Has necessary info'
                        #TODO: Make sure to include ? and [
                        info_split = mim_info.group(0).split(' ')
                        mim_disease_id = info_split[1].strip()
                        mim_phetype = info_split[2].strip()
                        if mim_phetype == limit_pheno:
                            #print 'Correct phenotype'
                            if not mimdiseases.has_key(mim_disease_id):
                                mimdiseases[mim_disease_id] = mim_disease()
                                mimdiseases[mim_disease_id].mimid = mim_disease_id
                                mimdiseases[mim_disease_id].phe_mm = mim_phetype
                            if '{' in d: mimdiseases[mim_disease_id].is_susceptibility = 1
                            if tuple_gid_status not in mimdiseases[mim_disease_id].genetuples:
                                mimdiseases[mim_disease_id].genetuples.append(tuple_gid_status)

out = arg.out
out.write('!Mapped file of Omim to Disease Ontology in GO gene association format')
out.write('\n!DB\tEntrez Gene ID\tIgnore\tIgnore\tDOID\tOMIMID\tPhenoMappingMethod-GeneStatus-Disease(D)orSusceptibilityto(S)\t7-17 are just blank tabs')

for i in doid_omim:
    doid = i
    omim_list = doid_omim[doid]
    for o in omim_list:
        omim_id = o
        if mimdiseases.has_key(omim_id):
            mim_entry = mimdiseases[omim_id]
            if mim_entry.is_susceptibility: d_or_s = 'S'
            else: d_or_s = 'D'
            
            for g in mim_entry.genetuples:
                out.write('\nOMIM-Entrez\t'+\
                        g[0] + '\t \t \t' +\
                        doid + '\t' +\
                        o + '\t' +\
                        mim_entry.phe_mm + '-' + g[1] + '-' + d_or_s + '\t' +\
                        '\t \t \t \t \t \t \t \t \t \t ')
