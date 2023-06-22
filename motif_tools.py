#!/usr/bin/env python

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import pandas as pd
import random
import csv
from scipy import stats
from scipy.stats import ttest_1samp
import logging
import sys
import python_codon_tables as pct
import numpy as np
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS ##########
parser = argparse.ArgumentParser(
    description='This package contains tools for testing motif enrichment and avoidance hypotheses.\n'
    'Input sequence(s) are shuffled n times (random or syn codon with usage bias),\n'
    'motifs are counted, and 1 sample t-test is performed.\n'
    'Individual tools can also be used stand alone for other purposes\n'
    'Single or mutli fasta files are accepted.\n'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-mod', '--mode', choices=['1','2','3','4','5','6','7'], 
                    help='1-random shuffle\n'
                    '2-synonomous codon shuffle\n'
                    '3-find and count motifs\n'
                    '4-map motifs on plus and minus strands \n'
                    '5-one sample t-test\n'
                    '6-complete analysis using random shuffle\n'
                    '7-complete analysis using synonomous codon shuffle\n'
                    ' \n')

parser.add_argument('-f', '--file', help='path to your fasta file, or first counts file for ttest\n'
                    ' \n')

parser.add_argument('-f2', '--file2', help='only required with tool 4; second motif counts file for stand alone t-test\n'
                    ' \n')

parser.add_argument('-n', '--number', type=int, default=1, help='how many different times you want to shuffle your sequence(s)\n'
                    ' \n')

parser.add_argument('-motif', '--motif_list', help='CSV file containing a single column with motifs.\n'
                     'For ambiguous bases regex format should be used\n'
                     ' \n')

parser.add_argument('-tid', '--ncbi_tax_id', help='NCBI taxon id. Required. Default = 139 (Borrelia burgdorferi)\n'
                    ' \n', default = 139, type = int)

parser.add_argument('-out', '--outfiles', 
                    help='creates an output file for shuffled sequences, counts, and t-test\n'
                    'Enter a pre-fix to use before generic file names\n'
                    ' \n', default=None)

args = parser.parse_args()


########## FUNCTIONS ##########
def seq_shuffler(file, n):
  seqList = list(SeqIO.parse(file, "fasta"))
  for i in range(len(seqList)):
    #k = 1
    for z in range(n):
      seq=list(seqList[i].seq)
      random.shuffle(seq)
      seq = ''.join(seq)
      #output.append({'id': seqList[i].id+'_shuffled_'+str(k), 'seq': seq}) #use for unique sequence names
      output.append({'id': seqList[i].id, 'seq': seq})
      #k = k+1

def motif_counter(motif, file):
    seqList = list(SeqIO.parse(file, "fasta"))
    for i in range(0, len(seqList)):
        cds = str(seqList[i].seq) #converting the seq to a string for use for regex
        rvc = str(seqList[i].seq.reverse_complement()) #making a string of the reverse complement of the seq
        foundCDS = re.findall(pattern = r'(?=(' + motif + '))', string = cds, flags = re.IGNORECASE) #searchng for a specified motif
        foundRVC = re.findall(pattern = r'(?=(' + motif + '))', string = rvc, flags = re.IGNORECASE) #searchng for a specified motif
        output.append({'id': seqList[i].id, 'seq_length': len(seqList[i].seq), 'motif': motif, 
                   'fwd_match':len(foundCDS), 
                   'rev_match': len(foundRVC),
                  'total_match': len(foundCDS) + len(foundRVC)})

def motif_mapper(motif, file):
    seqList = list(SeqIO.parse(file, "fasta"))
    for i in range(0, len(seqList)):
      plus_seq = str(seqList[i].seq)
      minus_seq = str(seqList[i].seq.reverse_complement()) 
      for matches in re.finditer(pattern = r'(?=(' + motif + '))', string = plus_seq, flags = re.IGNORECASE): 
        output.append({'id': seqList[i].id, 'seq_length': len(seqList[i].seq), 'motif': motif, 'strand': 'plus','pos_start': matches.start()})
      for matches in re.finditer(pattern = r'(?=(' + motif + '))', string = minus_seq, flags = re.IGNORECASE): 
        output.append({'id': seqList[i].id, 'seq_length': len(seqList[i].seq), 'motif': motif, 'strand': 'minus','pos_start': matches.start()})

def motif_count_ttest(f1, f2):
    df = pd.read_table(f1, sep="\t", header=None)
    df2 = pd.read_table(f2, sep="\t", header=None)
    tstats = []
    for motif in df[2]:
        group1 = df[df[2] == motif]
        group2 = df2[df2[2] == motif]
        x = group1[5].mean()
        y = group2[5].mean()
        #print(motif)
        #print(stats.ttest_1samp(group2[5],x))
        a = stats.ttest_1samp(group2[5],x)
        tstats.append({'motif': motif, 'test_stat': a[0], 'p_val': a[1], 'observed': x, 'shuffled_mean': y,'log2fc': np.around(np.log2(x/y),3)})
    for item in tstats:
        print(item['motif'],"\t",item['test_stat'],"\t",item['p_val'], item['observed'], item['shuffled_mean'], item['log2fc'])

def complement(seq):
    complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    bases = list(seq.lower()) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def draw_a_syn_codon(codon):
    aa = codon.translate() # aa residue as a seqRecord
    syn_codons = cd_table[str(aa.seq)]
    probs = list(syn_codons.values())
    #print(probs, end = ";")
    probs[-1] = 1 - sum(probs[0:-1]) # force sum to 1
    #print(probs)
    new_codon = np.random.choice(list(syn_codons.keys()), size = 1, replace = False, p = probs)
    return new_codon[0]
    #print(syn_codons)

def syn_codon_shuffle(fasta_file, n):
    for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
        aa_seq = seq_rec.translate()
        #if re.match(r'\*', str(aa_seq.seq)):
            #print(seq_rec.id, ": internal stop found. Skip")
            #continue
        for i in range(n):
            new_str = ''
            for j in range(0, len(seq_rec)-2, 3):
                codon = seq_rec[j:(j+3)]
                new_str += draw_a_syn_codon(codon)
            output.append({'id': seq_rec.id, 'seq': new_str})
        #for item in output:
            #print('>',item['id'],'\n',item['seq'],sep='')

########## Modules ##########

########## random shuffle ##########
if args.mode == '1':
    
    output = []
    seq_shuffler(file = args.file, n = args.number)
    for item in output:
        print('>',item['id'],'\n',item['seq'],sep='')


########## syn codon shuffle ##########
if args.mode == '2':

    cd_table = pct.get_codons_table(args.ncbi_tax_id)
    output = []
    syn_codon_shuffle(args.file, args.number)
    for item in output:
        print('>',item['id'],'\n',item['seq'],sep='')


########## count(3) or map(4) motifs ##########
if args.mode == '3' or args.mode == '4':

    motif_list = [] 
    with open(args.motif_list) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        motif_list_file = csv.reader(csv_file)
        for row in motif_list_file:
            motif_list.append(row[0])
    
    if args.mode == '3':

        output = []
        for motif in motif_list:
            motif_counter(motif, args.file)
        #print('id','\t','seq_len','\t','motif','\t','plus_count','\t','minus_count','\t','total_count')
        for item in output:
            print(item['id'],"\t",item['seq_length'],"\t",item['motif'],"\t",item['fwd_match'],"\t",item['rev_match'],"\t",item['total_match'])
    
    if args.mode == '4':

        output = []
        for motif in motif_list:
            motif_mapper(motif, args.file)
        #print('id','\t','seq_len','\t','motif','\t','strand','\t','pos')     
        for item in output:
            print(item['id'],"\t",item['seq_length'],"\t",item['motif'],'\t',item['strand'],"\t",item['pos_start'])


########## t-test ##########
if args.mode == '5':
    motif_count_ttest(args.file, args.file2)


########## Complete; random shuffle - count motifs - ttest ##########
if args.mode == '6':
    motif_list = [] 
    with open(args.motif_list) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        motif_list_file = csv.reader(csv_file)
        for row in motif_list_file:
            motif_list.append(row[0])

    output = []
    for motif in motif_list:
        motif_counter(motif, args.file)
    
    if args.outfiles != None:
        for item in output:
            print(item['id'],"\t",item['seq_length'],"\t",item['motif'],"\t",item['fwd_match'],"\t",item['rev_match'],"\t",item['total_match'], file= open(args.outfiles+'_wt_motif_counts.tsv','a'))
  
    df = pd.DataFrame(output)
    #print(df)

    output = []
    seq_shuffler(args.file, args.number)
    seqList = output
    
    if args.outfiles != None:
        for item in output:
            print('>',item['id'],'\n',item['seq'],sep='',file = open(args.outfiles+'_shuffled.fasta', "a"))

    output = []
    for motif in motif_list:
        for item in seqList:
            cds = item['seq'] #converting the seq to a string for use for regex
            rvc = reverse_complement(item['seq']) #making a string of the reverse complement of the seq
            foundCDS = re.findall(pattern = r'(?=(' + motif + '))', string = item['seq'], flags = re.IGNORECASE) #searchng for a specified motif
            foundRVC = re.findall(pattern = r'(?=(' + motif + '))', string = reverse_complement(item['seq']), flags = re.IGNORECASE) #searchng for a specified motif
            output.append({'id': item['id'], 'seq_length': len(item['seq']), 'motif': motif, 'fwd_match':len(foundCDS), 'rev_match': len(foundRVC), 'total_match': len(foundCDS) + len(foundRVC)})
    
    df2 = pd.DataFrame(output)

    if args.outfiles != None:
        for item in output:
            print(item['id'],"\t",item['seq_length'],"\t",item['motif'],"\t",item['fwd_match'],"\t",item['rev_match'],"\t",item['total_match'], file= open(args.outfiles+'_shuffled_motif_counts.tsv','a'))
    
    tstats = []
    for seq_id in df.id.unique():
        df1a = df[df['id'] == seq_id]
        df2a = df2[df2['id'] == seq_id]
        for motif in df1a['motif']:
            df1b = df1a[df1a['motif'] == motif]
            df2b = df2a[df2a['motif'] == motif]
            x = df1b['total_match'].mean()
            y = df2b['total_match'].mean()
            a = stats.ttest_1samp(df2b['total_match'],x)
            tstats.append({'id': seq_id,'motif': motif, 'test_stat': a[0], 'p_val': a[1], 'observed': x, 'shuffled_mean': y,'log2fc': np.around(np.log2(x/y),3)})

    for item in tstats:
        print(item['id'],'\t',item['motif'],"\t",item['test_stat'],"\t",item['p_val'], item['observed'], item['shuffled_mean'], item['log2fc'])
    
    if args.outfiles != None:
        for item in tstats:
            print(item['id'],'\t',item['motif'],"\t",item['test_stat'],"\t",item['p_val'], item['observed'], item['shuffled_mean'], item['log2fc'],file= open(args.outfiles+'_ttest.tsv','a'))


########## Complete; syn codon shuffle - count motifs - ttest ##########
if args.mode == '7':
    motif_list = [] 
    with open(args.motif_list) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        motif_list_file = csv.reader(csv_file)
        for row in motif_list_file:
            motif_list.append(row[0])

    output = []
    for motif in motif_list:
        motif_counter(motif, args.file)
    
    df = pd.DataFrame(output)
    #print(df)

    cd_table = pct.get_codons_table(args.ncbi_tax_id)
    output = []
    for seq_rec in SeqIO.parse(args.file, 'fasta'):
        aa_seq = seq_rec.translate()
        if re.match(r'\*', str(aa_seq.seq)):
            print(seq_rec.id, ": internal stop found. Skip")
            continue
        for i in range(args.number):
            new_str = ''
            for j in range(0, len(seq_rec)-2, 3):
                codon = seq_rec[j:(j+3)] # codon as a seqRecord
                #if str(codon.translate()) not in cd_table:
                    #continue
                new_str += draw_a_syn_codon(codon)
            output.append({'id': seq_rec.id, 'seq': new_str})
            #print(new_str)   
    
    seqList = output

    output = []
    for motif in motif_list:
        for item in seqList:
            cds = item['seq'] #converting the seq to a string for use for regex
            rvc = reverse_complement(item['seq']) #making a string of the reverse complement of the seq
            foundCDS = re.findall(pattern = r'(?=(' + motif + '))', string = item['seq'], flags = re.IGNORECASE) #searchng for a specified motif
            foundRVC = re.findall(pattern = r'(?=(' + motif + '))', string = reverse_complement(item['seq']), flags = re.IGNORECASE) #searchng for a specified motif
            output.append({'id': item['id'], 'seq_length': len(item['seq']), 'motif': motif, 'fwd_match':len(foundCDS), 'rev_match': len(foundRVC), 'total_match': len(foundCDS) + len(foundRVC)})
    
    df2 = pd.DataFrame(output)
    #print(df2)

    tstats = []
    for seq_id in df.id.unique():
        df1a = df[df['id'] == seq_id]
        df2a = df2[df2['id'] == seq_id]
        for motif in df1a['motif']:
            df1b = df1a[df1a['motif'] == motif]
            df2b = df2a[df2a['motif'] == motif]
            x = df1b['total_match'].mean()
            y = df2b['total_match'].mean()
            a = stats.ttest_1samp(df2b['total_match'],x)
            tstats.append({'id': seq_id,'motif': motif, 'test_stat': a[0], 'p_val': a[1], 'observed': x, 'shuffled_mean': y,'log2fc': np.around(np.log2(x/y),3)})

    for item in tstats:
        print(item['id'],'\t',item['motif'],"\t",item['test_stat'],"\t",item['p_val'], item['observed'], item['shuffled_mean'], item['log2fc'])
