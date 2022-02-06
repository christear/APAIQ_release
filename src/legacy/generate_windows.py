#!/usr/bin/env python
# coding: utf-8

import os,sys
import glob
from collections import OrderedDict,defaultdict
import pipes
import pprint
import tempfile
import re
import argparse
import gc
#from pyfaidx import Fasta

def get_genome_sequence(fa_file):
    f = open(fa_file,"r")
    line = f.readline()
    line = line.rstrip('\n')
    f.close()
    return line

def Get_block_position(root_dir,input_file,strand,window,block_length):
    chr_dict = dict()
    reference = dict()
    blocks = []
    pre_pos = 0
    nex_pos = 0
    start_pos = 0
    end_pos   = 0
    block_num = 0
    count = 0
    Chromosome = ''
    window /= 1.5
    window = int(window)
    with open(input_file, "r") as bg_file:
        lines = bg_file.readlines()
        for i,line in enumerate(lines):
            line = line.rstrip('\n')
            chromosome,start,end,val = line.split('\t')
            val = float(val)
            pos1 = int(start)+1
            pos2 = int(end)+1
            extend1 = window
            extend2 = window
            if(i>0 and i< len(lines)-1):
                pre_line = lines[i-1].rstrip('\n')
                nex_line = lines[i+1].rstrip('\n')
                pre_chr,_,pre_pos,_ = pre_line.split('\t')
                nex_chr,nex_pos,_,_ = nex_line.split('\t')
                pre_pos = int(pre_pos)+1
                nex_pos = int(nex_pos)+1
                if(pre_chr == chromosome and pos1-pre_pos<2*window):
                    if(pos1-pre_pos<window):
                        extend1 = 0
                    else:
                        extend1 = pos1-pre_pos-window
                if(nex_chr == chromosome and nex_pos-pos2<window):
                    extend2 = nex_pos-pos2

            if ('chr' not in chromosome):
                chromosome = 'chr'+chromosome

            if chromosome not in chr_dict.keys():
                end_pos = i
                if(i>0 and len(Chromosome)<=5 and 'Y' not in Chromosome and 'M' not in Chromosome):
                    blocks.append((Chromosome,strand,block_num,start_pos,end_pos))
                start_pos = i
                block_num = 0
                chr_dict[chromosome] = ''
                count =  pos2-pos1+extend1+extend2
            else:
                if(count>block_length and pos1-pre_pos>1000):
                    end_pos = i
                    #baseName = '%s_%s_%s'%(chromosome,strand,block_num)
                    if(len(Chromosome)<=5 and 'Y' not in Chromosome and 'M' not in Chromosome):
                        blocks.append((Chromosome,strand,block_num,start_pos,end_pos))
                        block_num += 1
                        start_pos = i
                        count = 0
                count += pos2-pos1+extend1+extend2
            Chromosome = chromosome
    #baseName = '%s_%s_%s'%(chromosome,strand,block_num)
    if(len(Chromosome)<=5 and 'Y' not in Chromosome and 'M' not in Chromosome):
        blocks.append((Chromosome,strand,block_num,start_pos,len(lines)))
    return blocks


def split_chr_bedGraph2(root_dir,input_file,chromosome,strand,window,reference,depth,start,end):
    window /= 1.5
    window = int(window)
    block = []
    pre_pos = 0
    nex_pos = 0
    print("Start generating block %d %d\n"%(start,end))
    #reference = get_genome_sequence('%s.%s.fa'%(fa_file,chromosome))
    with open(input_file, "r") as bg_file:
        lines = bg_file.readlines()
        for i,line in enumerate(lines):
            if(i<start or i>=end):
                continue
            line = line.rstrip('\n')
            chromosome,pos1,pos2,val = line.split('\t')
            if(len(chromosome)>5 or 'Y' in chromosome or 'M' in chromosome):
                continue
            val = float(val)
            pos1 = int(pos1)+1
            pos2 = int(pos2)+1

            extend1 = window
            extend2 = window
            if(i>0 and i< len(lines)-1):
                pre_line = lines[i-1].rstrip('\n')
                nex_line = lines[i+1].rstrip('\n')
                pre_chr,_,pre_pos,_ = pre_line.split('\t')
                nex_chr,nex_pos,_,_ = nex_line.split('\t')
                pre_pos = int(pre_pos)+1
                nex_pos = int(nex_pos)+1
                if(pre_chr == chromosome and pos1-pre_pos<2*window):
                    if(pos1-pre_pos<window):
                        extend1 = 0
                    else:
                        extend1 = pos1-pre_pos-window
                if(nex_chr == chromosome and nex_pos-pos2<window):
                    extend2 = nex_pos-pos2

            if ('chr' not in chromosome):
                chromosome = 'chr'+chromosome

            for pos in range(pos1-extend1,pos2+extend2):
                base = reference[pos-1]
                if(base == 'N'):
                    continue
                rpm = 0
                if(pos>=pos1 and pos<pos2):
                    rpm = val/depth
                block.append((pos,rpm,base))
                #block.append((pos,rpm))
    del reference
    gc.collect()
    print("Start generating block %d %d\n"%(start,end))
    return block

def split_chr_bedGraph(root_dir,input_file,strand,window,block_length,fa_file,depth,name):
    window /= 1.5
    window = int(window)
    chr_dict = dict()
    blocks = dict()
    block = []
    reference = dict()
    pre_pos = 0
    nex_pos = 0
    start_pos = 0
    block_num = 0
    count = 0
    chromosome = ''
    with open(input_file, "r") as bg_file:
        lines = bg_file.readlines()
        for i,line in enumerate(lines):
            line = line.rstrip('\n')
            chromosome,start,end,val = line.split('\t')
            if(len(chromosome)>5 or 'Y' in chromosome or 'M' in chromosome):
                continue
            val = float(val)
            pos1 = int(start)+1
            pos2 = int(end)+1

            extend1 = window
            extend2 = window
            if(i>0 and i< len(lines)-1):
                pre_line = lines[i-1].rstrip('\n')
                nex_line = lines[i+1].rstrip('\n')
                pre_chr,_,pre_pos,_ = pre_line.split('\t')
                nex_chr,nex_pos,_,_ = nex_line.split('\t')
                pre_pos = int(pre_pos)+1
                nex_pos = int(nex_pos)+1
                if(pre_chr == chromosome and pos1-pre_pos<2*window):
                    if(pos1-pre_pos<window):
                        extend1 = 0
                    else:
                        extend1 = pos1-pre_pos-window
                if(nex_chr == chromosome and nex_pos-pos2<window):
                    extend2 = nex_pos-pos2

            if ('chr' not in chromosome):
                chromosome = 'chr'+chromosome

            if chromosome not in chr_dict.keys():
                chr_dict[chromosome] = ''
                reference = get_genome_sequence('%s.%s.fa'%(fa_file,chromosome))
                block = []
                count = 0
                for pos in range(pos1-extend1,pos2+extend2):
                    base = reference[pos-1]
                    if(base == 'N'):
                        continue
                    rpm = 0
                    if(pos>=pos1 and pos<pos2):
                        rpm = val/depth
                    block.append((pos,rpm,base))
                    count += 1
            else:
                if(count > block_length and pos1-pre_pos >1000):
                    blocks['%s.%s_%s_%s'%(name,chromosome,strand,block_num)] = block
                    block_num += 1
                    block = []
                    start_pos = pos1
                    count = 0

                for pos in range(pos1-extend1,pos2+extend2):
                    base = reference[pos-1]
                    if(base == 'N'):
                        continue
                    rpm = 0
                    if(pos>=pos1 and pos<pos2):
                        rpm = val/depth
                    block.append((pos,rpm,base))
                    count += 1
            pre_pos = pos2
    blocks['%s.%s_%s_%s'%(name,chromosome,strand,block_num)] = block
    return blocks

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', default=None, help='out dir')
    parser.add_argument('--input_file', default=None, help='unstranded wig file')
    parser.add_argument('--input_plus', default=None, help='plus strand wig file')
    parser.add_argument('--input_minus', default=None, help='minus strand wig file')
    parser.add_argument('--fa_file',default=None,help='path to one line fa file')
    parser.add_argument('--keep_temp',default=None,help='if you want to keep temporary file, set to "yes"')
    parser.add_argument('--window', default=201, type=int, help='input length')
    parser.add_argument('--name', default='sample',help='sample name')
    parser.add_argument('--depth', default=1, type=float,help='total number of mapped reads( in millions)')
  
    argv = parser.parse_args()

    out_dir = argv.out_dir
    input_file = argv.input_file
    input_plus = argv.input_plus
    input_minus = argv.input_minus
    fa_file = argv.fa_file
    keep_temp =  argv.keep_temp
    window = argv.window
    name= argv.name
    depth = argv.depth

    return out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,depth


def Generate_windows(root_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,depth):
    block_length = 1e6
    
    if(root_dir[-1] == '/'):
        root_dir = root_dir[0:-1]
    if not os.path.exists(root_dir):
        os.makedirs(root_dir) 

    if(input_file is not None):
        strand = '+'
        if('bedgraph' in input_file.lower()):
            plus_blocks = split_chr_bedGraph(root_dir,input_plus,strand,window,block_length,fa_file,depth,name)
        else:
            sys.exit("input file extension should be wig or bedGraph")

    if(input_plus is not None):
        strand = '+'
        plus_blocks = dict()
        minus_blocks = dict()
        if('bedgraph' in input_plus.lower()):
            plus_blocks = split_chr_bedGraph(root_dir,input_plus,strand,window,block_length,fa_file,depth,name)
            for key,val in plus_blocks.items():
                key2 = key.replace('+','-')
                minus_blocks[key2] = val
        else:
            sys.exit("input file extension should be wig or bedGraph")
    if(input_minus is not None):
        strand = '-'
        mimus_blocks = dict()
        if('bedgraph' in input_minus.lower()):
            minus_block = split_chr_bedGraph(root_dir,input_minus,strand,window,block_length,fa_file,depth,name)
        else:
            sys.exit("input file extension should be wig or bedGraph")
    blocks = {**plus_blocks,**minus_block}
    print("Finished merging coverage and sequence information")
    return blocks
    #os.system('rm *%s/*.wig'%root_dir)



if __name__ == "__main__":
    #print(args())
    Generate_windows(*args())

