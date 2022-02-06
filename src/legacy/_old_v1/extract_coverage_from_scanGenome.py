#!/usr/bin/env python
# coding: utf-8


import argparse
from Bio.Seq import Seq
import re,sys
import pipes
import pprint
import tempfile


def read_pas(pas_file):
    f = open(pas_file, 'r')
    pas_dict = dict()
    f.readline() #sikp header
    for line in f.readlines():
        line = line.rstrip('\n')
        pas_id,pas_type,chromosome,pos,strand,symbol = line.split('\t')[0:6]
        pas_dict[pas_id] = [pas_type,symbol]
    f.close()
    return pas_dict

def read_gtf(pas_file,symbol_index,biotype_index):
    if re.findall('gz$',pas_file):
        p = pipes.Template()
        p.append("zcat", '--')
        f = p.open(pas_file, 'r')
    else:
        f = open(pas_file, 'r')
    pas_dict = dict()
    for line in f.readlines():
        if(re.findall('^#',line)):
            continue
        line = line.rstrip('\n')
        chromosome,_,feature,start,end,_,strand,_,attribute = line.split('\t')
        if(not re.findall('chr',chromosome)):
            chromosome = 'chr'+chromosome

        
        if(feature == 'transcript'):
            Attribute_array = attribute.split(' ')
            symbol = Attribute_array[symbol_index].strip(';').strip('\"')
            biotype = Attribute_array[biotype_index].strip(';').strip('\"')
        
            if(strand == "+"):
                pas_id=chromosome+':'+end+":"+strand
            elif(strand == "-"):
                pas_id=chromosome+':'+start+":"+strand
            else:
                raise Exception('strand='+strand+' not correct')
            
            pas_dict[pas_id] = [biotype,symbol]
    f.close()
    return pas_dict



def output(pas_dict,scan_file,out,window,species):
    
    extend  = int(window/2)
    
    f = open(scan_file,'r')
    lines = f.readlines()
    
    
    ww = open(out,'w')


    for i,line in enumerate(lines):
        line = line.rstrip('\n')
        pas_id,rpm,base = line.split('\t')
        if pas_id in pas_dict.keys():
            pas_type,symbol = pas_dict[pas_id]
            start = i-extend
            end   = i+extend
            if(start>0 and end+1<len(lines)):
                if(check(lines[start],lines[end],window)):
                    collpase(pas_id,pas_type,symbol,lines[start:end+1],ww,species)
    f.close()
    ww.close()
    
    
def get_motif(sequence,species):
    if(species == 'mm10'):
        motifs = ['AAGAAA', 'AATAAA', 'AATACA', 'AATAGA', 'AATATA','AATGAA',
                  'ACTAAA', 'AGTAAA', 'ATTAAA', 'CATAAA', 'GATAAA', 'TATAAA', 'AAAAAG']
    elif(species == 'hg38'):
        motifs = ['AAGAAA', 'AATAAA', 'AATACA', 'AATAGA', 'AATATA', 'ACTAAA',
                  'AGTAAA', 'ATTAAA', 'CATAAA', 'GATAAA', 'TATAAA', 'AAAAAG']
    else:
        raise Exception('species not correct')
        
    index1 = int(len(sequence)/2)-37
    index2 = index1+32
    subseq = sequence[index1:index2]
    motif_number = 0
    for motif in motifs:
        
        match = re.findall(motif,str(subseq))
        motif_number += len(match)
        
    
    return 'motif='+str(motif_number)
    
def check(line1,line2,window):
    line1 = line1.rstrip('\n')
    line2 = line2.rstrip('\n')
    _,pos1,_ = line1.split(':')
    _,pos2,_ = line2.split(':')
    pos1 = int(pos1)
    pos2 = int(pos2)
    if(pos2-pos1==window-1):
        return True
    else:
        return False
    
    
def collpase(pas_id,pas_type,symbol,array,ww,species):
    
    sequence = ''
    coverage = []
    contain_N = False
    for line in array:
        line = line.rstrip('\n')
        _,rpm,base = line.split('\t')
        base = base.upper()
        if(base=='N'):
            contain_N = True
            break
        sequence += base
        coverage.append(rpm)
    
    if(not contain_N):
        chromosome,pos,strand = pas_id.split(':')
        if(strand == "-"):
            sequence = Seq(sequence)
            sequence = sequence.reverse_complement()
            coverage.reverse()
        motif = get_motif(sequence,species)
        ww.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                 %(pas_id,pas_type,chromosome,pos,strand, symbol,motif,sequence,'\t'.join(coverage)))
        
    else:
        print("Discard PAS:%s containing N"%pas_id)
        
        
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--scan_file', default=None, help='scan transcriptom file')
    parser.add_argument('--pas_file', default=None, help='pas file')
    parser.add_argument('--window', default=201,type=int, help='window size')
    parser.add_argument('--species', default='hg38', help='pas file')
    parser.add_argument('--file_type', default='polyADB', choices={'polyADB', 'ensembl','gencode'},help='pas file')
    argv = parser.parse_args()

    scan_file = argv.scan_file
    pas_file = argv.pas_file
    window   = argv.window
    species = argv.species
    file_type = argv.file_type
    
    if(file_type == 'polyADB'):
        pas_dict = read_pas(pas_file)
        out = scan_file.replace('chr','dbcoverage_chr')
        
    elif(file_type == 'ensembl'):
        pas_dict = read_gtf(pas_file,9,19)
        out = scan_file.replace('chr','enscoverage_chr')
    elif(file_type == 'gencode'):
        pas_dict = read_gtf(pas_file,7,9)
        out = scan_file.replace('chr','gencoverage_chr')
    else:
        sys.exit("file type should be polyADB, ensembl and gencode")
    
    output(pas_dict,scan_file,out,window,species)
    
