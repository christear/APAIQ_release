from PolyAModel import *
import re
import os, sys, copy, getopt, re, argparse
import random
import pandas as pd 
import numpy as np
from Bio.Seq import Seq
from TrimmedMean import TrimmedMean
import gc

def check(line1,line2,window):
    pos1 = line1[0]
    pos2 = line2[0]
    if(pos2-pos1==window-1):
        return True
    else:
        return False
# 
def collpase(strand,array,rst):
    sequence = ''
    coverage = []
    contain_N = False

    for line in array:
        _,rpm,base = line
        base = base.upper()
        if(base=='N'):
            contain_N = True
            break
        sequence += base
        coverage.append(rpm)
    
    if(not contain_N):
        if(strand == "-"):
            sequence = Seq(sequence)
            sequence = sequence.reverse_complement()

            coverage = coverage[::-1]
        trimMean = TrimmedMean([float(coverage[i]) for i in range(int(len(coverage)/2))])
        if(trimMean>=rst):
            return sequence,coverage
        else:
            return 0,0
    else:
        print("Discard windows containing N")
        return 0,0

# block is the output from bedgraph_to_block
def dataProcessing(chromosome,strand,block,window,rst):
    extend  = int(window/2)
    alphabet = np.array(['A', 'T', 'C', 'G'])
    data1 = []
    data2 = []
    PASID = []
    for i,line in enumerate(block):
        pos,_,base = line
        
        if(base.upper()=='N'):
            continue
        start = i-extend
        end   = i+extend
        if(start>0 and end+1<len(block)):
            if(not check(block[start],block[end],window)):
                continue
            sequence,coverage = collpase(strand,block[start:end+1],rst)
            if(sequence!=0):
                pas_id = '%s:%s:%s'%(chromosome,pos,strand)
                sequence = list(sequence)
                seq = np.array(sequence, dtype = '|U1').reshape(-1,1)
                seq_data = (seq == alphabet).astype(np.float32)
                coverage = np.array(coverage).astype(np.float32)
                data1.append(seq_data)
                data2.append(coverage)
                PASID.append(pas_id)
    if PASID != []:
        data1 = np.stack(data1).reshape([-1, window, 4])
        data2 = np.stack(data2).reshape([-1, window, 1])
        PASID = np.array(PASID)
        return data1 , data2,  PASID
    else:
        return
        
def args():
    parser = argparse.ArgumentParser(description="Evaluate each locus with RNAseq coverage exceed threshold and return prediction score")
    parser.add_argument("--model", help="the model weights file", required=True)
    parser.add_argument('--chromosome',default = 'chr1', help='chromosome')
    parser.add_argument('--strand',default = '+', help='strand')
    parser.add_argument("--RNASeqRCThreshold",default=0.05,type=float,help="RNA-Seq Coverage Threshold")
    parser.add_argument('--window', default=201, type=int, help='input length')
    parser.add_argument('--block', help = 'blocks derived from bedgraph_to_block')
    args = parser.parse_args()
    
    chromosome = args.chromosome
    strand = args.strand
    block = args.block
    model = args.model
    rst = args.RNASeqRCThreshold
    window=args.window

def Evaluate(chromosome,strand,block,model,rst,window):
    print('### Evaluation data processing')
    processed_data = dataProcessing(chromosome,strand,block,window,rst)
    #print('### Size of the processed data: {}'.format(sys.getsizeof(processed_data)))
    if processed_data != None:
        seq_data,cov_data,pas_id = processed_data
        #print("Finish processing data")
        block = []
        processed_data = []
        print("### Start Evaluating")
        keras_Model = PolyA_CNN(window)
        keras_Model.load_weights(model)
        #print('### Size of the model: {}'.format(sys.getsizeof(keras_Model)))
        pred = keras_Model.predict({"seq_input": seq_data, "cov_input": cov_data})    
        pred_out = []
        for i in range(len(pas_id)):
            pred_out.append((pas_id[i],pred[i][0]))
        return pred_out
    else:
        return
            
if __name__ == "__main__":
    pred_out = Evaluate(*args())

