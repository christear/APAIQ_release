from RegressionModel import *
import os, sys, copy, re, argparse
import random
import pandas as pd 
import numpy as np
from collections import defaultdict
#
from pybedtools import BedTool

def normalize(data,factor_dict):
    data = np.log(data+factor_dict['pseudo_count'])
    data = (data-factor_dict['data_mean'])/factor_dict['data_std']
    return data

def inverse_normalize(rpm,factor_dict):
    rpm = np.exp(rpm*factor_dict['label_std']+factor_dict['label_mean'])
    return rpm

def read_factor(file_path):
    factor = dict()
    with open(file_path,'r') as f:
        for line in f.readlines():
            line = line.rstrip('\n')
            name,value = line.split('=')
            value  = float(value)
            factor[name] = value
    return factor
    
def read_bg(bg_file,depth,pas_bed,window,genome):
    bg_dict = defaultdict(dict)
    a_bed = BedTool(bg_file)
    b_bed = BedTool(pas_bed)
    b_bed = b_bed.slop(genome = genome,b = window + 100)
    a_with_b = a_bed.intersect(b_bed,u = True)
    for ab in a_with_b:    
        for i in range(ab.start + 1, ab.end + 1):
            bg_dict[ab.chrom][i] = float(ab.name)/depth
    return bg_dict
    
def dataProcessing(pas_pos_list,strand,window,bg_dict):
    data1 = []
    for _chr,_pos in pas_pos_list:
        coverage = np.zeros(window)
        extend = int(window/2)
        for i in range(window):
            if(strand == "+"):
                new_pos = _pos + i - extend
            elif(strand == "-"):
                new_pos = _pos -i + extend
            #
            if(new_pos in bg_dict[_chr]):
                coverage[i] = bg_dict[_chr][new_pos]
        coverage = np.array(coverage).astype(np.float32)
        data1.append(coverage)
    data1 = np.stack(data1).reshape([-1, window, 1]) 
    return data1            
        
def args():
    parser = argparse.ArgumentParser(description="Evaluate each locus with RNAseq coverage exceed threshold and return prediction score")
    parser.add_argument("--model", help="the model weights file", required=True)
    parser.add_argument("--factor_path", help="normalization file path", required=True)
    parser.add_argument('--input_file', default=None, help='unstranded bedGraph file')
    parser.add_argument('--input_plus', default=None, help='plus strand bedGraph file')
    parser.add_argument('--input_minus', default=None, help='minus strand bedGraph file')
    parser.add_argument('--pas_file', default=None, help='pAs location file to be predicted expression level')
    parser.add_argument("--out", default='predicted_rpm.txt',help="output file path")
    parser.add_argument('--threshold', default=12,type=int,help    ='peak length lower than threshold will be fiter out')
    parser.add_argument('--depth', default=1, type=float,help=    'total number of mapped reads( in millions)')
    parser.add_argument('--window', default=1001, type=int, help='input length')
    parser.add_argument('--genome', default='hg38', help='assembly name of the genome. i.e. hg19, hg38, mm10')
    args = parser.parse_args()
    
    model = args.model
    factor_path = args.factor_path
    input_file = args.input_file
    input_plus = args.input_plus
    input_minus = args.input_minus
    pas_file  = args.pas_file
    out = args.out
    threshold = args.threshold
    depth = args.depth
    window = args.window
    genome = args.genome
    return model,factor_path,input_file,input_plus,input_minus,pas_file,out,threshold,depth,window,genome

def Evaluate(model,factor_path,input_file,input_plus,input_minus,pas_file,out,threshold,depth,window,genome):
    out_dir = os.path.dirname(out)
    if (not os.path.exists(out_dir)) and (out_dir != ''):
        os.makedirs(out_dir)
    print('### Evaluation based on regression model')
    print('### Loading RNAseq coverage')
    if input_file is not None:
        bg_dict = read_bg(input_file,depth,pas_file,window,genome)
    else:
        bg_dict_plus = read_bg(input_plus,depth,pas_file,window,genome)
        bg_dict_minus = read_bg(input_minus,depth,pas_file,window,genome)
    print('### Reading facotrs from {}'.format(factor_path))
    factor_dict = read_factor(factor_path)
    print('### loading the PAS list from {}'.format(pas_file))
    sorted_pas_pos = {"+":[],"-":[]}
    sorted_pas_id = {"+":[],"-":[]}
    original_pas_id = []
    with open(pas_file,'r') as pb:
        for _line in pb:
            if (_line.startswith('#')):
                continue
            _line = _line.rstrip('\n')
            _chr,_,_end,_pas_id,_score,_strand = _line.split('\t')[0:6]
            _score = float(_score)
            _pos = int(_end)
            if(_score > threshold):
                original_pas_id.append(_pas_id)
                sorted_pas_id[_strand].append(_pas_id)
                sorted_pas_pos[_strand].append([_chr,_pos])
    print('### Running regression')
    rpms = {}
    for _strand in sorted_pas_pos:
        pas_pos_list = sorted_pas_pos[_strand]
        pas_id_list = sorted_pas_id[_strand]
        print('### Processing {} PAS on {} strand'.format(len(pas_id_list),_strand))
        if _strand == "+":
            if input_file is not None:
                cov_data = dataProcessing(pas_pos_list,_strand,window,bg_dict)
            else:
                cov_data = dataProcessing(pas_pos_list,_strand,window,bg_dict_plus)
        else:
            if input_file is not None:
                cov_data = dataProcessing(pas_pos_list,_strand,window,bg_dict)
            else:
                cov_data = dataProcessing(pas_pos_list,_strand,window,bg_dict_minus)
        cov_data = normalize(cov_data,factor_dict)
        keras_Model = Regression_CNN(window)
        keras_Model.load_weights(model)
        pred = keras_Model.predict({"cov_input":cov_data})
        for i in range(len(pred)):
            #print(pas_id_list[i],pred[i][0])
            predict_readCount = inverse_normalize(pred[i][0],factor_dict)
            rpms[pas_id_list[i]] = predict_readCount
    with open(out,'w') as w:
        w.write("pas_id\tpredict_readCount\n")
        for _pas_id in original_pas_id:
            #print(_pas_id,rpms[_pas_id])
            w.write('%s\t%s\n'%(_pas_id,rpms[_pas_id]))
    print('### results have been written to {}'.format(out))
    
            
if __name__ == "__main__":
    Evaluate(*args())