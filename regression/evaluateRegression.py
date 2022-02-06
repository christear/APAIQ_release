from RegressionModel import *
import os, sys, copy, re, argparse
import random
import pandas as pd 
import numpy as np
from collections import defaultdict


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

def read_bg(bg_file,depth):
    bg_dict = defaultdict(dict)
    with open(bg_file,'r') as bg:
        for line in bg.readlines():
            line = line.rstrip('\n')
            chromosome,start,end,rpm = line.split('\t')
            if ('chr' not in chromosome):
                chromosome = 'chr'+chromosome
            if(len(chromosome)>5 or 'Y' in chromosome or 'M' in chromosome):
                continue
            rpm = float(rpm)/depth
            start = int(start)+1
            end = int(end)+1
            for i in range(start,end):
                bg_dict[chromosome][i] = rpm
    return bg_dict



def dataProcessing(input_file,window,bg_dict_plus,bg_dict_minus,threshold):
    data1 = []
    PASID = []
    with open(input_file,'r') as inp:
        for line in inp.readlines():
            if (line.startswith('#')):
                continue
            line = line.rstrip('\n')
            chromosome,start,end,score,pasid,strand = line.split('\t')[0:6]
            pos = int(end)
            score = float(score)
            if(score >=threshold):
                coverage = np.zeros(window)
                extend  = int(window/2)
                for i in range(window):
                    if(strand == '+'):
                        new_pos = pos+i-extend
                        pos_dict_plus = bg_dict_plus[chromosome]
                        if(new_pos in pos_dict_plus):
                            coverage[i] = pos_dict_plus[new_pos]
                    elif(strand == '-'):
                        new_pos = pos-i+extend
                        pos_dict_minus = bg_dict_minus[chromosome]
                        if(new_pos in pos_dict_minus):
                            coverage[i] = pos_dict_minus[new_pos]
                coverage = np.array(coverage).astype(np.float32)
                data1.append(coverage)
                PASID.append(pasid)
    if PASID != []:
        data1 = np.stack(data1).reshape([-1, window, 1])
        PASID = np.array(PASID)
        return data1,PASID
    else:
        return
        
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
    window=args.window
    
    return model,factor_path,input_file,input_plus,input_minus,pas_file,out,threshold,depth,window

def Evaluate(model,factor_path,input_file,input_plus,input_minus,pas_file,out,threshold,depth,window):
    out_dir = os.path.dirname(out)
    if (not os.path.exists(out_dir)) and (out_dir != ''):
        os.makedirs(out_dir)

    print('### Evaluation data processing')
    factor_dict = read_factor(factor_path)
    if input_file is not None:
        bg_dict_plus = read_bg(input_plus,depth)
        bg_dict_minus = bg_dict_plus
    elif (input_plus is not None) and (input_minus is not None):
        bg_dict_plus = read_bg(input_plus,depth)
        bg_dict_minus = read_bg(input_minus,depth)
    cov_data,pasid  = dataProcessing(pas_file,window,bg_dict_plus,bg_dict_minus,threshold)
    cov_data = normalize(cov_data,factor_dict)
    print("### Start Evaluating")
    keras_Model = Regression_CNN(window)
    keras_Model.load_weights(model)
    pred = keras_Model.predict({"cov_input":cov_data})    
    OUT=open(out,'w')
    OUT.write("pas_id\tpredict_readCount\n")
    for i in range(len(pred)):
        predict = pred[i][0]
        predict_readCount = inverse_normalize(predict,factor_dict)
        OUT.write('%s\t%s\n'%(pasid[i],predict_readCount))
    OUT.close()
            
if __name__ == "__main__":
    Evaluate(*args())