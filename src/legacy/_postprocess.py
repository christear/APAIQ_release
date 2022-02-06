#!/usr/bin/env python -w

import argparse
import os
from collections import OrderedDict


def get_predict_score(scan_file,Threshold):
    pas_dict = OrderedDict()
    with open(scan_file,'r') as f:
        for line in f:
            if 'pas_id' in line:
                continue
            (pas_id,maxPoint,maxPos,start,end,peak) = line.split('\t')
            maxPoint = float(maxPoint)
            maxPos   = float(maxPos)
            peak     = float(peak)
            if(maxPoint < Threshold):
                continue
            pas_dict[maxPos] = [maxPoint,start,peak]
    return pas_dict

def merge_predict_pos(forward_pas_dict,backward_pas_dict):
    forward_pos = sorted(forward_pas_dict, key=lambda x:int(x))
    backward_pos = sorted(backward_pas_dict,key=lambda x:int(x))
    pas_dict = dict()
    i=0
    j=0
    while(i<len(forward_pos) and j < len(backward_pos)):
        pos1 = forward_pos[i]
        pos2 = backward_pos[j]
        maxPoint1 = forward_pas_dict[pos1][0]
        maxPoint2 = backward_pas_dict[pos2][0]
        start1 = forward_pas_dict[pos1][1]
        start2 = backward_pas_dict[pos2][1]
        if(pos2>pos1):
            i += 1
        elif(start1 > start2):
            j += 1
        else:
            pos = "%.f"%((pos1+pos2)/2)
            pas_dict[pos] = (maxPoint1+maxPoint2)/2
            i += 1
            j += 1

    return pas_dict        

def annotated(DB_file,pas_dict,chromosome,strand):
    start = min(pas_dict, key=lambda x:int(x))
    end = max(pas_dict, key=lambda x:int(x))
    start = int(start)-10000
    end   = int(end)  + 10000
    nearestID = dict()
    nearest    = dict()
    if('chr' not in chromosome):
        chromosome = 'chr'+chromosome
    with open(DB_file,'r') as f:
        for line in f:
            pas_id,pas_type,chro,pos,srd = line.split('\t')[0:5]
            if(not pos.isdigit()):
                continue
            pos = int(pos)
            if('chr' not in chro):
                chro = 'chr'+chro
            if(srd == strand and chro == chromosome and start <= pos and end >= pos):
                for key,val in pas_dict.items():
                    diff = int(key)-pos
                    if(key not in nearest.keys() or abs(diff)<abs(nearest[key])):
                        nearest[key] = diff
                        nearestID[key] = '%s:%d:%s'%(chro,pos,srd)
    return nearest,nearestID

def save_file(pas_dict,out,chro,srd,nearest=None,nearestID=None):
    if(nearest is not None):
        OUT=open(out,'w')
        for key,val in pas_dict.items():
            pasid = '%s:%s:%s'%(chro,key,srd)
            try:
                gtid = nearestID[key]
            except:
                gtid = 'NA' 
            try:
                gt_diff = nearest[key]
            except:
                gt_diff = 1e9
            #gtid  = nearestID[key]
            #gt_diff = nearest[key]
            OUT.write("%s\t%s\t%d\t%.1f\n"%(pasid,gtid,gt_diff,val))
        OUT.close()
    else:
        OUT=open(out,'w')
        for key,val in pas_dict.items():
            pasid = '%s:%s:%s'%(chro,key,srd)
            OUT.write("%s\t%.1f\n"%(pasid,val))
        OUT.close()
    return 0

def args():
    ### Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--DB_file', default=None, help='polyA database file')
    parser.add_argument('--out_dir', default=None, help='out dir')
    parser.add_argument('--baseName', help='baseName')
    parser.add_argument('--threshold', default=0,type=int,help='peak length lower than threshold will be fiter out')
    parser.add_argument('--penality', default=1,type=int,help='penality for prediction score lower than 0.5')
    args = parser.parse_args()

    DB_file = args.DB_file
    baseName = args.baseName
    threshold = args.threshold
    penality  = args.penality
    out_dir=args.out_dir

    return DB_file,baseName,threshold,penality,out_dir

def Postprocess(DB_file,baseName,threshold,penality,out_dir):
    print("Start postprocessing"+baseName)
    _,block = baseName.split('.')
    chromosome,strand,_ = block.split('_')

    if(out_dir[-1] == '/'):
        out_dir = out_dir[0:-1]
    forward_file=out_dir+"/maxSum/%s.forward.%d.%d.txt"%(baseName,threshold,penality)
    forward_pas_dict = get_predict_score(forward_file,threshold)
    backward_file=out_dir+"/maxSum/%s.backward.%d.%d.txt"%(baseName,threshold,penality)
    backward_pas_dict = get_predict_score(backward_file,threshold)
    pas_dict = merge_predict_pos(forward_pas_dict,backward_pas_dict)
    out=out_dir+"/maxSum/%s.bidirection.%d.%d.txt"%(baseName,threshold,penality)
    if(DB_file is None):
        save_file(pas_dict,out,chromosome,strand)
    else:
        nearest,nearestID = annotated(DB_file,pas_dict,chromosome,strand)
        save_file(pas_dict,out,chromosome,strand,nearest,nearestID)
    #print("Finish postprocessing"+baseName)
    return 0

if __name__ == "__main__":
    Postprocess(*args())
