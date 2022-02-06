#!/usr/bin/env python -w
#
#
import argparse
import os
from collections import OrderedDict
from pybedtools import BedTool


def get_predict_score(scan_out,threshold):
    pas_dict = OrderedDict()
    for pas_id,maxPoint,maxPos,start,end,peak in scan_out:
        if(maxPoint < threshold):
            continue
        pas_dict[maxPos] = [maxPoint,start,peak]
    return pas_dict

def merge_predict_pos(forward_pas_dict,backward_pas_dict):
    forward_pos = sorted(forward_pas_dict, key=lambda x:int(x))
    backward_pos = sorted(backward_pas_dict,key=lambda x:int(x))
    pas_dict = dict()
    #pas_out = []
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
            #pas_out.append((pos,(maxPoint1+maxPoint2)/2))
            i += 1
            j += 1
    return pas_dict
    #return pas_out        

# convert output in generator Class from multipleprocessing to dict 
def generator_to_dict(pas_generator):
    if type(pas_generator) is dict:
        return pas_generator
    else:
        pas_dict = {}
        for pg in pas_generator:
            if pg is not None:
                for _pos in pg:
                    pas_dict[_pos] = pg[_pos]                
        return pas_dict
        
# DB_file should be in sorted bed format 
def annotatePAS(DB_file,pas_generator,chromosome,strand):
    if DB_file is not None:
        long_bed_str = ''
        if pas_generator is not None:
            pas_dict = generator_to_dict(pas_generator)
            _i = 1
            for _pos in pas_dict:
                _bl = '\t'.join(str(e) for e in [chromosome,int(_pos) - 1,_pos,pas_dict[_pos],chromosome + ":" + strand + ":" + str(_i),strand]) + '\n'
                long_bed_str += _bl
                _i += 1
            pas_bed = BedTool(long_bed_str,from_string=True)
            pas_bed = pas_bed.sort()
            anno_pas_bed = pas_bed.closest(DB_file,s = True, D = 'b')
            annotated_pas_out = []
            for _apb in anno_pas_bed:
                annotated_pas_out.append((_apb[0],_apb[1],_apb[2],_apb[3],_apb[4],_apb[5],_apb[9],_apb[10],_apb[12]))
            return annotated_pas_out
        else:
            return
    else:
        if pas_generator is not None:
            pas_out = []
            pas_dict = generator_to_dict(pas_generator)
            _i = 1
            for _pos in pas_dict:
                pas_out.append((chromosome,int(_pos) - 1,_pos,sub_pas[_pos],chromosome + ":" + strand + ":" + str(_i),strand))
                _i += 1
            return pas_out
        else:
            return
            
            
def args():
    ### Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--forward_scan', default=None, help='output from forword scanning')
    parser.add_argument('--backward_scan', default=None, help='output from backward scanning')
    args = parser.parse_args()
    forward_scan = args.forward_scan
    backward_scan = args.backward_scan
    return forward_scan,backward_scan

def Postprocess(forward_scan, backward_scan,threshold):
    forward_score = get_predict_score(forward_scan,threshold)
    backward_score = get_predict_score(backward_scan,threshold)
    if len(forward_score) > 0 and len(backward_score) > 0:
        pas_dict = merge_predict_pos(forward_score,backward_score)
        return pas_dict
    else:
        return
    
if __name__ == "__main__":
    Postprocess(*args())
