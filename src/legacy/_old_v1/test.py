#!/usr/bin/env python
# coding: utf-8

import sys
#from pybedtools import BedTool
#from bedgraph_to_blocks import get_block_position,Bedgraph_to_blocks
#from evaluate import Evaluate
#from postprocess import get_predict_score,merge_predict_pos
from concurrent.futures import ProcessPoolExecutor

with open('RNAseq.depth.bedGraph','r') as b:
    bg = b.readlines()
    print(sys.getsizeof(bg))

def make_dict(a):
    o_dict = dict()
    o_dict[str(a)] = a
    if a != 4:
        return o_dict
    else:
        return
        
def make_tuple(a):
    if a != 4:
        return (str(a),a)
    else:
        return None
        
test_a = [1,2,3,4,5]
#for _a in test_a:
#    _ta = do_something(_a)
#    print(_ta)
with ProcessPoolExecutor(max_workers=4) as executor:
    test_out1 = executor.map(make_dict,test_a)
    test_out2 = executor.map(make_tuple,test_a)
    
print(type(test_out1))
print(type(test_out2))
for _t in test_out1:
    print(_t)

for _t in test_out2:
    print(_t)    

a = [1,2]
def test(a):
    print(a)

t = test("yes")
a.append(t)


model='model/snu398_model.ckpt'
rst=0.05
window=201
keep_temp='yes'
block_length=1e5
depth=1
out_dir='out_dir/test100000'
input_bg='chr8.depth.bedGraph'
baseName='test100000.chr12_+_1060'
fa_file='/home/zhanb0d/c2066/TCGA/annotation/hg38.intg.virg.fa'
threshold=0
penality=1

#blocks_pos = get_block_position(input_bg,window,block_length)
#blocks = Bedgraph_to_blocks(input_bg,fa_file,window,depth,blocks_pos[-1])

# 
#for b in blocks:
#    print(b)
    
#
"""
test_file='out_dir/chr12-13100000/chr12-13100000.predicted.txt'
long_bed_str = ''
with open(test_file,'r') as t:
    for _line in t.readlines():
        if 'pasid' in _line:
            continue
        #print(_line)
        _pas_id,_,_,_score = _line.rstrip().split('\t')
        #print(_pas_id)
        _chr,_loci,_strand = _pas_id.split(':')
        #print(_chr,_loci,_strand)
        _loci = int(_loci)
        _bl = '\t'.join(str(e) for e in [_chr,_loci-1,_loci,_pas_id,_score,_strand]) + '\n'
        #print(_bl)
        long_bed_str += _bl
        #print(_bl)
        #break
print(sys.getsizeof(long_bed_str))
pred_bed = BedTool(long_bed_str,from_string=True)
pred_bed = pred_bed.sort()
_i = 0
for _b in pred_bed:
    print(_b)
    _i += 1
    if _i > 10:
        break
        
#anno_bed = BedTool('polyADB3_gencode.pAs.bed')
anno_pred_bed = pred_bed.closest('polyADB3_gencode.pAs.bed',s = True,D = 'a')
for _ab in anno_pred_bed:
    print(_ab[0],_ab[1],_ab[2],_ab[4],_ab[3],_ab[5],_ab[9],_ab[10],_ab[11])
    _i += 1
    if _i > 20:
        break

#Evaluate(baseName,blocks,model,out_dir,rst,window,keep_temp)  
forward_file=out_dir+"/maxSum/%s.forward.%d.%d.txt"%(baseName,threshold,penality)
backward_file=out_dir+"/maxSum/%s.backward.%d.%d.txt"%(baseName,threshold,penality)
forward_pas_dict = get_predict_score(forward_file,threshold)
backward_pas_dict = get_predict_score(backward_file,threshold)
print(len(forward_pas_dict))
print(len(backward_pas_dict))
if len(forward_pas_dict) > 0 and len(backward_pas_dict) > 0:
    print('list is not empty')
else:
    print('list is empty')
"""
 