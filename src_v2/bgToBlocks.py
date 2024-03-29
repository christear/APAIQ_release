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

from pybedtools import BedTool
#
from Bio import SeqIO
#

def get_block_position(lines,win_width,block_length,keep_temp):
    _ext_width = int(win_width/1.5)
    blocks = []
    _ini_chr = ''
    _start_pos = 0
    _end_pos = 0
    _block_num = 0
    _count = 0
    _start_line = 0
    _skipped = 0
    #with open(input_bg, "r") as bg_file:
        #lines = bg_file.readlines()
    for _i,_line in enumerate(lines):
        _line = _line.rstrip('\n')
        _chr,_start,_end,_ = _line.split('\t')
        # ignore the Y chromosome and other fragments  
        if 'Y' in _chr or 'M' in _chr or '_' in _chr:
            _chr = _ini_chr
            _skipped += 1
            continue
        _pos1 = int(_start)+1
        _pos2 = int(_end)+1
        if _ini_chr == '':
            _start_pos = _pos1 - _ext_width
            _count = _ext_width + _pos2 - _pos1
            _ini_chr = _chr
            _start_line = _i
        elif _chr != _ini_chr:
            blocks.append([_block_num,_start_line,_i - 1,_ini_chr,_start_pos,_end_pos + _ext_width])
            _start_pos = _pos1 - _ext_width
            _count = _ext_width + _pos2 - _pos1
            _ini_chr = _chr
            _start_line = _i
            _block_num += 1
        else:
            _count += _pos2 - _pos1
            if _count > block_length:
                ### get error when processing the last line ..lines[_i + 1] with error: list index out of range
                if _i < len(lines) - 1:
                    _next_line = lines[_i + 1].rstrip('\n')
                    _next_chr,_next_s,_,_ = _next_line.split('\t')
                    if int(_next_s) + 1 - _pos2 > 2 * _ext_width:
                        #print("change block: {}\t{}\t{}".format(_chr,_pos2,_next_s))
                        _end_pos = _pos2 + _ext_width
                        blocks.append([_block_num,_start_line,_i,_chr,_start_pos,_end_pos])
                        _count = 0
                        _block_num += 1
                        _start_line = _i + 1
                        _start_pos = int(_next_s) + 1 - _ext_width
                    else:
                        continue
                else:
                    _end_pos = _pos2 + _ext_width
                    #_end_line = 
            else:
                _end_pos = _pos2 + _ext_width
        # add the last block 
    blocks.append([_block_num,_start_line,_i - _skipped,_chr,_start_pos,_end_pos])
    if keep_temp == 'yes':
        _tempf = 'temp.blocks.' + str(_chr) + ':' + str(_block_num)
        with open(_tempf,'w') as w:
            for _b in blocks:
                w.write('\t'.join(str(e) for e in _b) + '\n')
    #for _b in blocks:
    #    print(_b)                
    return blocks
    
#def Bedgraph_to_blocks(input_bg,fa_file,win_width,depth,block_pos,chromosome):
def Bedgraph_to_blocks(lines,fa_file,win_width,depth,block_pos,chromosome):    
    _ext_width = int(win_width/1.5)
    blocks = []
    #_block_num,_start_line,_end_line,_chr,_block_start,_block_end = block_pos
    #_start_line,_end_line,_block_start,_block_end,_block_num = block_pos
    _block_num,_start_line,_end_line,_,_block_start,_block_end = block_pos
    _chr = chromosome
    #print('Start generating block: {}:{}-{}'.format(_chr,_block_start,_block_end))
    # add one more nucleotide to the end of _ref to aviod errors 
    #_ref = BedTool.seq((_chr,_block_start,_block_end),fa_file)
    ### checking the chromosome size
    max_len = _block_end + 1
    for rec in SeqIO.parse(fa_file,'fasta'):
        if rec.id == _chr and len(rec.seq) < max_len:
            max_len = len(rec.seq)
    #_ref = BedTool.seq((_chr,_block_start,_block_end+1),fa_file)
    _ref = BedTool.seq((_chr,_block_start,max_len),fa_file)
    #with open(input_bg,"r") as bg:
        #lines = bg.readlines()
    _n = 0
    for _i,_f in enumerate(lines):
        if _i < _start_line:
            continue
        if _i > _end_line:
            break
        _n += 1
        #print(_i,_f)
        #_f = _f.rstrip('\n')
        _l_chr,_start,_end,_cov = _f.split('\t')
        if _l_chr != _chr:
            continue
        _start = int(_start)
        _end = int(_end)
        _cov = float(_cov)
        #print(_chr,_start,_end,_cov,_i)
        if _i == _start_line:
            _pos_start = _block_start
            ### some block only have one line 
            if _i != _end_line:
                next_f = lines[_i + 1].rstrip('\n')
                _,_next_start,_,_ = next_f.split('\t')
                _next_start = int(_next_start)
                if _next_start - _end < _ext_width:
                    _pos_end = _next_start
                else:
                    _pos_end = _end + _ext_width
            else:
                _pos_end = _block_end - 1
        elif _i == _end_line:
            _pos_start = _start - _ext_width
            _pos_end = _block_end - 1
        else:
            _pos_start = _start - _ext_width
            _pos_end = _end + _ext_width
            last_f = lines[_i - 1].rstrip('\n')
            next_f = lines[_i + 1].rstrip('\n')
            _,_,_last_end,_ = last_f.split('\t')
            _,_next_start,_,_ = next_f.split('\t')
            _last_end = int(_last_end)
            _next_start = int(_next_start)
            if _start - _last_end < 2 * _ext_width:
                if _start - _last_end < _ext_width:
                    _pos_start = _start
                else:
                    _pos_start = _last_end + _ext_width
                if _next_start - _end < _ext_width:
                    _pos_end = _next_start      
    ######
        for _pos in range(_pos_start,_pos_end):
            _res_pos = _pos - _block_start
            if _pos >= _block_end:
                print(_block_num,_chr,_block_start,_block_end,_pos)
                if _pos > _block_end:
                    break
            ### debuging 
            if _res_pos > len(_ref) - 1:
                print('### Error: index out of sequence boundary')
                print(len(_ref),_res_pos,_pos,_pos_start,_pos_end)
            _base = _ref[_res_pos]
            if _base == 'N':
                continue
            _rpm = 0.0
            if _pos >= _start and _pos < _end:
                _rpm = _cov/depth
            blocks.append((_pos,_rpm,_base))
                #print(_pos,_rpm,_base)
            #if _n > 10:
            #    break
    ####
    #lines = []
    _ref = ''            
    #print(blocks)
    return blocks

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_bg', default=None, help='input bedgraph file')
    parser.add_argument('--fa_file',default=None,help='path to genome fasta file')
    parser.add_argument('--window', default=201, type=int, help='length to window for PAS prediction')
    parser.add_argument('--depth', default=1, type=float,help='total number of mapped reads( in millions)')
    parser.add_argument('--block_length', default=1e5, type=int,help='block length')
    parser.add_argument('--keep_temp',default=None,help='if you want to keep temporary file, set to "yes"')
    argv = parser.parse_args()
    input_bg = argv.input_bg
    fa_file = argv.fa_file
    window = argv.window
    depth = argv.depth
    block_length = argv.block_length
    keep_temp =  argv.keep_temp
    return input_bg,fa_file,window,depth,block_length,keep_temp
    
if __name__ == "__main__":
    input_bg,fa_file,window,depth,block_length,keep_temp = args()
    bg_lines = []
    with open (input_bg,'r') as r:
        for _line in r:
            _line = _line.rstrip('\n')
            bg_lines.append(_line)
    blocks_pos = get_block_position(bg_lines,window,block_length,keep_temp)
    print(blocks_pos)
    blocks = Bedgraph_to_blocks(bg_lines,fa_file,window,depth,blocks_pos[-1],'chrX')
    print(sys.getsizeof(blocks))

    


