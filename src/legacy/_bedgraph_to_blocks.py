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

from pybedtools import BedTool

block_length = 1000000 
# input file should be bedgraph 
def get_block_position(input_bg,win_width,block_length):
    _ext_width = int(win_width/1.5)
    blocks = []
    _ini_chr = ''
    _start_pos = 0
    _end_pos = 0
    _block_num = 0
    _count = 0
    _start_line = 1
    with open(input_bg, "r") as bg_file:
        lines = bg_file.readlines()
        for _i,_line in enumerate(lines):
            _line = _line.rstrip('\n')
            _chr,_start,_end,_ = _line.split('\t')
            _pos1 = int(_start)+1
            _pos2 = int(_end)+1
            if _ini_chr == '':
                _start_pos = _pos1 - _ext_width
                _count = _ext_width + _pos2 - _pos1
                _ini_chr = _chr
            elif _chr != _ini_chr:
                blocks.append([_block_num,_start_line,_i,_ini_chr,_start_pos,_end_pos + _ext_width])
                _start_pos = _pos1 - _ext_width
                _count = _ext_width + _pos2 - _pos1
                _ini_chr = _chr
                _start_line = _i + 1
                _block_num += 1
            else:
                _count += _pos2 - _pos1
                if _count > block_length:
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
        # ignore the Y chromosome 
        #blocks.append([_block_num,_start_line,_i,_chr,_start_pos,_end_pos])
    for _b in blocks:
        print(_b)                
    return blocks
    
#def bedgraph_to_block(input_bg,strand,win_width,fa_file,depth,chr,start,end):
def bedgraph_to_block(input_bg,fa_file,win_width,depth,block_pos):
    _ext_width = int(win_width/1.5)
    blocks = []
    _block_num,_start_line,_end_line,_chr,_block_start,_block_end = block_pos
    print('Start generating block: {}:{}-{}'.format(_chr,_block_start,_block_end))
    #bg = BedTool(input_bg)
    #region = BedTool(_chr + '\t' + str(_block_start) + '\t' + str(_block_end),from_string=True)
    #sub_bg = bg.intersect(region,u = True)
    sub_bg = bg[_start_line:_end_line]
    _ref = BedTool.seq((_chr,_block_start,_block_end),fa_file)
    print(_ref[1:50])
    _i = 0
    for f in sub_bg:
        #_chr = f.chrom
        _start = f.start
        _end = f.stop
        _cov = float(f.name)
        _rpm = _cov/depth
        print(f.chrom,_start,_end,_rpm)
        if _i == 0:
            _pos_start = _block_start
            next_f = bg[_start_line + _i + 1]
            if next_f.start - _end < _ext_width:
                _pos_end = next_f.start
            else:
                _pos_end = _end + _ext_width
        elif _i == _end_line - _start_line:
            _pos_start = _start - _ext_width
            _pos_end = _block_end
        else:
            last_f = bg[_start_line + _i - 1]
            next_f = bg[_start_line + _i + 1]  
            if _start - last_f.end < 2 * _ext_width:
                if _start - last_f.end < _ext_width:
                    _pos_start = _start
                else:
                    _pos_start = last_f.end + _ext_width
            if next_f.start - _end < _ext_width:
                _pos_end = next_f.start      
######
        for _pos in range(_pos_start,_pos_end):
            _res_pos = _pos - _block_start
            _base = _ref[_res_pos]
            if _base == 'N':
                continue
            _rpm = 0
            if _pos >= _start and _pos < _end:
                _rpm = _cov/depth
            blocks.append((_pos,_rpm,_base))
        _i += 1
        if _i > 10:
            break
    #reference = get_genome_sequence('%s.%s.fa'%(fa_file,chromosome))
    print(blocks)
    return blocks

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_bg', default=None, help='input bedgraph file')
    parser.add_argument('--fa_file',default=None,help='path to genome fasta file')
    parser.add_argument('--keep_temp',default=None,help='if you want to keep temporary file, set to "yes"')
    parser.add_argument('--window', default=201, type=int, help='length to window for PAS prediction')
    parser.add_argument('--depth', default=1, type=float,help='total number of mapped reads( in millions)')
    argv = parser.parse_args()
    input_bg = argv.input_bg
    fa_file = argv.fa_file
    keep_temp =  argv.keep_temp
    window = argv.window
    depth = argv.depth
    block_len = 1000000
    #return input_bg,fa_file,keep_temp,window,depth
    #return input_bg,window,block_len
    return input_bg,fa_file,window,depth,[8,1574072,1746262,'chr1',95199112,113701084]
    
if __name__ == "__main__":
    #print(args())
    #get_block_position(*args())
    bedgraph_to_block(*args())

