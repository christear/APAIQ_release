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

#block_length = 1000000 
# input file should be bedgraph 
def get_block_position(input_bg,win_width,block_length):
    _ext_width = int(win_width/1.5)
    blocks = []
    _ini_chr = ''
    _start_pos = 0
    _end_pos = 0
    _block_num = 0
    _count = 0
    _start_line = 0
    _skipped = 0
    with open(input_bg, "r") as bg_file:
        lines = bg_file.readlines()
        for _i,_line in enumerate(lines):
            _line = _line.rstrip('\n')
            _chr,_start,_end,_ = _line.split('\t')
            # ignore the Y chromosome and other patches 
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
        # add the last block 
        #blocks.append([_block_num,_start_line,_i,_chr,_start_pos,_end_pos])
        blocks.append([_block_num,_start_line,_i - _skipped,_chr,_start_pos,_end_pos])
    #for _b in blocks:
    #    print(_b)                
    return blocks
    
#def bedgraph_to_block(input_bg,strand,win_width,fa_file,depth,chr,start,end):
def Bedgraph_to_blocks(input_bg,fa_file,win_width,depth,block_pos):
    _ext_width = int(win_width/1.5)
    blocks = []
    _block_num,_start_line,_end_line,_chr,_block_start,_block_end = block_pos
    #print(block_pos)
    print('Start generating block: {}:{}-{}'.format(_chr,_block_start,_block_end))
    _ref = BedTool.seq((_chr,_block_start,_block_end),fa_file)
    #print(_ref[1:50])
    with open(input_bg,"r") as bg:
        lines = bg.readlines()
        _n = 0
        for _i,_f in enumerate(lines):
            #if _i < _start_line or _i > _end_line:
            if _i < _start_line:
                continue
            if _i > _end_line:
                break
            _n += 1
            _f = _f.rstrip('\n')
            _l_chr,_start,_end,_cov = _f.split('\t')
            if _l_chr != _chr:
                continue
            _start = int(_start)
            _end = int(_end)
            _cov = float(_cov)
            #print(_chr,_start,_end,_cov,_i)
            if _i == _start_line:
                _pos_start = _block_start
                next_f = lines[_i + 1].rstrip('\n')
                _,_next_start,_,_ = next_f.split('\t')
                _next_start = int(_next_start)
                if _next_start - _end < _ext_width:
                    _pos_end = _next_start
                else:
                    _pos_end = _end + _ext_width
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
            #print(_pos_start,_pos_end)
            for _pos in range(_pos_start,_pos_end):
                _res_pos = _pos - _block_start
                if _pos >= _block_end:
                    print(_block_num,_chr,_block_start,_block_end,_pos)
                    if _pos > _block_end:
                        break
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
    argv = parser.parse_args()
    input_bg = argv.input_bg
    fa_file = argv.fa_file
    window = argv.window
    depth = argv.depth
    block_length = argv.block_length
    return input_bg,fa_file,window,depth,block_length
    
if __name__ == "__main__":
    input_bg,fa_file,window,depth,block_length = args()
    blocks_pos = get_block_position(input_bg,window,block_length)
    blocks = Bedgraph_to_blocks(input_bg,fa_file,window,depth,blocks_pos[-1])
    print(sys.getsizeof(blocks))
    for b in blocks:
        print(b)
        #print(len(b))
    


