#!/usr/bin/python

import sys,os
import argparse
# 
from bgToBlocks import get_block_position,Bedgraph_to_blocks
from evaluateBlock import Evaluate
from scanPredctions import Scan
from postScan import Postprocess,annotatePAS
#
#from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
import datetime
import gc
#
# change multiprocessing to concurrent.futures to check the memory usage problem...
# 
def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', default='out_dir', help='out dir')
    parser.add_argument('--input_file', default=None, help='unstranded wig file')
    parser.add_argument('--input_plus', default=None, help='plus strand wig file')
    parser.add_argument('--input_minus', default=None, help='minus strand wig file')
    parser.add_argument('--fa_file',default=None,help='genome fasta file')
    parser.add_argument('--keep_temp',default=None,help='if you want to keep temporary file, set to "yes"')
    parser.add_argument('--window', default=201, type=int, help='input length')
    parser.add_argument('--name', default='sample',help='sample name')
    parser.add_argument("--model", help="the model weights file", required=True)
    parser.add_argument("--RNASeqRCThreshold",default=0.05,type=float,help="RNA-Seq Coverage Threshold")
    parser.add_argument('--threshold', default=0,type=int,help='peak length lower than threshold will be fiter out')
    parser.add_argument('--penality', default=1,type=int,help='penality for prediction score lower than 0.5')
    parser.add_argument('--DB_file', default=None, help='polyA database file')
    parser.add_argument('--depth', default=1, type=float,help='total number of mapped reads( in millions)')
    parser.add_argument('--t', default = 30, type = int, help='number of process/thread/cpu')
    parser.add_argument('--block_length',default = 1e5,type = int, help='scanned length of each block')
    
    argv = parser.parse_args()
    out_dir = argv.out_dir
    input_file = argv.input_file
    input_plus = argv.input_plus
    input_minus = argv.input_minus
    fa_file = argv.fa_file
    keep_temp =  argv.keep_temp
    window   = argv.window
    name     = argv.name
    model    = argv.model
    rst      = argv.RNASeqRCThreshold
    threshold = argv.threshold
    penality  = argv.penality
    DB_file = argv.DB_file
    depth   = argv.depth
    thread = argv.t
    block_length = argv.block_length
    return out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread,block_length
    
def run_single_block(input_list):
    #global window
    #global keep_temp
    #global threshold
    #global penality
    #global DB_file
    #global model
    #global input_file
    #global fa_file 
    #global chromosome
    #global strand
    #print(input_list)
    #for e in [window,keep_temp,threshold,penality,DB_file,model,input_file,fa_file,chromosome,strand]:
        #print(e)
    #_start_line,_end_line,_start_pos,_end_pos,_block_num = input_list
    _block_num,_start_line,_end_line,_,_start_pos,_end_pos = input_list
    baseName = chromosome + '_' + strand + '_' + str(_block_num)
    print('### Generating block {}:{}-{}'.format(baseName,_start_pos,_end_pos))
    ####Generate blocks
    gw_start_time = datetime.datetime.now()
    #block = Bedgraph_to_blocks(input_file,fa_file,window,depth,input_list,chromosome)
    block = Bedgraph_to_blocks(lines,fa_file,window,depth,input_list,chromosome)
    gw_end_time = datetime.datetime.now()
    print('### Generating block {} used time: {}'.format(baseName,gw_end_time - gw_start_time))
    print('### Evaluating block {}:{}-{}'.format(baseName,_start_pos,_end_pos))
    ev_start_time = datetime.datetime.now()
    pred_out = Evaluate(chromosome,strand,block,model,rst,window)
    del block #destroyed the block reference
    gc.collect() #manually run garbage collection process
    ev_end_time = datetime.datetime.now()
    print('### Evaluating block {} used time: {}'.format(baseName,ev_end_time - ev_start_time))
    if pred_out is not None:
        print('### Scanning predictions in {}:{}-{}'.format(baseName,_start_pos,_end_pos))
        sc_start_time = datetime.datetime.now()
        forward_scan = Scan(pred_out,threshold,penality,method = 'forward')
        backward_scan = Scan(pred_out,threshold,penality,method = 'reverse')
        pas_dict = Postprocess(forward_scan,backward_scan,threshold)
        sc_start_time = datetime.datetime.now()
        print('### Scanning predictions {} used time: {}\n'.format(baseName,sc_start_time - sc_start_time))
        return pas_dict
    else:
        return
                
if __name__ == '__main__':
    out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread,block_length = args()
    start_time = datetime.datetime.now()
    # prepare output directory 
    if not os.path.exists(out_dir):
        print('### Creating output directory:{}'.format(out_dir))
        os.makedirs(out_dir)
    if input_file is not None:
        input_plus = input_file
        input_minus = input_file
    files = (input_plus,input_minus)
    strands = ('+','-')
    log = open('%s/%s.log'%(out_dir,name),'w')
    pas_out_list = []
    for i in range(2):
        input_file = files[i]
        strand = strands[i]
        print('### loading bedgraph file {}'.format(input_file))
        all_lines = {}
        with open(input_file,'r') as bg:
            #lines = bg.readlines()
            #for _f in lines:
            for _f in bg.readlines():
                _f = _f.rstrip('\n')
                _chr,_start,_end,_cov = _f.split('\t')
                if 'chrY' in _chr or '_' in _chr:
                    continue
                try:
                    all_lines[_chr].append(_f)
                except:
                    all_lines[_chr] = []
                    all_lines[_chr].append(_f)
            #lines = bg.readlines()
        #blocks_pos = get_block_position(lines,window,block_length)
        #blocks_input_list = {}
        #for bp in blocks_pos:
        #    _block_num,_start_line,_end_line,_chr,_start_pos,_end_pos = bp
        #    if 'chrY' in _chr:
        #        continue
        #    try:
        #        blocks_input_list[_chr].append([_start_line,_end_line,_start_pos,_end_pos,_block_num])
        #    except:
        #        blocks_input_list[_chr] = []
        #        blocks_input_list[_chr].append([_start_line,_end_line,_start_pos,_end_pos,_block_num])
        #    log.write('Blocks inf:{}\n'.format('\t'.join(str(e) for e in bp)))
        #for chromosome in blocks_input_list:
        for chromosome in all_lines:
            _chr_start_time = datetime.datetime.now()
            #if 'chr2' not in chromosome:
            #    continue
            print('### Processing {} in {} strand'.format(chromosome,strand))
            lines = all_lines[chromosome]
            #print(len(lines))
            #break
            chr_blocks = get_block_position(lines,window,block_length)
            log.write('Blocks inf:{}_{}\t{}\n'.format(chromosome,strand,'\t'.join(str(e) for e in chr_blocks)))
            #chr_blocks = blocks_input_list[chromosome]
            print('### {} blocks in {}_{} strand'.format(len(chr_blocks),chromosome,strand))
            with ProcessPoolExecutor(max_workers=thread) as executor:
                chr_out_pas = executor.map(run_single_block,chr_blocks)
                chr_anno_pas = annotatePAS(DB_file,chr_out_pas,chromosome,strand)
                pas_out_list += chr_anno_pas
            ### test run in single thread/cpu/core/process 
            #for cb in chr_blocks:
            #    each_po = run_single_block(cb)
            #    _i = 0
            #    for _pos in each_po:
            #        print(_pos,each_po[_pos])
            #        if _i > 10:
            #            break
            #        _i += 1
            #    break
            _chr_end_time = datetime.datetime.now()
            print('### Process {}_{} used time: {}'.format(chromosome,strand,_chr_end_time - _chr_start_time))
            #break
    log.close()
    
    out_file = '%s/%s.predicted.txt' %(out_dir,name)
    ww = open(out_file,'w')
    if(DB_file is not None): 
        ww.write('#chromosme\tstart\tend\tscore\tid\tstrand\tanno_id\tanno_source\tdistance\n')
    else:
        ww.write('#chromosme\tstart\tend\tscore\tid\tstrand\n')
    for p in pas_out_list:
        _line = '\t'.join(str(e) for e in p) + '\n'
        ww.write(_line)
    ww.close()
    end_time = datetime.datetime.now()    
    print("### Job Done and the time used: {}".format(end_time - start_time))
    
#if __name__ == '__main__':
#    main(*args())
