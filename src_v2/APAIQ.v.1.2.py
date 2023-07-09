#!/usr/bin/python

import sys,os
import argparse
# 
from bgToBlocks import get_block_position,Bedgraph_to_blocks
from evaluateBlock import Evaluate
from evaluateBlock import EvaluateDeepPASS
from evaluateBlock import Evaluate_with_GPU
from scanPredctions import Scan
from postScan import Postprocess,annotatePAS
#
#from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
import datetime
import gc
import pickle
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
    parser.add_argument('--use_GPU',default = 'no', help='whether use GPU, can only be run with single thread')
    
    global keep_temp
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
    use_GPU = argv.use_GPU
    
    return out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread,block_length,use_GPU
    
def run_single_block(input_list):
    if(len(input_list) != 6):
        print('Error: input list is not correct',input_list)
        return
    else:
        _block_num,_start_line,_end_line,use_GPU,_start_pos,_end_pos = input_list
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
    if '.h5' in model:
        print('### Using DeepPASS from SCAPTURE')
        pred_out = EvaluateDeepPASS(chromosome,strand,block,model,rst,window,keep_temp)
    else:
        if use_GPU == 'yes':
            pred_out = Evaluate_with_GPU(chromosome,strand,block,model,rst,window,keep_temp)
            #pred_out = Evaluate(chromosome,strand,block,model,rst,window,keep_temp)
        else:
            os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
            pred_out = Evaluate(chromosome,strand,block,model,rst,window,keep_temp)
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
    out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread,block_length,use_GPU = args()
    start_time = datetime.datetime.now()
    # prepare output directory 
    if not os.path.exists(out_dir):
        print('### Creating output directory:{}'.format(out_dir))
        os.makedirs(out_dir)
    if fa_file is None:
        print('### error: the reference genome fasta file is required')
        sys.exit(1)
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
        #for chromosome in blocks_input_list:
        for chromosome in all_lines:
            _chr_start_time = datetime.datetime.now()
            print('### Processing {} in {} strand'.format(chromosome,strand))
            lines = all_lines[chromosome]
            o_chr_blocks = get_block_position(lines,window,block_length,keep_temp)
            # add use_GPU options to blocks 
            chr_blocks = []
            for oe in o_chr_blocks:
                oe[3] = use_GPU
                chr_blocks.append(oe)
            #chr_blocks[3] = use_GPU
            log.write('Blocks inf:{}_{}\t{}\n'.format(chromosome,strand,'\t'.join(str(e) for e in chr_blocks)))
            print('### {} blocks in {}_{} strand'.format(len(chr_blocks),chromosome,strand))
            if use_GPU == 'yes':
                print('### Run with GPU with single thread'.format(thread))
                thread = 1
            else:
                print('### {} threads would be run in paralle with CPU'.format(thread))
            # determin the number of blocks and defined thread 
            if len(chr_blocks) > thread:
                thread = len(chr_blocks)
            #    
            with ProcessPoolExecutor(max_workers=thread) as executor:
                temp_out = out_dir + '/temp.' + name + '.' + chromosome + '.' + strand + '.pickle'
                if os.path.exists(temp_out) and os.path.getsize(temp_out) > 0:
                    print('### Resuming the process by loading data from the temporary file {}'.format(temp_out))
                    chr_anno_pas = pickle.load(open(temp_out,'rb'))
                else:
                    #chr_out_pas = executor.map(run_single_block,chr_blocks,chunksize = 2)
                    chr_out_pas = executor.map(run_single_block,chr_blocks)
                    chr_anno_pas = annotatePAS(DB_file,chr_out_pas,chromosome,strand)
                    if keep_temp == 'yes':
                        print('### Saving the results to from the temporary file {}'.format(temp_out))
                        pickle.dump(chr_anno_pas, open(temp_out,"wb"))    
                pas_out_list += chr_anno_pas
            #for each_block in chr_blocks:
            #    eachout_pas = run_single_block(each_block)
            #    pas_out_list += eachout_pas
            _chr_end_time = datetime.datetime.now()
            print('### Process {}_{} used time: {}'.format(chromosome,strand,_chr_end_time - _chr_start_time))
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
