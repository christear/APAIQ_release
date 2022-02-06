#!/usr/bin/python

import sys,os
import argparse
import glob
from generate_windows import Generate_windows,Get_block_position,split_chr_bedGraph2
from evaluate import Evaluate
from scanTranscriptome_forward import Scan_Forward
from scanTranscriptome_reverse import Scan_Backward
from postprocess import Postprocess
#
from multiprocessing import Pool
import datetime
import logging as log
import gc
#
def get_ref(fa_file):
    ref = dict()
    fa = open(fa_file,'r')
    skip = False
    chro = ''
    seq  = ''
    for line in fa.readlines():
        line = line.rstrip('\n')
        if('>' in line):
            chro = line.split(' ')[0]
            chro = chro[1:]
            if('chr' not in chro):
                chro = 'chr'+chro

            if(len(chro)>6 or 'M' in chro):
                skip = True
                continue
            else:
                skip = False
        else:
            if(skip):
                continue
            if chro in ref.keys():
                ref[chro] += line
            else:
                ref[chro] = line
    fa.close()
    return ref

def get_genome_sequence(fa_file):
    f = open(fa_file,"r")
    line = f.readline()
    line = line.rstrip('\n')
    f.close()
    return line

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', default='out_dir', help='out dir')
    parser.add_argument('--input_file', default=None, help='unstranded wig file')
    parser.add_argument('--input_plus', default=None, help='plus strand wig file')
    parser.add_argument('--input_minus', default=None, help='minus strand wig file')
    parser.add_argument('--fa_file',default=None,help='path to one line fa file')
    parser.add_argument('--keep_temp',default=None,help='if you want to keep temporary file, set to "yes"')
    parser.add_argument('--window', default=201, type=int, help='input length')
    parser.add_argument('--name', default='sample',help='sample name')
    parser.add_argument("--model", help="the model weights file", required=True)
    parser.add_argument("--RNASeqRCThreshold",default=0.05,type=float,help="RNA-Seq Coverage Threshold")
    parser.add_argument('--threshold', default=0,type=int,help='peak length lower than threshold will be fiter out')
    parser.add_argument('--penality', default=1,type=int,help='penality for prediction score lower than 0.5')
    parser.add_argument('--DB_file', default=None, help='polyA database file')
    parser.add_argument('--depth', default=1, type=float,help='total number of mapped reads( in millions)')
    parser.add_argument('--t', default = 30, type = int, help='number of thread')
    
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
    return out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread
    
def run_single_block(input_list):
    global ref
    #print(ref.keys())
    baseName,model,out_dir,rst,window,keep_temp,threshold,penality,DB_file,input_file,chromosome,strand,depth,start,end = input_list

    #log_dir = out_dir+'/log'
    #if not os.path.exists(log_dir):
    #    os.makedirs(log_dir)
    #log.basicConfig(filename='%s/%s.log'%(log_dir,baseName), level=log.INFO)
    print("Generating blocks ...%s %d %s"%(baseName,start,end))
    ####Generate sliding windlows
    gw_start_time = datetime.datetime.now()
    block = split_chr_bedGraph2(out_dir,input_file,chromosome,strand,window,ref[chromosome],depth,start,end)
    ww = open(baseName,'w')
    for a,b,c in block:
        ww.write('%s\t%s\t%s\n'%(a,b,c))
    ww.close()
    gw_end_time = datetime.datetime.now()
    print("Generate blocks used time: {}\n".format(gw_end_time - gw_start_time))

    print("Evaluating blocks ...%s %d %s"%(baseName,start,end))
    ev_start_time = datetime.datetime.now()
    Evaluate(baseName,block,model,out_dir,rst,window,keep_temp)
    
    del block #destroyed the block reference
    gc.collect() #manually run garbage collection process 

    ev_end_time = datetime.datetime.now()
    print("Evaluated blocks used time: {}\n".format(ev_end_time - ev_start_time))

    print("Postprocessing blocks ...%s %d %s"%(baseName,start,end))
    ps_start_time = datetime.datetime.now()
    Scan_Forward(baseName,threshold,penality,out_dir)
    Scan_Backward(baseName,threshold,penality,out_dir)
    if(keep_temp != 'yes'):
        predict_file = out_dir+'/predict/'+baseName+'.txt'
        os.system('rm %s'%predict_file)
    Postprocess(DB_file,baseName,threshold,penality,out_dir)
    ps_end_time = datetime.datetime.now()
    print("Postprocessed blocks used time: {}\n".format(ps_end_time - ps_start_time))

    if(keep_temp != 'yes'):
        forward_file=out_dir+"/maxSum/%s.forward.%d.%d.txt"%(baseName,threshold,penality)
        backward_file=out_dir+"/maxSum/%s.backward.%d.%d.txt"%(baseName,threshold,penality)
        os.system('rm %s %s'%(forward_file,backward_file))
    #print('Finished postprocessing...%s\n'%baseName)
    return [gw_end_time-gw_start_time,ev_end_time-ev_start_time,ps_end_time-ps_start_time]
            
            
#def main(out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread):
if __name__ == '__main__':
    out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread = args()
    if(out_dir[-1] == '/'):
        out_dir = out_dir[0:-1]
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_dir = out_dir+'/'+name
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if input_file is not None:
        input_plus = input_file
        input_minus = input_file

    files = (input_plus,input_minus)
    strands = ('+','-')

    print("Load reference")
    ref = dict()
    #ref = get_ref(fa_file)#
    for root, ds, fs in os.walk(fa_file):
        #print(root,ds,fs)
        for f in fs:
            fil = root+'/'+f
            chro_seq = get_genome_sequence(fil)
            chro = f.split('.')[1]
            ref[chro] = chro_seq

    print("Finished Load reference")

    log = open('%s/%s.log'%(out_dir,name),'w')
    for i in range(2):
        input_file = files[i]
        strand = strands[i]
        print("Processing %s strand"%strand)
        blocks = Get_block_position(out_dir,input_file,strand,window,1e6)
        block_input_list = []
        for chromosome,strand,block_num,start,end in blocks:
            baseName = '%s.%s_%s_%s'%(name,chromosome,strand,block_num)
            print('%s\t%d\t%d'%(baseName,start,end))
            block_input_list.append([baseName,model,out_dir,rst,window,keep_temp,threshold,penality,DB_file,input_file,chromosome,strand,depth,start,end])
        print("Predicting results ...")
        pred_start_time = datetime.datetime.now()

        #block_out_indic = []




        with Pool(thread) as p:
            #p.map(run_single_block,block_input_list)
            time_lists = p.map(run_single_block,block_input_list)
            #p.close()
            p.terminate()
            p.join()
            
            for i,input_list in enumerate(block_input_list):
                baseName = input_list[0]
                gw_time,ev_time,ps_time = time_lists[i]
                log.write('%s\t%s\t%s\t%s\n'%(baseName,str(gw_time),str(ev_time),str(ps_time)))
            #    print('%s\t%s\t%s\t%s\n'%(baseName,str(gw_time),str(ev_time),str(ps_time)))

        pred_end_time = datetime.datetime.now()
        print("Prediction used time: {}".format(pred_end_time - pred_start_time))

    log.close()
    out_file = '%s/%s.predicted.txt' %(out_dir,name)
    ww = open(out_file,'w')
    if(DB_file is not None): 
        ww.write('predicted_pasid\tdb_pasid\tdb_diff\tscore\n')
    else:
        ww.write('predicted_pasid\tscore\n')
    ww.close()
    os.system('cat %s/maxSum/*bidirection* >>%s'%(out_dir,out_file))
    if(keep_temp != 'yes'):
        os.system('rm -rf  %s/predict %s/maxSum'%(out_dir,out_dir))
        
    print("Job Done!")
    
#if __name__ == '__main__':
#    main(*args())
