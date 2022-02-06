#!/usr/bin/python

import sys,os
import argparse
import glob
from generate_windows import Generate_windows
from evaluate import Evaluate
from scanTranscriptome_forward import Scan_Forward
from scanTranscriptome_reverse import Scan_Backward
from postprocess import Postprocess
#
from multiprocessing import Pool
import datetime
#
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
    data = input_list[0]
    model = input_list[1]
    out_dir = input_list[2]
    rst = input_list[3]
    window = input_list[4]
    keep_temp = input_list[5]
    threshold = input_list[6]
    penality = input_list[7]
    DB_file = input_list[8]
    if 'wig' not in data:
        baseName = data.split('/')[-1]
        Evaluate(model,out_dir,rst,window,baseName,keep_temp)
        Scan_Forward(baseName,threshold,penality,out_dir)
        Scan_Backward(baseName,threshold,penality,out_dir)
        if(keep_temp != 'yes'):
            predict_file = out_dir+'/predict/'+baseName+'.txt'
            os.system('rm %s'%predict_file)
        Postprocess(DB_file,baseName,threshold,penality,out_dir)
        if(keep_temp != 'yes'):
            forward_file=out_dir+"/maxSum/%s.forward.%d.%d.txt"%(baseName,threshold,penality)
            backward_file=out_dir+"/maxSum/%s.backward.%d.%d.txt"%(baseName,threshold,penality)
            os.system('rm %s %s'%(forward_file,backward_file))
        print('Finished postprocessing...\n')
            
            
def main(out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,model,rst,threshold,penality,DB_file,depth,thread):
    
    if(out_dir[-1] == '/'):
        out_dir = out_dir[0:-1]
        
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_dir = out_dir+'/'+name
	####Generate sliding windlows
    print("Generating blocks ...")
    gw_start_time = datetime.datetime.now()
    #Generate_windows(out_dir,input_file,input_plus,input_minus,fa_file,keep_temp,window,name,depth)
    gw_end_time = datetime.datetime.now()
    print("Gerate blocks used time: {}".format(gw_end_time - gw_start_time))
    data_dir = out_dir+'/data'
    data_files = glob.glob(data_dir+"/*")
    block_input_list = []
    for data in data_files:
        block_input_list.append([data,model,out_dir,rst,window,keep_temp,threshold,penality,DB_file])
    print("Predicting results ...")
    pred_start_time = datetime.datetime.now()
    with Pool(thread) as p:
        p.map(run_single_block,block_input_list)
    pred_end_time = datetime.datetime.now()
    print("Prediction used time: {}".format(pred_end_time - pred_start_time))
    out_file = '%s/%s.predicted.txt' %(out_dir,name)
    ww = open(out_file,'w')
    if(DB_file is not None): 
        ww.write('predicted_pasid\tdb_pasid\tdb_diff\tscore\n')
    else:
        ww.write('predicted_pasid\tscore\n')
    ww.close()
    os.system('cat %s/maxSum/*bidirection* >>%s'%(out_dir,out_file))
    if(keep_temp != 'yes'):
        os.system('rm -rf %s/data %s/predict %s/maxSum'%(out_dir,out_dir,out_dir))
        
    print("Job Done!")
    
    
if __name__ == '__main__':
	main(*args())
            
            
            
            