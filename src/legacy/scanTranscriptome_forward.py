#!/usr/bin/env python -w
#
#
#Update 2021/04/28 Check peak

import argparse
import os


def maxSum(file,threshold,penality,out):

    OUT=open(out,'w')
    OUT.write('pas_id\tmaxPoint\tmaxPos\tstart\tend\n')
    predict = dict()
    coor2pas = dict()
    with open(file,'r') as f:
        for line in f:
            pas_id,score =line.rstrip('\n').split('\t')
            score = float(score)
            chromosome,coor,strand = pas_id.split(':')
            coor = int(coor)
            predict[coor] = score
            coor2pas[coor] = pas_id
    
    predict[100000000000] = 1
    coor2pas[100000000000] = 'chr'
    coor_list = list(coor2pas.keys())
    coor_list.sort()
    sum=0
    maxPos=0   ## position of peak
    maxPoint=0 ## score of peak
    start=0       ## peak start
    end = 0    ## peak end
    peak = 0 ## peak coor
    peak_score = 0  ##peak score
    for coor in coor_list:
        score = predict[coor]
        if(coor-end>1):
            if(maxPoint>threshold and sum>0):
                #newpas_id = chromosome+":"+str(maxPos)+":"+strand
                newpas_id = coor2pas[maxPos]
                #start += 1
                #OUT.write('%s\t%d\t%d\t%d\t%.3f\t%.3f\n'%(newpas_id,start,end,length,maxPoint,area))
                OUT.write('%s\t%.3f\t%d\t%d\t%d\t%d\n'%(newpas_id,maxPoint,maxPos,start,end,peak))

            start = coor
            end   = coor
            if(score>0.5):
                maxPos = coor
                maxPoint = score
                peak_score = score
                sum = score
            else:
                maxPos = coor
                maxPoint = 0
                peak_score = 0
                sum  = 0

        elif(score < 0.5):
            sum -= penality
            if(sum <= 0):
                if(maxPoint > threshold):
                    #newpas_id = chromosome+":"+str(maxPos)+":"+strand
                    #OUT.write('%s\t%d\t%d\t%d\t%.3f\t%.3f\n'%(newpas_id,start,end,length,maxPoint,area))
                    newpas_id = coor2pas[maxPos]
                    OUT.write('%s\t%.3f\t%d\t%d\t%d\t%d\n'%(newpas_id,maxPoint,maxPos,start,end,peak))
                start = coor
                sum=0
                maxPoint = 0
                peak_score = 0
            end = coor
        else:
            sum += score
            if(peak_score < score):
                peak_score = score
                peak = coor
            if(maxPoint < sum):
                maxPoint = sum
                maxPos   = coor
            if(sum<1):
                start = coor
                maxPoint = sum
                maxPos   = coor
            end=coor

    OUT.close()
    return 0


def args():
    ### Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', default=None, help='out dir')
    parser.add_argument('--baseName', help='baseName')
    parser.add_argument('--threshold', default=0,type=int,help='peak length lower than threshold will be fiter out')
    parser.add_argument('--penality', default=1,type=int,help='penality for prediction score lower than 0.5')
    args = parser.parse_args()

    baseName = args.baseName
    threshold = args.threshold
    penality = args.penality
    out_dir=args.out_dir
    return baseName,threshold,penality,out_dir


def Scan_Forward(baseName,threshold,penality,out_dir):
    if(out_dir[-1] == '/'):
        out_dir = out_dir[0:-1]
    new_dir = out_dir+'/maxSum'
    if not os.path.exists(new_dir):
        os.makedirs(new_dir) 

    predict=out_dir+'/predict/'+baseName+'.txt'
    print("Start forward scaning %s"%predict)
    out=out_dir+"/maxSum/%s.forward.%d.%d.txt"%(baseName,threshold,penality)
    maxSum(predict,threshold,penality,out)
    print("End forward scaning %s\n"%predict)
    return 0

if __name__ == "__main__":
    Scan_Forward(*args())
