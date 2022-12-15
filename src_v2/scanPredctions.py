#!/usr/bin/env python -w
#
#
import argparse
import os

# method should be forward or reverse 
def Scan(pred_out,threshold,penality,method):
    print("### Start {} scaning".format(method))
    scan_out = []
    predict = dict()
    coor2pas = dict()
    for pas_id,score in pred_out:
        chromosome,coor,strand = pas_id.split(':')
        coor = int(coor)
        predict[coor] = score
        coor2pas[coor] = pas_id
    if 'forward' in method:
        predict[100000000000] = 1
        coor2pas[100000000000] = 'chr'
        coor_list = list(coor2pas.keys())
        coor_list.sort()
        
        maxPos=0   ## position of peak
        start=0       ## peak start
        end = 0    ## peak end
        peak = 0 ## peak coor
    elif 'reverse' in method:
        predict[-1] = 1
        coor2pas[-1] = 'chr'
        coor_list = list(coor2pas.keys())
        coor_list.sort(reverse=True)
        
        maxPos= coor_list[0] ## position of peak
        start= coor_list[0]      ## peak start
        end = coor_list[0]    ## peak end
        peak = coor_list[0] ## peak coor
        
    else:
        print('### Error: {} has not been defined'.format(method))   
    sum=0
    maxPoint=0 ## score of peak
    peak_score = 0  ##peak score
    
    for coor in coor_list:
        score = predict[coor]
        if(abs(coor-end)>1):
            if(maxPoint>threshold and sum>0):
                newpas_id = coor2pas[maxPos]
                scan_out.append((newpas_id,maxPoint,maxPos,start,end,peak))

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
                    newpas_id = coor2pas[maxPos]
                    scan_out.append((newpas_id,maxPoint,maxPos,start,end,peak))
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
    print("### End {} scaning".format(method))
    return scan_out
    
def args():
    ### Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--pred_out', default=None, help='output list from evaluation')
    parser.add_argument('--threshold', default=0,type=int,help='peak length lower than threshold will be fiter out')
    parser.add_argument('--penality', default=1,type=int,help='penality for prediction score lower than 0.5')
    parser.add_argument('--method', default='forward',type=str,help='method for scanning, forward or reverse')
    args = parser.parse_args()
    
    pred_out = args.pred_out
    threshold = args.threshold
    penality = args.penality
    method = args.method
    return pred_out,threshold,penality,method

if __name__ == "__main__":
    scan_out = Scan(*args())
