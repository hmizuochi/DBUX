import rasterio; import sys; import os; import math
import numpy as np; import matplotlib.pyplot as plt

#NVALUE=-20000; TPRANGE=10000; TNRANGE=-10000; SPRANGE=10000; SNRANGE=-10000
#PAIRSIZE=9; PREDSIZE=366; COL=180; ROW=180; LUT_MWSIZE=1; STEP=420; MAP_MWSIZE=3

pred_listname='../input/predlist.txt'
s_listname='../input/spatial_pairlist.txt'
t_listname='../input/temporal_pairlist.txt'

argvs=sys.argv; argc=len(argvs)
if(argc != 13):
    print('12 arguments are required: DBUX.py NVALUE TPRANGE TNRANGE SPRANGE SNRANGE PAIRSIZE PREDSIZE COL ROW LUT_MWSIZE STEP MAP_MWSIZE')
    quit()
NVALUE=int(argvs[1]); TPRANGE=int(argvs[2]); TNRANGE=int(argvs[3])
SPRANGE=int(argvs[4]); SNRANGE=int(argvs[5]); PAIRSIZE=int(argvs[6])
PREDSIZE=int(argvs[7]); COL=int(argvs[8]); ROW=int(argvs[9])
LUT_MWSIZE=int(argvs[10]); STEP=int(argvs[11]); MAP_MWSIZE=int(argvs[12])

def CropWindow(x,y,t,data):
    window=np.ones((SW*2+1,SW*2+1,TW*2+1))*NV
    for i in range(-SW,SW+1):
        for j in range(-SW,SW+1):
            for k in range(-TW,TW+1):
                if((x+i)>=0 and (x+i)<MAPCOL and (y+j)>=0 and (y+j)<MAPROW and (t+k)>=0 and (t+k)<MAPNUM):
                    window[j,i,k]=data[y+j][x+i][t+k]
    return window

def ReadST(listname):
    with open(listname,mode='r') as fp1:
        filename=fp1.readline().replace('\n','').replace('\r','')
        with open('../input/'+filename,mode='rb') as fp2: #map is raw binary with int16
            pair=np.fromfile(fp2,np.int16,-1).reshape(ROW,COL,1) #prototype reading, y[0:ROW-1],x[0:COL-1],t[0]
        for k in range(2,PAIRSIZE+1):
            filename=fp1.readline().replace('\n','').replace('\r','')
            with open('../input/'+filename,mode='rb') as fp2:
                tmp=np.fromfile(fp2,np.int16,-1).reshape(ROW,COL,1)
                pair=np.concatenate((pair,tmp),axis=2) #ROW,COL,PAIRSIZE
    return pair

def SearchMatchup(t_x,t_y,s_x,s_y,t_pair,s_pair,LEVEL,count,UPPERLEVEL): #t_y,t_x: pixel of interest s_y,s_x: searching matchup
    if(LEVEL==0):
    #s_image=np.where((TNRANGE<=t_pair[:,:,count]) & (t_pair[:,:,count]<=tmin_input+(LEVEL+1)*STEP),s_pair[:,:,count],NVALUE)
        if(TNRANGE<=t_pair[t_y,t_x,count] and t_pair[t_y,t_x,count]<=tmin_input[t_y,t_x]+(LEVEL+1)*STEP):
            s_image=s_pair[s_y,s_x,count]
        else:
            s_image=NVALUE
    elif(LEVEL==UPPERLEVEL[t_y,t_x]):
        #s_image=np.where((tmin_input+LEVEL*STEP<t_pair[:,:,count]) & (t_pair[:,:,count]<=TPRANGE),s_pair[:,:,count],NVALUE)
        if(tmin_input[t_y,t_x]+LEVEL*STEP<t_pair[t_y,t_x,count] and t_pair[t_y,t_x,count]<=TPRANGE):
            s_image=s_pair[s_y,s_x,count]
        else:
            s_image=NVALUE
    else:
        #s_image=np.where((tmin_input+LEVEL*STEP<t_pair[:,:,count]) & (t_pair[:,:,count]<=tmin_input+(LEVEL+1)*STEP),s_pair[:,:,count],NVALUE)
        if(tmin_input[t_y,t_x]+LEVEL*STEP<t_pair[t_y,t_x,count] and t_pair[t_y,t_x,count]<=tmin_input[t_y,t_x]+(LEVEL+1)*STEP):
            s_image=s_pair[s_y,s_x,count]
        else:
            s_image=NVALUE
    return s_image

def Gaussian(x,y):
    x_offset=x+SW; y_offset=y+SW; dim=MAP_MWSIZE-1
    Cx=math.factorial(dim)/(math.factorial(x_offset)*math.factorial(dim-x_offset))
    Cy=math.factorial(dim)/(math.factorial(y_offset)*math.factorial(dim-y_offset))
    #total=(2**dim)**2
    #return 1.0*Cx*Cy/total
    return 1.0*Cx*Cy

def GenLUT(tmin_input,tmax_input,t_pair,s_pair,MAXLEVEL,SW):
    #-----------initialization----------------
    lookup=np.zeros((ROW,COL,MAXLEVEL+1))
    second=np.zeros((ROW,COL,MAXLEVEL+1))
    snum=np.zeros((ROW,COL,MAXLEVEL+1),np.int)
    uncertainty=np.ones((ROW,COL,MAXLEVEL+1))*NVALUE
    #-----------------------------------------
    UPPERLEVEL=((tmax_input-tmin_input)/STEP).astype(np.int16) #ROW,COL

    print('Initial estimation of uncertainty...')
    for x in range(COL):            
        for y in range(ROW):
            for LEVEL in range(UPPERLEVEL[y,x]+1):
                for count in range(PAIRSIZE):
                    matchup=NVALUE
                    matchup=SearchMatchup(x,y,x,y,t_pair,s_pair,LEVEL,count,UPPERLEVEL)
                    if(SNRANGE<=matchup and matchup<=SPRANGE):
                        snum[y,x,LEVEL]=snum[y,x,LEVEL]+1
                        lookup[y,x,LEVEL]=lookup[y,x,LEVEL]+matchup
                        second[y,x,LEVEL]=second[y,x,LEVEL]+matchup**2
            for LEVEL in range(UPPERLEVEL[y,x]+1,MAXLEVEL+1):
                lookup[y,x,LEVEL]=lookup[y,x,UPPERLEVEL[y,x]]
                snum[y,x,LEVEL]=snum[y,x,UPPERLEVEL[y,x]]
                second[y,x,LEVEL]=second[y,x,UPPERLEVEL[y,x]]
    with np.errstate(divide='ignore',invalid='ignore'):
        lookup=np.where(snum==0,NVALUE,(1.0*lookup/snum).astype(int))
        second=np.where((snum<=1) |(lookup==NVALUE), NVALUE, ((1.0*snum/(snum-1))*(second/snum-lookup**2))) #unbiased variance

    #-------re-initialization------------------
    lookup=np.zeros((ROW,COL,MAXLEVEL+1))
    snum=np.zeros((ROW,COL,MAXLEVEL+1))
    #-----------------------------------------0

    print('LUT generating with MAP MovingWindow size = %d ...' % MAP_MWSIZE)
    for x in range(COL):            
        if(x%(COL/5)==0):
            sys.stdout.write('\rprocessed pixels: %d percent' % (x/COL*100))
            sys.stdout.flush()
        for y in range(ROW):
            for LEVEL in range(UPPERLEVEL[y,x]+1):
                totalweight=0
                for count in range(PAIRSIZE):
                    for i in range(-SW,SW+1):
                        for j in range(-SW,SW+1):
                            matchup=NVALUE
                            if((x+i>=0 and (x+i)<COL and (y+j)>=0 and (y+j)<ROW)):
                                matchup=SearchMatchup(x,y,x+i,y+j,t_pair,s_pair,LEVEL,count,UPPERLEVEL)
                            if(SNRANGE<=matchup and matchup<=SPRANGE):
                                snum[y,x,LEVEL]=snum[y,x,LEVEL]+1
                                lookup[y,x,LEVEL]=lookup[y,x,LEVEL]+matchup*Gaussian(i,j)**2
                                totalweight=totalweight+Gaussian(i,j)**2
                if(totalweight==0):
                    lookup[y,x,LEVEL]=NVALUE
                    uncertainty[y,x,LEVEL]=NVALUE
                else:
                    lookup[y,x,LEVEL]=lookup[y,x,LEVEL]/totalweight
                    uncertainty[y,x,LEVEL]=second[y,x,LEVEL]*Gaussian(0,0)**2/totalweight #standard error^2 of estimated lookup value
            for LEVEL in range(UPPERLEVEL[y,x]+1,MAXLEVEL+1):
                lookup[y,x,LEVEL]=lookup[y,x,UPPERLEVEL[y,x]]
                snum[y,x,LEVEL]=snum[y,x,UPPERLEVEL[y,x]]
                uncertainty[y,x,LEVEL]=uncertainty[y,x,UPPERLEVEL[y,x]]
    return lookup.astype(int),snum,uncertainty.astype(int)

def AveLUT(lookup_input,snum_input,second_input,MAXLEVEL):
    #-----initialization-------------
    lookup_ave=np.zeros((ROW,COL,MAXLEVEL+1),np.int)
    second_ave=np.zeros((ROW,COL,MAXLEVEL+1),np.int)
    snum_ave=np.zeros((ROW,COL,MAXLEVEL+1),np.int16)
    #--------------------------------
    for LEVEL in range(MAXLEVEL+1):
        for x in range(COL):
            for y in range(ROW):
                if(int(LUT_MWSIZE/2)<=LEVEL and LEVEL<=MAXLEVEL-int(LUT_MWSIZE/2)):
                    count=0
                    for m in range(LUT_MWSIZE):
                        if(lookup_input[y,x,LEVEL+int(m-LUT_MWSIZE/2)]<SNRANGE or SPRANGE<lookup_input[y,x,LEVEL+int(m-LUT_MWSIZE/2)]):
                            count=count+1
                        else:
                            lookup_ave[y,x,LEVEL]=lookup_ave[y,x,LEVEL]+lookup_input[y,x,LEVEL+int(m-LUT_MWSIZE/2)]
                            snum_ave[y,x,LEVEL]=snum_ave[y,x,LEVEL]+snum_input[y,x,LEVEL+int(m-LUT_MWSIZE/2)]
                            second_ave[y,x,LEVEL]=second_ave[y,x,LEVEL]+second_input[y,x,LEVEL+int(m-LUT_MWSIZE/2)]
                    if(count==LUT_MWSIZE):
                        lookup_ave[y,x,LEVEL]=NVALUE
                        snum_ave[y,x,LEVEL]=0
                        second_ave[y,x,LEVEL]=NVALUE
                    else:
                        lookup_ave[y,x,LEVEL]=lookup_ave[y,x,LEVEL]/(LUT_MWSIZE-count)
                        second_ave[y,x,LEVEL]=second_ave[y,x,LEVEL]/(LUT_MWSIZE-count)**2
                else:
                    lookup_ave[y,x,LEVEL]=lookup_input[y,x,LEVEL]
                    snum_ave[y,x,LEVEL]=snum_input[y,x,LEVEL]
                    second_ave[y,x,LEVEL]=second_input[y,x,LEVEL]
    for LEVEL in range(MAXLEVEL+1):
        lookup_output_filename='../output/LOOKUP'+str(LEVEL)+'_ave.raw'
        with open(lookup_output_filename,mode='wb') as fp1:
            lookup_ave[:,:,LEVEL].astype('int16').tofile(fp1)
        lookup_output_filename='../output/LOOKUP'+str(LEVEL)+'_ave_rel.raw'
        with open(lookup_output_filename,mode='wb') as fp1:
            snum_ave[:,:,LEVEL].astype('int16').tofile(fp1)
        lookup_output_filename='../output/LOOKUP'+str(LEVEL)+'_ave_uncertainty.raw'
        with open(lookup_output_filename,mode='wb') as fp1:
            second_ave[:,:,LEVEL].astype('int16').tofile(fp1)
    return lookup_ave,snum_ave,second_ave

def PredDBUX(pred_listname,tmin_input,lookup,snum,second,MAXLEVEL):

    with open(pred_listname,mode='r') as fp1:
        for i in range(PREDSIZE):
            #------map initialization------------------
            image=np.ones((ROW,COL),np.int)*NVALUE
            snum_output=np.ones((ROW,COL),np.int)*NVALUE
            second_output=np.ones((ROW,COL),np.int)*NVALUE
            level_output=np.ones((ROW,COL),np.int16)*NVALUE
            #------------------------------------------

            if(i%(PREDSIZE/20)==0):
                sys.stdout.write('\rprocessing date: TOTAL %d / CURRENT %d' % (PREDSIZE,i))
                sys.stdout.flush()
            date=fp1.readline().replace('\n','').replace('\r','')
            if(os.path.exists('../input/spatial_'+date+'.raw')):
                with open('../input/spatial_'+date+'.raw',mode='rb') as fp2: #map is raw binary with int16
                    s_input=np.fromfile(fp2,np.int16,-1).reshape(ROW,COL) #[0:ROW-1],x[0:COL-1]
            else:
                s_input=np.ones((ROW,COL),np.int16)*NVALUE
            if(os.path.exists('../input/temporal_'+date+'.raw')):
                with open('../input/temporal_'+date+'.raw',mode='rb') as fp2:
                    t_input=np.fromfile(fp2,np.int16,-1).reshape(ROW,COL)
            else:
                t_input=np.ones((ROW,COL),np.int16)*NVALUE
    
            for x in range(COL):
                for y in range(ROW):
                    if(t_input[y,x]<TNRANGE or t_input[y,x]>TPRANGE):
                        image[y,x]=s_input[y,x] #no temporal base map.
                        snum_output[y,x]=-88
                    else:
                        if(s_input[y,x]<SNRANGE or s_input[y,x]>SPRANGE):
                            #exist temporal-base map, no spatial-base map -> gap-filling
                            for LEVEL in range(MAXLEVEL+1):
                                if(LEVEL==0):
                                    if(TNRANGE<=t_input[y,x] and t_input[y,x]<=tmin_input[y,x]+(LEVEL+1)*STEP):
                                        image[y,x]=lookup[y,x,LEVEL]
                                        snum_output[y,x]=snum[y,x,LEVEL]
                                        second_output[y,x]=second[y,x,LEVEL]
                                        level_output[y,x]=LEVEL
                                elif(LEVEL==MAXLEVEL):
                                    if(tmin_input[y,x]+STEP*LEVEL<t_input[y,x] and t_input[y,x]<=TPRANGE):
                                        image[y,x]=lookup[y,x,LEVEL]
                                        snum_output[y,x]=snum[y,x,LEVEL]
                                        second_output[y,x]=second[y,x,LEVEL]
                                        level_output[y,x]=LEVEL
                                else:
                                    if(tmin_input[y,x]+STEP*LEVEL<t_input[y,x] and t_input[y,x]<=tmin_input[y,x]+STEP*(LEVEL+1)):
                                        image[y,x]=lookup[y,x,LEVEL]
                                        snum_output[y,x]=snum[y,x,LEVEL]
                                        second_output[y,x]=second[y,x,LEVEL]
                                        level_output[y,x]=LEVEL
                        else:
                            #exist temporal-base map, exist spatial-base map -> spatial-base map
                            image[y,x]=s_input[y,x]
                            snum_output[y,x]=-99
                            second_output[y,x]=0
            #output gap-filled map.
            with open('../output/spatial_'+date+'_comp.raw',mode='wb') as fp2:
                image.astype('int16').tofile(fp2)
                #fig=plt.figure(figsize=(13,7))
                #im=plt.imshow(image,cmap="bwr")
                #fig.colorbar(im); figname='../output/spatial_'+date+'_comp.png'; plt.savefig(figname);
                #plt.close(fig)
            with open('../output/spatial_'+date+'_rel.raw',mode='wb') as fp2:
                snum_output.astype('int16').tofile(fp2)
            with open('../output/spatial_'+date+'_uncertainty.raw',mode='wb') as fp2:
                with np.errstate(invalid='ignore'):
                    second_output=np.where(second_output==NVALUE,NVALUE,(second_output**0.5).astype('int16'))
                second_output.astype('int16').tofile(fp2)
            with open('../output/spatial_'+date+'_level.raw',mode='wb') as fp2:
                level_output.astype('int16').tofile(fp2)
    print('\n')

#main

#----------initial parameter check and data load------------
if(MAP_MWSIZE%2 == 0 or LUT_MWSIZE%2 ==0):
    print('MAP_MWSIZE and LUT_MWSIZE must be odd numbers'); quit()
SW=int(MAP_MWSIZE/2)
if(TPRANGE<TNRANGE):
    print('TPRANGE must be greater than TNRANGE'); quit()
if(SPRANGE<SNRANGE):
    print('SPRANGE must be greater than SNRANGE'); quit()
if((TNRANGE<=NVALUE and NVALUE<=TPRANGE) or (SNRANGE<=NVALUE and NVALUE<=SPRANGE)):
    print('NVALUE must NOT be between NRANGE and PRANGE'); quit()

print('input data loading...')
s_pair=ReadST(s_listname); t_pair=ReadST(t_listname);  #ROW,COL,PAIRSIZE
#-----------------------------------------------------------

print('doStatistics...')
#-------------doStatistics--------------
tmin_input=np.nanmin(np.where(t_pair!=NVALUE,t_pair,np.nan),axis=2)
tmax_input=np.nanmax(np.where(t_pair!=NVALUE,t_pair,np.nan),axis=2)
if(np.isnan(tmin_input).sum()!=0 or np.isnan(tmax_input).sum()!=0):
    print('at least one historical pixel must be available in for each location')
    quit()
MAXLEVEL=int(((tmax_input-tmin_input)/STEP).max())
print('MAXLEVEL=%d' % MAXLEVEL)
#---------------------------------------

#---------matrix initialization----------
lookup=np.empty((ROW,COL,MAXLEVEL+1),np.int)
uncertainty=np.empty((ROW,COL,MAXLEVEL+1),np.int)
snum=np.empty((ROW,COL,MAXLEVEL+1),np.int)

s_image=np.ones((ROW,COL),np.int16)*NVALUE
#----------------------------------------

#generate LUT
lookup,snum,uncertainty=GenLUT(tmin_input,tmax_input,t_pair,s_pair,MAXLEVEL,SW)

#---------output lookup maps----------------
for LEVEL in range(MAXLEVEL+1):
    lookup_output_filename='../output/LOOKUP'+str(LEVEL)+'.raw'
    with open(lookup_output_filename,mode='wb') as fp1:
        lookup[:,:,LEVEL].astype('int16').tofile(fp1)
    lookup_output_filename='../output/LOOKUP'+str(LEVEL)+'_rel.raw'
    with open(lookup_output_filename,mode='wb') as fp1:
        snum[:,:,LEVEL].astype('int16').tofile(fp1)
    lookup_output_filename='../output/LOOKUP'+str(LEVEL)+'_uncertainty.raw'
    with open(lookup_output_filename,mode='wb') as fp1:
        uncertainty[:,:,LEVEL].astype('int16').tofile(fp1)
#-------------------------------------------

#---------switch depending on LUT moving average will be applyed or not----------
if(LUT_MWSIZE==1):
    print('LUT generated! No moving average for LUT.')
    #if you do not want to gapfill LUT, comment out from here to CO
    print('LUT gapfilling by lower LUT.')
    for LEVEL in range(1,MAXLEVEL+1):
        lookup[:,:,LEVEL]=np.where(lookup[:,:,LEVEL]==NVALUE,lookup[:,:,LEVEL-1],lookup[:,:,LEVEL])
        snum[:,:,LEVEL]=np.where(lookup[:,:,LEVEL]==NVALUE,snum[:,:,LEVEL-1],snum[:,:,LEVEL])
        uncertainty[:,:,LEVEL]=np.where(lookup[:,:,LEVEL]==NVALUE,uncertainty[:,:,LEVEL-1],uncertainty[:,:,LEVEL])
    #CO
    print('do prediction...')
    PredDBUX(pred_listname,tmin_input,lookup,snum,uncertainty,MAXLEVEL)
else:
    print('LUT generated! applying moving average with size %d' % LUT_MWSIZE)
    lookup_ave,snum_ave,uncertainty_ave=AveLUT(lookup,snum,uncertainty,MAXLEVEL)
    #if you do not want to gapfill LUT, comment out from here to CO
    print('LUT gapfilling by lower LUT.')
    for LEVEL in range(1,MAXLEVEL+1):
        lookup_ave[:,:,LEVEL]=np.where(lookup_ave[:,:,LEVEL]==NVALUE,lookup_ave[:,:,LEVEL-1],lookup_ave[:,:,LEVEL])
        snum_ave[:,:,LEVEL]=np.where(lookup_ave[:,:,LEVEL]==NVALUE,snum_ave[:,:,LEVEL-1],snum_ave[:,:,LEVEL])
        uncertainty_ave[:,:,LEVEL]=np.where(lookup_ave[:,:,LEVEL]==NVALUE,uncertainty_ave[:,:,LEVEL-1],uncertainty_ave[:,:,LEVEL])
    #CO
    print('do prediction...')
    PredDBUX(pred_listname,tmin_input,lookup_ave,snum_ave,uncertainty_ave,MAXLEVEL)
#---------------------------------------------------------------------------------

