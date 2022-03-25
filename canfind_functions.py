#CANFind functions 
#Katie Fasbender 

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import matplotlib #necessary for Hyalite, don't ask questions you don't want the answers to! 
from dl import queryClient as qc
import healpy as hp
import numpy as np
from numpy import arange,array,ones,linalg
from astropy import utils, io
from astropy.io import fits
from astropy.table import Table, vstack, join, Column
from astropy.coordinates import SkyCoord, ICRS, Galactic
import math as mat
import itertools as it
import math as m
from statistics import mean
from scipy.stats import pearsonr
import numpy.ma as ma
from astropy.stats import sigma_clip,SigmaClip
from collections import Counter
import sys
from sklearn.cluster import DBSCAN
from sklearn import linear_model, datasets
from dl import queryClient as qc

#-----------------------------------------------------------------------------
# Functions 
#-----------------------------------------------------------------------------
def query_fix(query,dtypes): 
    ''' fix a query in csv format 
   Arguments:
       query = the query (str)
       dtypes = data types for each column queried (array-typr)
   '''
    bc=[]
    bd=query.splitlines()
    for j in bd:
        bc.append(j.split(",")) 
    return(Table(np.array(bc[1:]),names=np.array(bc[0]),dtype=dtypes)) #'f8' makes it a float64
#----------------------------------------------------------------------------
def removal(cluster,out_table):
    '''Remove a cluster if invalid
    Arguments:
        cluster = the cluster you want to analyze (astropy table)
        out_table = table with columns "measid", "cluster_label" (astropy table)
    '''
    for kt in range(0,len(cluster)): #for each point in the cluster
        clonck=out_table['measid']==cluster['measid'][kt] #using measid, find the same point in out_table
        it=clonck.tolist().index(True) #index of bad point
        out_table[it]['cluster_label']=-1 #give outlier label
#----------------------------------------------------------------------------
def labeling(label,out_table,members,pos):
    '''Give mmts a track label in out_table for either track formation method 
    Arguments:
        label = some number for the track label 
        out_table = table with columns "measid", "track_h", "track_p"
        members = track_members function outpit
        pos = "p" for pred_pos, "h" for hyp_pos
    '''
    for kt in range(0,len(members)): #for each point
        clonck=out_table['measid']==members[kt] #find matching point in out_table
        it=clonck.tolist().index(True) #index of bad point
        if pos=='h': #if hyp_pos was used to predict position,
            out_table[it]['track_h']=label #give appropriate track label 
        if pos=='p':  #if pred_pos was used to predict position,  
            out_table[it]['track_p']=label #give appropriate track label
#----------------------------------------------------------------------------
def time_dup_removal(cluster,out_table):
    '''Remove duplicate 'mjd's
    Arguments:
        cluster = the cluster you want to analyze
        out_table = table with columns "measid", "cluster_label", "mjd"
    ''' 
    my_list=np.array(cluster['mjd']) #list of mjds
    dups=[k for k,v in Counter(my_list).items() if v>1] #array of duplicate mjds
    for d in dups: #for each duplicate mjd,
        bad=cluster['mjd']==d #bad points = ones whose mjd matches a duplicate mjd
        for j in cluster[bad]['measid']: #for every measurement with a duplicate mjd,
            get_out=out_table['measid']==j #remove the measurement from the cluster
            indo=get_out.tolist().index(True)
            out_table[indo]['cluster_label']=-1
#----------------------------------------------------------------------------
def sig_clip_t(cluster,out_table):
    '''Sigma clip on the "mjd"s
    Arguments:
        cluster = the cluster you want to analyze
        out_table = table with columns "mjd"
    '''
    filtered_data=sigma_clip(cluster['mjd'],sigma=2) #masked array of data
    clip_mask=[] 
    for i in filtered_data: #for each data point in the cluster,
        if i is ma.masked: #if the data point is masked,
            clip_mask.append(True) #append a "True" to clip_mask so it can be removed
        else:
            clip_mask.append(False)        
    removal(cluster[clip_mask],out_table) #then remove the bad data points
#----------------------------------------------------------------------------
def peacc(Cluster,spacetime):
    '''Calculate pearson correlation coefficient of cluster 
    Arguments:
        spacetime = "s" if you want the PCC in RA/Dec
                  = "t" if you want the PCC in RA/mjd 
    '''
    x=Cluster['ra'] #RA coord.s of cluster members (X)
    if spacetime=='s': #if spacial,
        y=Cluster['dec'] #Dec coord.s of cluster members  (Y_space)
    if spacetime=='t': #if temporal,
        y=Cluster['mjd'] #time of measurement of cluster members (Y_time)
    pearson=pearsonr(x,y)[0] #calculate PCC of cluster 
    return(pearson)            
#----------------------------------------------------------------------------
def ranslap(cluster,out_table):
    '''RANSAC analysis to remove linearity outliers of cluster in ra-mjd space, 
    also calculates velocity of cluster (at this point, a tracklet) 
    Arguments:
        cluster = the cluster you want to analyze
        out_table = table with columns "ra", "dec", "mjd"
    ''' 
    ra=np.array(cluster['ra']) #RA coord.s of cluster members 
    dec=np.array(cluster['dec']) #Dec coords of cluster members
    Xra=np.reshape(ra,(len(ra),1)) #reshape RA coord.s array (X_ra)
    Xdec=np.reshape(dec,(len(dec),1)) #reshape Dec coord.s array  (X_dec)
    y=np.array(cluster['mjd']) #time of measurement of cluster members (Y)
    if len(np.unique(y))<3: #if there are fewer than 3 unique measurement times (fewer than 3 moving object detections)
        inlier_mask=np.zeros(len(ra), dtype=bool) #list of zeros
        outlier_mask=np.ones(len(ra), dtype=bool) #list of ones
        clu_t=cluster[inlier_mask] #time cluster inliers 
        no_t=cluster[outlier_mask] #time cluster outliers
        clu_s=clu_t #space cluster inliers
        no_s=no_t #space cluster outliers
        return(clu_t,no_t,0,0,clu_s,no_s) #does not RANSAC but gives all members "outlier" labels
    else: #if there are at least 3 independent moving object detections in the cluster,
        # Robustly fit linear model with RANSAC algorithm-----------------------------------------
        ransac_ra = linear_model.RANSACRegressor(residual_threshold=.001)
        ransac_ra.fit(Xra, y) #RANSAC on RA and MJD
        ransac_dec = linear_model.RANSACRegressor(residual_threshold=.001)
        ransac_dec.fit(Xdec, y) #RANSAC on Dec and MJD
        ransac_s = linear_model.RANSACRegressor(residual_threshold=.001)
        ransac_s.fit(Xra, Xdec) #RANSAC on RA and Dec
        # Predict data of estimated models--------------------------------------------------------
        line_Xra = np.reshape(np.arange(Xra.min(), Xra.max()+.0001,step=((Xra.max()+.0001)-Xra.min())/20),(20,1))
        line_y_ransac_ra = ransac_ra.predict(line_Xra) #line for RANSAC fit, ra     
        line_Xdec = np.reshape(np.arange(Xdec.min(), Xdec.max()+.0001,step=((Xdec.max()+.0001)-Xdec.min())/20),(20,1))
        line_y_ransac_dec = ransac_dec.predict(line_Xdec) #line for RANSAC fit, dec    
        xsra=np.concatenate(line_Xra) #x values of RA,MJD RANSAC line
        xsdec=np.concatenate(line_Xdec) #x values of Dec,MJD RANSAC line
        ysra=line_y_ransac_ra #y values of RA,MJD RANSAC line
        ysdec=line_y_ransac_dec #y values of Dec,MJD RANSAC line
        #---------------------------------------------------------------------------------------
        mra = (ysra[-1]-ysra[0])/(xsra[-1]-xsra[0]) # 1/slope of RA,MJD RANSAC (1/velocity in RA)
        mdec = (ysdec[-1]-ysdec[0])/(xsdec[-1]-xsdec[0]) # 1/slope of Dec,MJD RANSAC (1/velocity in Dec)
        #---------------------------------------------------------------------------------------
        inlier_mask_t = ransac_ra.inlier_mask_
        outlier_mask_t = np.logical_not(inlier_mask_t)
        clu_t=cluster[inlier_mask_t] #time cluster after outlier removal (cluster inliers)
        no_t=cluster[outlier_mask_t] #time cluster outliers
        #---------------------------------------------------------------------------------------
        inlier_mask_s = ransac_s.inlier_mask_
        outlier_mask_s = np.logical_not(inlier_mask_s)
        clu_s=cluster[inlier_mask_s] #space cluster after outlier removal (cluster inliers)
        no_s=cluster[outlier_mask_s] #spacecluster outliers
        #---------------------------------------------------------------------------------------
        if mra!=0 and mdec!=0: #if the slopes are not zero,
            return (clu_t,no_t,1/mra,1/mdec,clu_s,no_s)  #This is the velocity of the object that the tracklet (cluster) represents in ra & dec directions
        else: #if the slopes ARE zero,
            return(clu_t,no_t,0,0,clu_s,no_s) #return "0" as velocities 
#----------------------------------------------------------------------------
def validate_it(X,db,out_table,labels,min_ps,min_pt):
    '''Validate each cluster using functions "peacc" and "ranslap", removes the outliers
    Arguments:
        X = outliers from first dbscan
        db = second dbscan output
        labels = cluster labels from second dbscan output
        min_ps = PCC lower cutoff in RA,Dec 
        min_pt = PCC lower cutoff in RA,mjd
    '''
    for i in range(0,max(labels)+1): #for each cluster i, except the outliers
        clust=out_table['cluster_label']==i #define the cluster
        cluster=out_table[clust]
        if len(cluster)>1:   #if there's more than 1 member in the cluster,
            time_dup_removal(cluster,out_table) #remove duplicate mjds
            clusta=out_table['cluster_label']==i  #re-define cluster
            clustera=out_table[clusta]
            sig_clip_t(clustera,out_table) #sigma clip the mjds
            clustb=out_table['cluster_label']==i #re-define the cluster
            clusterb=out_table[clustb]
            rans=ranslap(clusterb,out_table) #ransac the cluster, rans[1] = time outliers, rans[5] = space outliers
        #Space
            if len(rans[5])>0: #if there are spacial outliers,
                removal(rans[5],out_table) #remove the outliers
            clustc=out_table['cluster_label']==i #re-define the cluster
            clusterc=out_table[clustc]
            if len(clusterc)>2:
                pp=peacc(clusterc,'s') #calculate PCC of cluster, spacially
                if abs(pp)<min_ps: #if spacial PCC is too low, get rid of cluster!
                    removal(clusterc,out_table) #gets rid of cluster (gives points "outlier" label)
                else: #if PCC is high enough,
        #Time
                    if len(rans[1])>0: #If there are outliers from the RANSAC in mjd/ra,
                        removal(rans[1],out_table)  #remove the outliers        
                    clustd=out_table['cluster_label']==i #re-define the cluster
                    clusterd=out_table[clustd]
                    if len(np.unique(clusterd['mjd']))>2: #if there are more than 2 unique mjds,
                        pt=peacc(clusterd,'t') #calculate PCC of cluster, temporally
                        if abs(pt)<min_pt: #if temporal PCC is too low, get rid of cluster!
                            removal(clusterd,out_table) #gets rid of cluster (gives points "outlier" label)
                    else: #if there are fewer than 2 unique mjds,
                        removal(clusterd,out_table) #gets rid of cluster (gives points "outlier" label)
        #Tracklets need to have 3 or more measurements. If so, give them their appropriate velocities.
        new=out_table['cluster_label']==i
        new_cluster=out_table[new] #new cluster
        if len(new_cluster)<3 or (rans[2]==0 and rans[3]==0): #if the length of the cluster is less than 3, or tracklet velocity=0, make 'em all outliers!
            removal(new_cluster,out_table) #give cluster "outlier" label
        else: #if the length of the cluster is greater than 3, give them the appropriate velocities 
            v_ra=rans[2] #define tracklet velocity in RA
            v_dec=rans[3] #define tracklet velocity in dec
            for plu in range(0,len(out_table)): #for every measurement
                if out_table['cluster_label'][plu]==i: #if the measurement corresponds to cluster i,
                    out_table['v_ra'][plu]=v_ra #give it the tracklet velocity in RA
                    out_table['v_dec'][plu]=v_dec #and give it the tracklet velocity in Dec 
#----------------------------------------------------------------------------
def pred_pos(table,t=[],to=[]): 
    '''Predict tracklet member positions using corresponding tracklet velocities
    Arguments:
        table = out-table
        t = time of desired pred_pos.(see below for specifications)
    '''
    mid=[] #empty array for measurement id 
    lab=[] #empty array for cluster label
    ra_ps=[] #empty array for predicted ra
    dec_ps=[] #empty array for predicted dec
    for p in table: #for every entry (measurement, point) in out-table
        if t==0:  #do this if you just want to advance all the times in out-table by amount "to"
            too=p['mjd']+to #in this case you'd set t=0, and to=however much time you want to advance by
        else: #this is when you want to see predicted positions at some specific mjd "t", and set to=0
            too=t
        if p['cluster_label']!=-1: #if cluster label is NOT outlier (-1),
            ra_pred=p['v_ra']*(too-p['mjd'])+p['ra'] #calculate predicted RA
            dec_pred=p['v_dec']*(too-p['mjd'])+p['dec'] #calculate predicted DEC
            dev=abs(pow((pow(p['v_ra'],2)+pow(p['v_dec'],2)),0.5)*7) #distance a point in a cluster (tracklet) could have traveled in 7 days 
            devv=abs(pow(pow((ra_pred-p['ra']),2)+pow((dec_pred-p['dec']),2),0.5)) #distance actually traveled
            if devv<dev: #if the distance traveled is less than the distance the tracklet COULD have traveled in 7 days,
                mid.append(p['measid'])  #add the point's measurement id to mid
                lab.append(p['cluster_label']) #add cluster label to lab
                ra_ps.append(ra_pred) #add predicted RA to ra_ps
                dec_ps.append(dec_pred) #add predicted Dec to dec_ps
    return(mid,ra_ps,dec_ps,lab) #returns measurement id's, predicted RA & Dec, and cluster label
#----------------------------------------------------------------------------
def hyp_pos(table,cluster,t=[],to=[]): 
    '''Predict hypothetical point positions using one tracklet's velocity 
    Arguments:
        table = out_table
        t = time of desired pred_pos.(see below for specifications) 
        cluster = cluster label (#)
    '''
    mid=[] #empty array for measurement id 
    lab=[] #empty array for cluster label
    ra_ps=[] #empty array for predicted RA
    dec_ps=[] #empty array for predicted Dec
    clustie=table['cluster_label']==cluster 
    tab=table[clustie] #define the cluster whose velocities you're using 
    vevra=tab['v_ra'][0] #cluster velocity in RA
    vevdec=tab['v_dec'][0] #cluster velocity in Dec
    dev=abs(pow((pow(vevra,2)+pow(vevdec,2)),0.5)*7) #distance a point could have traveled in 7 days at cluster (tracklet) velocity 
    for p in table: #for every entry (measurement, point) in out-table
        if t==0:  #do this if you just want to advance all the times in out-table by amount "to"
            too=p['mjd']+to #in this case you'd set t=0, and to=however much time you want to advance by
        else: #this is when you want to see predicted positions at some specific mjd "t", and set to=0
            too=t
        ra_pred=vevra*(too-p['mjd'])+p['ra'] #calculate predicted RA
        dec_pred=vevdec*(too-p['mjd'])+p['dec'] #calculate predicted DEC
        devv=abs(pow(pow((ra_pred-p['ra']),2)+pow((dec_pred-p['dec']),2),0.5)) #distance actually traveled
        if devv<dev: #if the distance traveled is less than the point COULD have traveled at cluster's velocity in 7 days, 
            ra_ps.append(ra_pred) #add predicted RA to ra_ps
            dec_ps.append(dec_pred) #add predicted Dec to dec_ps
            mid.append(p['measid'])  #add the measurement id to mid
            lab.append(p['cluster_label']) #add cluster label to lab
    return(mid,ra_ps,dec_ps,lab) #returns measurement id's, hypothetical RA & Dec, and cluster label
#----------------------------------------------------------------------------
def track_members(cluster,pos,out_table): 
    '''Find new track members after projecting measurement positions to common time using pred_pos or hyp_pos
    Arguments:
        cluster = the cluster whose velocity you used to project point positions. 
        pos = the hyp_pos or pred_pos you ran using cluster's velocity
    '''
    cl_id=[] #empty array for measurement id
    clo=out_table['cluster_label']==cluster 
    cc=out_table[clo] #define the cluster whose velocity you used to project point positions 
    center_ra=cc['ra'][0] #RA of the cluster at the time you projected the positions to 
    center_dec=cc['dec'][0] #Dec of the cluster at the time you projected the positions to 
    for i in range(0,len(pos[1])): #for every point projected,
        position=pow(pow((pos[1][i]-center_ra),2)+pow((pos[2][i]-center_dec),2),0.5) #find the distance from RA,Dec of cluster
        if position<0.005: #if that distance is less than my decided max distance,
            #cl_ra.append(hypo_pos[1][i])
            #cl_dec.append(hypo_pos[2][i])
            cl_id.append(pos[0][i]) #append the point's measurement id to the output list
    return(cl_id) #returns track members' measurement id's