# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 13:27:57 2022




@author: Braden
"""

import MDAnalysis as mda
from MDAnalysis.analysis import distances
import pandas as pd
import numpy as np
import sys, glob, os
from matplotlib import pyplot as plt
import statistics

print('for the purposes of finding the angle and distances in a potential amide-pi bond')
print('these are the args: {}'.format(sys.argv[0]))
os.chdir(sys.argv[1])
dir_out = sys.argv[1]
tpr = 'tpr-noh.pdb'
trj = 'trj-noh.xtc'
print('this script also takes a fname descriptor in arg 3')
fname = 'xylanase'
temp = sys.argv[2]
tset_num = sys.argv[3]
uni = mda.Universe(tpr, trj)
print('Universe made...')
"""
make the hardcoded indexes to send frame results to
we are using the center of mass for the ring system of N13F and the amide
region center of mass of Y17/F18
"""
traj = uni.trajectory
amide = [] # center of mass positions
ring = [] # center of mass positions
frames = [] # index for trajectory reference
tprs = [] # index of data frame of xyz and atomids
anchor = [] # pos of CG for angle calculation

print('Gathering center of masses from trajectory...')

for i, frame_open in enumerate(traj):
    "get the ring positions"
    a,b,c,d,e,f = frame_open.positions[[78,79,80,81,82,83]] # ring atoms
    ring_skel= np.zeros([6,3])
    zone= [a,b,c,d,e,f] # make an iterable list for assigning to skeleton array
    for u in range(0,6):
        ring_skel[u]=zone[u]
    ring_pos = np.average(ring_skel[:,:3], axis=0, weights=ring_skel[:,:3])
    ring.append(ring_pos)
    "get the amide positions"
    carb, oxy, nit = frame_open.positions[[117,118,119]] #index -1, start at 0
    noc_skel = np.zeros([3,3])
    noc_skel[0] = nit; noc_skel[1]=oxy; noc_skel[2] = carb
    amide_pos = np.average(noc_skel[:,:3], axis = 0, weights = noc_skel[:,:3])
    amide.append(amide_pos)
    "set the CG of the ring as the third point for angle calc"
    cg_pos = zone[0]
    anchor.append(cg_pos)
    "make dataframes and pass them to a list"
    tpr_trj = uni.trajectory[i]
    tpr_df=pd.DataFrame()
    tpr_df['atomno']=uni.atoms.ids
    tpr_df[['x','y','z']] = tpr_trj.positions
    tprs.append(tpr_df)
    frames.append(i)

print('Center of massses calculated. Finding angles and distances...')
"Make a dataframe of the 3-point coordinate systems"
system = pd.DataFrame(columns = ['distance', 'angle', 'cg',
                                 'ring','amide'])
system['cg'] = anchor; system['ring'] = ring; system['amide'] = amide

for index, row in system.iterrows():
    "get equalization for angle calculation"
    x = row['cg']; y = row['ring']; z = row['amide']
    xy = x-y
    zy = z-y
    cosine_angle = np.dot(xy, zy) / (np.linalg.norm(xy) * np.linalg.norm(zy))
    angle_pre = np.arccos(cosine_angle)
    angle=np.degrees(angle_pre)
    system.at[index,'angle'] = angle
    "get distances"
    distance = distances.distance_array(row['ring'],row['amide'])
    distance = distance[0][0]
    system.at[index, 'distance'] = distance
avg_d = np.average(system['distance']); var_d = np.var(system['distance']); std_dis = np.std(system['distance'])
avg_ang = np.average(system['angle']); var_ang = np.var(system['angle']); std_ang = np.std(system['angle'])
print('the average distance of the system is: {}'.format(avg_d))
print('the variance for distance is: {}'.format(var_d))
print('the average angle of the system is: {}'.format(avg_ang))
print('The variation in angle is: {}'.format(var_ang))
print('the standard deviation in angle is: {}'.format(std_ang))
print('the standard deviation in distance is: {}'.format(std_dis))
in_dist_pocket = system['distance'][system['distance']<=4.5]
perc_dist = (len(in_dist_pocket)/len(system['distance'])) *100
print('percent of sim distance less than 4 angstroms is {}%'.format(int(perc_dist)))
output=pd.DataFrame(system); output.to_csv(fname+'-pi-amide_RESULTS.csv', header = True,
                                          index = None, sep = ',')
original_stdout = sys.stdout # save a reference to the original stdout
with open(fname + '-pi-amide_dat.txt', 'w') as f:
    sys.stdout = f # switch the output to a file.
    print('the average distance of the system is: {}'.format(avg_d))
    print('the variance for distance is: {}'.format(var_d))
    print('the average angle of the system is: {}'.format(avg_ang))
    print('The variation in angle is: {}'.format(var_ang))
    print('percent of sim distance less than 4 angstroms is {}%'.format(int(perc_dist)))
    print('the standard deviation in angle is: {}'.format(std_ang))
    print('the standard deviation in distance is: {}'.format(std_dis))
    print(f'{avg_d}  {avg_ang}  {std_dis}  {std_ang}  {perc_dist}  {var_d}  {var_ang}')
    sys.stdout = original_stdout
