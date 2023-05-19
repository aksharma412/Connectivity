# -*- coding: utf-8 -*-

import numpy as np
import h5py
import math
import os
import scipy.signal
from scipy.signal import argrelextrema
import matplotlib
import time

# program to plot the radial distribution of water molecules surorunding the solute
# this can be useful to validate our choice of the solvent shell

molecule_name = input("Enter the molecule name : ")
spectroscopy = input("Enter the spectroscopy technique : ")
snapshots = input("Enter the number of snapshots : ")
D = input("Enter the cutoff distance : ")
D = float(D)
snapshots = int(snapshots)

frames = 10000
flag = frames//int(snapshots)
temp_indx = flag
distribution = []
conf_list = []
final = []

# solvent identifiers
O = 349
H = 350
t1 = time.time()
# Set the file to be read.
for index in range(1, snapshots+1):
    with h5py.File(f"{molecule_name}_MD.h5", "r") as f:
        print(f"Obtaining coordinate information for frame {temp_indx} \n")
        atom_sym = []
        intmdt = f"frame{temp_indx}/"
        box = np.array(f[intmdt]['dims_box'])
        atom_number = np.array(f[intmdt]['indx'])
        coords = np.array(f[intmdt]['coord'])
        atom_type = np.array(f[intmdt]['at_type'])
        connect_array = np.array(f[intmdt]['connect'])
        for i in range(len(connect_array)):
            if connect_array[i][0] and connect_array[i][1] and connect_array[i][2] and connect_array[i][3]:
                atom_sym.append('C')
            elif connect_array[i][0] and connect_array[i][1]:
                atom_sym.append('O')
            elif connect_array[i][0]:
                atom_sym.append('H')
    
        atom_info = [atom_number, np.array(atom_sym), coords, atom_type]
        solute_info = [[],[],[]] # has all the solute information in order : number, symbol, x,y,z coords
        solvent_info = [[],[],[]]
        avg_coord = []
        
       
        for i in range(3):
            for j in range(len(atom_info[0])):
                if atom_info[-1][j] != O and atom_info[-1][j] != H:
                    solute_info[i].append(atom_info[i][j])
                else:
                    solvent_info[i].append(atom_info[i][j])
            solute_info[i] = np.array(solute_info[i])
            solvent_info[i] = np.array(solvent_info[i])
            
        
        print(f'Number of Solute Atoms: {len(solute_info[0])}')
        
        #calculate the average corodinates
        for i in range(3):
            avg_coord.append(np.mean(solute_info[2][:,i]))
        
        # Set conditions to center the solvent atoms around the solute.
        for i in range(3):
            for j in range(len(solvent_info[2])):
                if (solvent_info[2][j][i]-avg_coord[i]) > box[i]/2:
                    solvent_info[2][j][i] -= box[i]
                elif (solvent_info[2][j][i]-avg_coord[i]) < -box[i]/2:
                    solvent_info[2][j][i] += box[i]

        # initializing intermediate arrays for faster computation
        int_info = [[],[],[]]
        
        max_dist = 0
        for n in range(len(solute_info[0])):
            dist = 0
            for i in range(3):
                dist += (solute_info[2][n][i]- avg_coord[i])**2
            dist = np.sqrt(np.sum(dist))
            if dist >  max_dist:
                max_dist = dist
         
        Dist_thresh = max_dist + float(D)
        for p in range(3):
            for o in range(len(solvent_info[0])):
                solvent_dist = 0
                for i in range(3):
                    solvent_dist += (solvent_info[2][o][i]- avg_coord[i])**2
                solvent_dist = np.sqrt(np.sum(solvent_dist))
                if solvent_dist <= Dist_thresh:
                    int_info[p].append(solvent_info[p][o])
            int_info[p] = np.array(int_info[p])
        
        print(f"Number of Intermediate Solvent Atoms: {len(int_info[0])}")
        
        
        final_info = [[],[],[]]
        for i in range(3):
            final_info[i].extend(solute_info[i])
        
        for j in range(len(solute_info[0])):
            for k in range(len(int_info[0])):
                d = 0
                for i in range(3):
                    d += (solute_info[2][j][i]- int_info[2][k][i])**2
                d = np.sqrt(np.sum(d))
                if d<float(D):
                    no_duplicates = True
                    for l in range(len(final_info[0])):
                        if int_info[0][k] == final_info[0][l]:
                            no_duplicates = False
                            break
                    if no_duplicates:
                        for m in range(3):
                            final_info[m].append(int_info[m][k])
        print(f"Number of Atoms before solvent check : {len(final_info[0])}")
            
        # Confirming the presence of full solvent molecules
        
        solvent_complete = False
        while not solvent_complete:
            new_number = []
            con = []
            mod_connect = np.column_stack((atom_info[0], connect_array))
            temp_ctr = 0
            for p in range(len(final_info[0])):
                ln_num = final_info[0][p]
                for q in range(len(atom_info[0])):
                    if ln_num == atom_info[0][q]:
                        for g in range(len(connect_array)):
                            if mod_connect[g][0] == ln_num:
                                con.append(mod_connect[g][1:])
                                break
            con = np.array(con)
            connect = [con[:,0], con[:,1], con[:,2], con[:,3]]
            
            # Checking the status of atoms connectivity
            for q in range(len(final_info[0])):
                lb = [False, False, False, False]
                for r in range(len(connect[0])):
                    for y in range(4):
                        if connect[y][q] == final_info[0][r]:
                            lb[y] = True
                            break
                    for y in range(4):
                        if connect[y][q] == 0:
                            if y<4:
                                lb[y] = True
                            if y<3:
                                lb[y] = True
                            if y<2:
                                lb[y] = True
                            if y<1:
                                lb[y] = True
                for y in range(4):
                    if lb[y] == False:
                        no_duplicates = True
                        for t in range(len(new_number)):
                            #print(len(new_number))
                            if connect[y][q] == new_number[t]:
                                no_duplicates = False
                                break
                        if no_duplicates:
                            new_number.append(connect[y][q])
                            break
                        
            #Checking if solvent is complete.
            print("Number of New Atoms:", len(new_number))
            if not len(new_number):
                print("Solvent Complete.")
                solvent_complete = True
                break
            else:
                print("Solvent Incomplete. Appending and obtaining connectivity of new atoms.")
        
            # Appending new atoms to final atom arrays.
            for s in range(len(new_number)):
                ind = new_number[s] - len(solute_info[0]) -1
                ind = int(ind)
                for i in range(3):
                    final_info[i].append(solvent_info[i][ind])
            
        print("Final Number of Atoms:", len(final_info[0]))
        conf_list.append(len(final_info[0])-len(solute_info[0]))
        #D += 0.25
        print(D)
        
        final.append(conf_list)    
        
    temp_indx += flag
    
t2 = time.time()
print(t2-t1)

  