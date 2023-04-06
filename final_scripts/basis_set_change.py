# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 09:10:41 2023

@author: Aparna K
"""

import os
import numpy as np
import time

# cwd=os.getcwd()

# molecule_name = input("Enter the molecule name : ")
# desired_bs = input("Enter the new basis set  : ")
# dft_bs = "CAM-B3LYP/"

# t1 = time.time()

# number_of_conformers = len(os.listdir(f"{molecule_name}_MD"))
# print(f"Number of Conformers: {number_of_conformers}")
# print()

# conf = 1
# while conf <= number_of_conformers:
#     # Enters the directory with the input file.
#     os.chdir(f"{cwd}/{molecule_name}_MD/cmpd_{conf}")
    
#     #replace the basis set with the desired basis set
#     with open("input.dat","r") as f:
#         data = f.readlines()
        
#     data[3] = data[3].replace(data[3][3:19], dft_bs+desired_bs + ' ')

    
#     with open("input.dat","w") as f:
#         f.writelines(data)
#     conf += 1
#     os.chdir(f"{cwd}/{molecule_name}_MD")
    
# os.chdir(cwd)

# t2 = time.time()
# print(t2-t1)

import os
from multiprocessing import Pool

cwd = os.getcwd()
molecule_name = input("Enter the molecule name : ")
desired_bs = input("Enter the new basis set  : ")
dft_bs = "CAM-B3LYP/"

def modify_input_file(conf):
    # Enters the directory with the input file.
    os.chdir(f"{cwd}/{molecule_name}_MD/cmpd_{conf}")
    
    #replace the basis set with the desired basis set
    with open("input.dat","r") as f:
        data = f.readlines()
        
    data[3] = data[3].replace(data[3][3:19], dft_bs+desired_bs + ' ')
    
    with open("input.dat","w") as f:
        f.writelines(data)
    
    os.chdir(f"{cwd}/{molecule_name}_MD")
    

if __name__ == '__main__':
    t1 = time.time()
    number_of_conformers = len(os.listdir(f"{molecule_name}_MD"))
    print(f"Number of Conformers: {number_of_conformers}")
    print()
    
    # create a pool of worker processes
    pool = Pool()
    
    # apply the modification function to each conformer in parallel
    pool.map(modify_input_file, range(1, number_of_conformers+1))
    
    # close the worker processes
    pool.close()
    pool.join()
    
    os.chdir(cwd)
    
    t2 = time.time()
    print(t2-t1)



