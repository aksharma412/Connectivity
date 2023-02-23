
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 21:57:40 2023

@author: Aparna K
"""
import os
import shutil
import numpy as np
import math
import subprocess

# User Input prompts
molecule_name = input('Molecule Name: ')
spectroscopy = input('Spectroscopy: ')
basis_set = input('Basis Set: ')

if basis_set == 'Mixed':
    solute_basis_set = input('Solute Basis Set: ')
    solvent_basis_set = input('Solvent Basis Set: ')

snapshots = input('Number of Snapshots: ')
snapshots = int(snapshots)
D = input('Distance Threshold (Angstroms): ')

if spectroscopy != 'VCD' and spectroscopy != 'ROA':
    print(f"{spectroscopy} is not a valid Gaussian spectroscopy supported in this script.")
    print("Valid spectroscopies currently implemented are VCD and ROA.")
    exit(1)

if spectroscopy == 'ROA':
    frequency = input('Frequency of Incident Radiation (nm): ')

cwd = os.getcwd()

#_____________________________________________________________
# ==> Generate Data for Single Snapshot from MD Trajectory <==
#_____________________________________________________________

# Read the number of frames in the MD trajectory.

with open(f"{molecule_name}_MD.arc", "r") as f:
    lines = len(f.readlines())  #word count in bash
with open(f"{molecule_name}_MD.arc", "r") as f:
    atoms = int(f.readline().split()[0])  #head in bash
    
    

frames = (lines // (atoms + 2))
print(f"Frames: {frames}")
print()

# Obtains the specific frames associated with the number of snapshots requested.
flag = frames//int(snapshots)
counter = 1
for snap in range(flag, frames+flag, flag):
    print(f"Obtaining coordinate information for frame {snap}.")
    ends = snap * (atoms + 2) - 1
    #print(ends)
    begins = ends - atoms -1
    #print(begins)
    with open(f"cmpd_{counter}", "w") as out:
        with open(f"{molecule_name}_MD.arc", "r") as inp:
            for i, line in enumerate(inp): #sed command in bash
                if i >= begins and i <= ends:
                    out.write(line)
                elif i > ends:
                    break
    counter += 1
        #A_glucosecounter += 1
print()

#__________________________________________________________________________________
# ==> Generate Directories and Modify SLURM Submission and Gaussian Input Files <==
#__________________________________________________________________________________


os.mkdir(f"{molecule_name}_MD")
print("Generating directory structure and modifying input and submission files.")
os.chdir(f"{molecule_name}_MD")

conformer_count = 1
while conformer_count <= snapshots:
    #os.chdir(cwd)

    # Makes the directory for a specific conformer.
    os.mkdir(f"cmpd_{conformer_count}")

    # Copies the original submission script into the conformer's directory.
    shutil.copy(f"{cwd}/G09_sub_SLURM.sh", f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/")

    # Makes the submission script executable.
    os.chmod(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/G09_sub_SLURM.sh", 0o755)

    # Copies the original input file into the conformer's directory.
    shutil.copy(f"{cwd}/input.dat", f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/")

    # Modifies the input file with conformer, spectroscopy specific data and inputs the basis set specified by the user.
    with open(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/input.dat", "r+") as file:
        content = file.read()
        content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{spectroscopy}")
        content = content.replace("SPECTROSCOPY", spectroscopy)
        content = content.replace("BASIS_SET", basis_set if basis_set != "Mixed" else "Gen")
        file.seek(0)
        file.write(content)
        file.truncate()

    #os.chdir(cwd)
    conformer_count += 1

print()
os.chdir(cwd)


#_____________________________________________________________________
# ==> Modify Snapshots to Include Only Specified Solvent Molecules <==
#_____________________________________________________________________

print("Generating data for Gaussian input file.")

#___________________________
# ==> Initial Parameters <==
#___________________________

# Set solvent atom types.
O=349
H=350

#_____________________________
# ==> Loop Over Compounds  <==
#_____________________________

# Set the file to be read.
for index in range(1, snapshots+1):
    file=f"cmpd_{index}"
    print(f"Compound {index}")

    #________________________________________________
    # ==> Setup Solute, Solvent, and Total Arrays <==
    #________________________________________________

    # Initializing atom number, atom symbol, X, Y, Z, and atom type arrays for all atoms.
    atom_number=[]
    atom_sym=[]
    X=[]
    Y=[]
    Z=[]
    atom_type=[]

    # Reads the lines in the file and appends to the arrays.
    with open(file, 'r') as f:
        lines = f.readlines()
        #print(os.getcwd())
        for line in lines[2:len(lines)]: #skipping the first two lines
            temp = line.split()  #goes through each field/column in a line
            #print(temp)
            atom_number.append(int(temp[0]))
            atom_sym.append(temp[1])
            X.append(float(temp[2]))
            Y.append(float(temp[3]))
            Z.append(float(temp[4]))
            atom_type.append(int(temp[5]))

    # Initializes solute X, Y, and Z arrays.
    solute_number = []
    solute_sym = []
    solute_X = []
    solute_Y = []
    solute_Z = []

    # Initializes solvent X, Y, and Z arrays.
    solvent_number = []
    solvent_sym = []
    solvent_X = []
    solvent_Y = []
    solvent_Z = []

    # Appends solute and solvent arrays.
    for i in range(len(atom_number)):
        if atom_type[i] != O and atom_type[i] != H:
            solute_number.append(atom_number[i])
            solute_sym.append(atom_sym[i])
            solute_X.append(X[i])
            solute_Y.append(Y[i])
            solute_Z.append(Z[i])
        else:
            solvent_number.append(atom_number[i])
            solvent_sym.append(atom_sym[i])
            solvent_X.append(X[i])
            solvent_Y.append(Y[i])
            solvent_Z.append(Z[i])
    print(f'Number of Solute Atoms: {len(solute_number)}')

    # ___________________________________________
    # ==> Center Solvent Atoms Around Solute <==
    # ___________________________________________

    # Read the size of the box in the X, Y, and Z dimensions.
    box = np.genfromtxt(file, skip_header=1, max_rows=1)
    # with open(f"cmpd_{index}") as f:
    #     box = f.readlines()
    #     sp_box = box.split()
    box_X = float(box[0]) #x-coordinate from the file
    box_Y = float(box[1]) #y-coordinate from the file
    box_Z = float(box[2]) #z-coordinate from the file

    # Compute average coordinate of solute molecule.
    avg_X = np.mean(solute_X)
    avg_Y = np.mean(solute_Y)
    avg_Z = np.mean(solute_Z)

    # Set conditions to center the solvent atoms around the solute.
    for a in range(len(solvent_number)):
        if (solvent_X[a]-avg_X) > box_X/2:
            solvent_X[a] -= box_X
        elif (solvent_X[a]-avg_X) < -box_X/2:
            solvent_X[a] += box_X
        if (solvent_Y[a]-avg_Y) > box_Y/2:
            solvent_Y[a] -= box_Y
        elif (solvent_Y[a]-avg_Y) < -box_Y/2:
            solvent_Y[a] += box_Y
        if (solvent_Z[a]-avg_Z) > box_Z/2:
            solvent_Z[a] -= box_Z
        elif (solvent_Z[a]-avg_Z) < -box_Z/2:
            solvent_Z[a] += box_Z
 
    
    #______________________________________________
    # ==> Setup and Compute Intermediate Arrays <==
    #______________________________________________
    
    # Initialize intermediate solvent arrays for faster computation.
    int_number=[]
    int_sym=[]
    int_X=[]
    int_Y=[]
    int_Z=[]
    
    # Compute maximum distance of atoms in solute from average coordinate.
    max_dist = 0
    for n in range(len(solute_number)):
        dist = math.sqrt((solute_X[n] - avg_X)**2 + (solute_Y[n] - avg_Y)**2 + (solute_Z[n] - avg_Z)**2)
        if dist > max_dist:
            max_dist = dist
    
    # Set the distance threshold for the retaining solvent atoms.
    Dist_thresh = max_dist + float(D)
    
    # Compute intermediate solvent arrays.
    for o in range(len(solvent_number)):
        solvent_dist = math.sqrt((avg_X - solvent_X[o])**2 + (avg_Y - solvent_Y[o])**2 + (avg_Z - solvent_Z[o])**2)
        if solvent_dist <= Dist_thresh:
            int_number.append(solvent_number[o])
            int_sym.append(solvent_sym[o])
            int_X.append(solvent_X[o])
            int_Y.append(solvent_Y[o])
            int_Z.append(solvent_Z[o])
    
    print(f"Number of Intermediate Solvent Atoms: {len(int_number)}")
    
    # ==> Setup and Compute Gaussian Input Arrays <==

    # Initialize arrays for Gaussian input.
    final_number = []
    final_sym = []
    final_X = []
    final_Y = []
    final_Z = []
    
    # Appending arrays for Gaussian input with solute data.
    final_number.extend(solute_number)
    final_sym.extend(solute_sym)
    final_X.extend(solute_X)
    final_Y.extend(solute_Y)
    final_Z.extend(solute_Z)
    
    # Determines solvent molecules distance from solute atoms and appends.
    for j in range(len(solute_number)):
        #if solute_sym[j] != 'H':
        for k in range(len(int_number)):
            d = math.sqrt((solute_X[j] - int_X[k])**2 + (solute_Y[j] - int_Y[k])**2 + (solute_Z[j] - int_Z[k])**2)
            if d <= float(D):
                no_duplicates = True
                for l in range(len(final_number)):
                    if int_number[k] == final_number[l]:
                        no_duplicates = False
                        break
                if no_duplicates:
                    final_number.append(int_number[k])
                    final_sym.append(int_sym[k])
                    final_X.append(int_X[k])
                    final_Y.append(int_Y[k])
                    final_Z.append(int_Z[k])
    
    print(f"Number of Atoms before Solvent Check: {len(final_number)}")
    
    #_________________________________________________
    # ==> Confirm Presence of Full Solvent Molecules <==
    #_________________________________________________
    
    
    solvent_complete = False
    poorni = 0

    while not solvent_complete:
        #print(os.getcwd())
        # Initializing solvent connectivity arrays.
        connect_1 = []
        connect_2 = []
        connect_3 = []
        connect_4 = []
        col = []
        columns = 0
        
        # with open (f"cmpd_{index}") as f:
        #     content = f.readlines()
        #     for line in content:
        #         splitted_line = line.split()
        #         for p in range(len(final_number)):
        #             ln_num = final_number[p]
        #             if splitted_line[0] == str(ln_num) and splitted_line[1] != "Great":
        #                 columns = len(splitted_line)
        #                 break
        #         if columns == 6:
        #             connect_1.append(0)
        #             connect_2.append(0)
        #             connect_3.append(0)
        #             connect_4.append(0)
        #             col.append(columns)
        #         elif columns == 7:
        #             connect_1.append(int(splitted_line[6]))
        #             connect_2.append(0)
        #             connect_3.append(0)
        #             connect_4.append(0)
        #             col.append(columns)
        #         elif columns == 8:
        #             connect_1.append(int(splitted_line[6]))
        #             connect_2.append(int(splitted_line[7]))
        #             connect_3.append(0)
        #             connect_4.append(0)
        #             col.append(columns)
        #         elif columns == 9:
        #             connect_1.append(int(splitted_line[6]))
        #             connect_2.append(int(splitted_line[7]))
        #             connect_3.append(int(splitted_line[8]))
        #             connect_4.append(0)
        #             col.append(columns)
        #         elif columns == 10:
        #             connect_1.append(int(splitted_line[6]))
        #             connect_2.append(int(splitted_line[7]))
        #             connect_3.append(int(splitted_line[8]))
        #             connect_4.append(int(splitted_line[9]))
        #             col.append(columns)
        #         elif columns!=0:
        #             col.append(columns)
        
        
        
        lb= [True, True, True, True]
        l = [connect_1, connect_2, connect_3, connect_4]
    
        for p in range(len(final_number)):
            ln_num = final_number[p]
            with open(f"cmpd_{index}") as f:
                #print("hi")
                content = f.readlines()
                for line in content:
                    splitted_line  = line.split()
                    if splitted_line[0] == str(ln_num) and splitted_line[1] != "Great": #string
                        columns = len(splitted_line)
                        break
                
                l[0].append(int(splitted_line[6]))
                l[1].append(0 if columns<8 else int(splitted_line[7]))
                l[2].append(0 if columns<9 else int(splitted_line[8]))
                l[3].append(0 if columns<10 else int(splitted_line[9]))
                col.append(columns)
                
                # for i in range(min(4, columns - 6)):
                #     for j in range(4):
                #         l[j][i] = int(splitted_line[i+6+j] if i+6+j < columns else 0)
                        

                    # if columns == 6:
                    #     connect_1.append(0)
                    #     connect_2.append(0)
                    #     connect_3.append(0)
                    #     connect_4.append(0)
                    # elif columns == 7:
                    #     connect_1.append(int(splitted_line[6]))
                    #     connect_2.append(0)
                    #     connect_3.append(0)
                    #     connect_4.append(0)
                    # elif columns == 8:
                    #     connect_1.append(int(splitted_line[6]))
                    #     connect_2.append(int(splitted_line[7]))
                    #     connect_3.append(0)
                    #     connect_4.append(0)
                    # elif columns == 9:
                    #     connect_1.append(int(splitted_line[6]))
                    #     connect_2.append(int(splitted_line[7]))
                    #     connect_3.append(int(splitted_line[8]))
                    #     connect_4.append(0)
                    # elif columns == 10:
                    #     connect_1.append(int(splitted_line[6]))
                    #     connect_2.append(int(splitted_line[7]))
                    #     connect_3.append(int(splitted_line[8]))
                    #     connect_4.append(int(splitted_line[9]))
                    # col.append(columns)

    
    
        # Initialize new atom array.
    
        new_number = []
        for q in range(len(final_number)):
            for r in range(len(l[0])):
                for y in range(4):
                    if l[y][q] == final_number[r]:
                        lb[y] = True
                for y in range(4):
                    if l[y][q] == 0:
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
                            if l[y][q] == new_number[t]:
                                no_duplicates = False
                                break
                        if no_duplicates:
                            new_number.append(l[y][q])
        
        # Checking status of atoms connectivity.
        # for q in range(len(final_number)):
        
        #     # Initializing status of connectivity.
        #     connect_1_present = False
        #     connect_2_present = False
        #     connect_3_present = False
        #     connect_4_present = False
        #     for r in range(len(connect_1)):
        #         if connect_1[q] == final_number[r]:
        #             connect_1_present = True
        #         elif connect_2[q] == final_number[r]:
        #             connect_2_present = True
        #         elif connect_3[q] == final_number[r]:
        #             connect_3_present = True
        #         elif connect_4[q] == final_number[r]:
        #             connect_4_present = True
        #         if connect_1[q] == 0:
        #             connect_1_present = True
        #             connect_2_present = True
        #             connect_3_present = True
        #             connect_4_present = True
        #         elif connect_2[q] == 0:
        #             connect_2_present = True
        #             connect_3_present = True
        #             connect_4_present = True
        #         elif connect_3[q] == 0:
        #             connect_3_present = True
        #             connect_4_present = True
        #         elif connect_4[q] == 0:
        #             connect_4_present = True
        
        #     # Checking for connectivity and duplicates.
        #         if not connect_1_present:
        #             no_duplicates = True
        #             for t in range(len(new_number)):
        #                 if connect_1[q] == new_number[t]:
        #                     no_duplicates = False
        #                     break
        #             if no_duplicates:
        #                 new_number.append(connect_1[q])
        #         elif not connect_2_present:
        #             no_duplicates = True
        #             for t in range(len(new_number)):
        #                 if connect_2[q] == new_number[t]:
        #                     no_duplicates = False
        #                     break
        #             if no_duplicates:
        #                 new_number.append(connect_2[q])
        #         elif not connect_3_present:
        #             no_duplicates = True
        #             for t in range(len(new_number)):
        #                 if connect_3[q] == new_number[t]:
        #                     no_duplicates = False
        #                     break
        #             if no_duplicates:
        #                 new_number.append(connect_3[q])
        #         elif not connect_4_present:
        #             no_duplicates = True
        #             for t in range(len(new_number)):
        #                 if connect_4[q] == new_number[t]:
        #                     no_duplicates = False
        #                     break
        #             if no_duplicates:
        #                 new_number.append(connect_4[q])
            
            # Checking if solvent is complete.
            print("Number of New Atoms:", len(new_number))
            if len(new_number) == 0:
                print("Solvent Complete.")
                solvent_complete = True
                break
            else:
                print("Solvent Incomplete. Appending and obtaining connectivity of new atoms.")
    
            # Appending new atoms to final atom arrays.
            for s in range(len(new_number)):
                ind = new_number[s] - len(solute_number)
                final_number.append(solvent_number[ind])
                final_sym.append(solvent_sym[ind])
                final_X.append(solvent_X[ind])
                final_Y.append(solvent_Y[ind])
                final_Z.append(solvent_Z[ind])
            poorni += 1
            print(poorni)
        
    print("Final Number of Atoms:", len(final_number))
    print(" ")
    
    #______________________________________________
    # ==> Writing Gaussian Input Arrays to File <==
    #______________________________________________
    
        # import os
        # print("hi")
        # Overwriting compound files with final data.
    for i in range(len(final_number)):
        with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
            f.write("{:<2} \t {:<10} \t {:<10} \t {:<10}\n".format(final_sym[i], final_X[i], final_Y[i], final_Z[i]))
    
        with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
            f.write(" \n")
    
        if basis_set == "Mixed":
            with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
                f.write("1 - {} 0\n".format(len(solute_number)))
                f.write("{}\n".format(solute_basis_set))
                f.write("****\n")
                f.write("{} - {} 0\n".format(len(solute_number)+1, len(final_number)))
                f.write("{}\n".format(solvent_basis_set))
                f.write("****\n")
    
        with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
            f.write(" \n")
            f.write(" \n")
        
os.chdir(cwd)
for index in range(1,snapshots+1):
    os.remove(f"cmpd_{index}")

        
                
                
            
