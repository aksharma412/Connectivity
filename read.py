import numpy as np
import h5py


no_atoms = []
indx =[]
coord = []
at_type = []
connect = []
dims_box = np.zeros((3))
frame_no = 1

with h5py.File("dumb.h5" , "w") as h5:
    with open("test1" , "r") as f:
        content = f.readlines()
        content.pop(0)
        for i , line in enumerate(content):
            splt = line.split()
            if splt[1] != "Great":
                if splt[-1] == "90.000000":
                    for a in range(3):
                        dims_box[a] = float(splt[a])
                elif len(splt)>6:
                    connect_line = np.zeros((4))
                    xyz_line = np.zeros((3))
                    indx.append(int(splt[0]))
                    at_type.append(int(splt[5]))
                    
                    for a in range(3):
                        xyz_line[a]  = float(splt[a+2])
                    coord.append(xyz_line)

                    for a in range(len(splt)-6):
                        connect_line[a] = int(splt[6+a])
                    connect.append(connect_line)
            
            elif splt[1]=="Great" or len(splt) ==0:
                grp = h5.create_group(f"frame{frame_no}")
                grp.create_dataset("no_atoms" , data = np.array(no_atoms))
                grp.create_dataset("dims_box" , data = np.array(dims_box))
                grp.create_dataset("indx" ,     data = np.array(indx))
                grp.create_dataset("coord" ,    data = np.array(coord))
                grp.create_dataset("at_type" ,  data = np.array(at_type))
                grp.create_dataset("connect" ,  data = np.array(connect))
                frame_no +=1
                no_atoms = []
                indx =[]
                coord = []
                at_type = []
                connect = []
                dims_box = np.zeros((3))

