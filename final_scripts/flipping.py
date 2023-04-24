import numpy as np
import os

molecule_name = input("Enter the molecule name :")
snapshot = input("Enter the number of snapshots : ")

cwd = os.getcwd()
for i in range(1,len(snapshots)+1):
	at_name = []
	at_x = []
	at_y = []
	at_z = []
	os.chdir(f"{molecule_name}_MD"/f"cmpd_{i}")
	with open("input.dat", "r") as f:
		content = f.readlines()
	temp = content
	content = content[8:122]
	for line in content:
		splitted_line = line.split()
		at_name.append(splitted_line[0])
		at_x.append(float(splitted_line[1])*-1)
		at_y.append(float(splitted_line[2])*-1)
		at_z.append(float(splitted_line[3])*-1)
	os.remove("input.dat")
	with open("input.dat","a") as f:
		f.writelines(temp[:8])
		for i in range(len(at_x)):
			f.write("{:<2} \t {:<10} \t {:<10} \t {:<10}\n".format(at_name[i], at_x[i], at_y[i], at_z[i]))
		f.writelines(temp[122:])

	

