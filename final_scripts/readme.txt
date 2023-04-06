The set of scripts required to a simulate a VCD/ROA spectrum from a Molecular Dynamics trajectory file. 
The snapshot_hdf5.py creates separate input files for the Gaussain program to calculate the parameters for the soecific spectrum. It does so by selecting frames in a semi-random manner (to ensure linear independence between the frames) and creates input files for each frame. 
The program makes use of hdf5 datastructure which stores the relevant columns as arrays and fasten up the code. 
Once the input files are ready, they are submitted to run a Gaussian calculation. This is doen using the submit.py program. 
data_MD.py filters out the required data from each of the Gaussian output files and creates csv files for plotting. 
The corresponding spectrum is then plotted by using the final_converged jupyter notebook to check for convergence with the number of snapshots. The code is optimized by using vectorized operations which gives a considerabl speedup. The final spectrum is animated with the number of snapshots, with overlap(error) plotted in a subgraph. 
