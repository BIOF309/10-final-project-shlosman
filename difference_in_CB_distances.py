# Name: difference_in_CB_distances.py
# Author: Antoniya Aleksandrova & Irina Shlosman
# Date: 4 Dec 2017
# Description: Given the two PDB files of the same protein, it calculates all distances between their CB atoms, outputs them and then outputs a file with the difference between the distances in the two files (diff_CB_dist.txt)

import os
import numpy as np

def output_distances(input,output): # input is the path to the in file, output is the path to the out file
    f=open(input,"r")
    coord=[] # stores the residue index and coordinates of the CB atom
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip()=='CB':
            x,y,z=float(line[30:38]),float(line[38:46]),float(line[46:54])
            coord.append((line[22:26].strip(),x,y,z))

    f.close()
    dist1=[]  # list of the form [(2,[0.0,5.43,2.41,...]),(3,[1.7,0.0,14.23,5.6,...),...]; i.e. each element is (#residue_id, [distance to CB, distance to CB, ...])
    for i, elem in enumerate(coord):
        tmp=[]
        for j in coord:
            dist=(np.sqrt(((elem[1]-j[1])**2+(elem[2]-j[2])**2+(elem[3]-j[3])**2)))  # stores the distance between the current residue's CB atom and the CB atoms of all other residues in the protein
            tmp.append(dist)
        dist1.append((elem[0],tmp))

    out=open(output,"w")
    out.write("ID\t")
    for i, elem in enumerate(dist1):
        out.write(elem[0]+'\t')
    out.write('\n')
    for i, elem in enumerate(dist1):
        out.write(elem[0]+'\t')
        for d in elem[1]:
            out.write(str(d)+'\t')
        out.write('\n')
    out.close()

    return dist1


# Main
input1="Aligned_Model.pdb"
output1=input1[:-4]+"_CB_dist.txt"
input2="clean_struct.pdb"
output2=input2[:-4]+"_CB_dist2.txt"
general_out="diff_CB_dist.txt"

dist1=output_distances(input1,output1)
dist2=output_distances(input2,output2)
out=open(general_out,"w")
out.write("ID\t")
for i, elem in enumerate(dist1):
    out.write(elem[0]+'\t')
out.write('\n')
for i, elem in enumerate(dist1):
    out.write(elem[0]+'\t')
    for j, d1 in enumerate(elem[1]):
        assert elem[0]==dist2[i][0], "Not the same protein! Check residue ids!"
        out.write(str(float('%.3f'% abs(d1-dist2[i][1][j])))+'\t')
    out.write('\n')
out.close()
