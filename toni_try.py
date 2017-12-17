import os
import numpy as np

def output_distances(input,output):
  f=open(input,"r")
  coord=[]
  for line in f:
      if line.startswith("ATOM") and line[12:16].strip()=='CB':
          x,y,z=float(line[30:38]),float(line[38:46]),float(line[46:54])
          coord.append((line[22:26].strip(),x,y,z))

  dist1=[]
  for i, elem in enumerate(coord):
      tmp=[]
      for j in coord:
	  dist=(coord[i][0],np.sqrt(((elem[1]-j[1])**2+(elem[2]-j[2])**2+(elem[3]-j[3])**2)))
	  tmp.append(dist)
      dist1.append(tmp)

  out=open(output,"w")
  out.write("ID\t")
  for i, elem in enumerate(coord):
      out.write(elem[0]+'\t')
  out.write('\n')
  for i, elem in enumerate(dist1):
      out.write(coord[i][0]+'\t')
      for d in elem:
          out.write(str(d[1])+'\t')
      out.write('\n')
  out.close()
  f.close()

return dist1

input1=""
output1=""
input2=""
output2=""
general_out=""

dist1=output_distances(input1,output1)
dist2=output_distances(input2,output2)
out=open(general_out,"w")
out.write("ID\t")
for i, elem in enumerate(dist1):
    out.write(elem[i][0]+'\t')
out.write('\n')
for i, elem in enumerate(dist1):
    out.write(dist1[i][0]+'\t')
    for d in elem:
        out.write(str(d[1])+'\t')
    out.write('\n')
out.close()
f.close()
