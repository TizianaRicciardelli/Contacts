#script to retrieve interchain from a pdb res-res interaction between heavy atoms
#usage python contacts.py XXXX AB CDE 5
import os
import sys
from math import sqrt

#read pdb and append only ATOM and HETATM lines
pdb=[]
with open(sys.argv[1]+".pdb",'r') as input_file:
  for line in input_file:
    if line.startswith("ATOM"): # or line.startswith("HETATM"):   #activate the second part if you want to retrieve also contacts with non-heavy atoms
      pdb.append(line)
#create list of input chains ex.: AB vs CDE
chain1=[]
chain2=[]
for elem1 in sys.argv[2]:
	chain1.append(elem1)

for elem2 in sys.argv[3]:
        chain2.append(elem2)

#dist=int(sys.argv[4])

dicts_chain1={}
dicts_chain2={}
for p in range(0,len(pdb)):
	if pdb[p][21:22] in chain1: #21:22 is chain ID position in pdb format
		for m in range(p+1,len(pdb)):
			if pdb[m][21:22] in chain2:
				x1=float(pdb[p][30:38])
				y1=float(pdb[p][38:46])
				z1=float(pdb[p][46:54])
				x2=float(pdb[m][30:38])
				y2=float(pdb[m][38:46])
				z2=float(pdb[m][46:54])
				x=(x1-x2)**2
				y=(y1-y2)**2
				z=(z1-z2)**2
				my_dist= sqrt(x+y+z) #distance between 2 points in 3D
        if my_dist <=dist: #activate this line when you want a different distance threshold
#				if my_dist <=5:
					chain_ID1=str(pdb[p][21:22])
					chain_ID2=str(pdb[m][21:22])
					resi_name1=str(pdb[p][17:20])
					resi_name2=str(pdb[m][17:20])
					resi_num1=str(pdb[p][22:27])
					resi_num2=str(pdb[m][22:27])
					atom_name1=str(pdb[p][13:16])
					atom_name2=str(pdb[m][13:16])
					key=chain_ID1+"_"+resi_name1+"_"+str(resi_num1)+"-"+chain_ID2+"_"+resi_name2+"_"+str(resi_num2)
					if key not in dicts_chain1:
						dicts_chain1[key]=[]
						dicts_chain1[key].append(my_dist)
						dicts_chain1[key].append(atom_name1)
						dicts_chain1[key].append(atom_name2)
					else:
						if my_dist < dicts_chain1[key][0]:
							dicts_chain1[key]=[my_dist,atom_name1,atom_name2]
						#elif my_dist == dicts_chain1[key][0]: for now, let's skeep atom differenxes in distances
	elif pdb[p][21:22] in chain2:
		for m in range(p+1,len(pdb)):
			if pdb[m][21:22] in chain1:
				x1=float(pdb[p][30:38])
				y1=float(pdb[p][38:46])
				z1=float(pdb[p][46:54])
				x2=float(pdb[m][30:38])
				y2=float(pdb[m][38:46])
				z2=float(pdb[m][46:54])
				x=(x1-x2)**2
				y=(y1-y2)**2
				z=(z1-z2)**2
				my_dist= sqrt(x+y+z) #distance between 2 points in 3D
				if my_dist <=dist: #activate this line when you want a different distance threshold
#				if my_dist <=5:
					chain_ID2=str(pdb[p][21:22])
					chain_ID1=str(pdb[m][21:22])
					resi_name2=str(pdb[p][17:20])
					resi_name1=str(pdb[m][17:20])
					resi_num2=str(pdb[p][22:27])
					resi_num1=str(pdb[m][22:27])
					atom_name2=str(pdb[p][13:16])
					atom_name1=str(pdb[m][13:16])
					key=chain_ID1+"_"+resi_name1+"_"+str(resi_num1)+"-"+chain_ID2+"_"+resi_name2+"_"+str(resi_num2)
					if key not in dicts_chain1:
						dicts_chain1[key]=[]
						dicts_chain1[key].append(my_dist)
						dicts_chain1[key].append(atom_name1)
						dicts_chain1[key].append(atom_name2)
					else:
						if my_dist < dicts_chain1[key][0]:
							dicts_chain1[key]=[my_dist,atom_name1,atom_name2]
			                        #elif my_dist == dicts_chain1[key][0]: for now, let's skeep atom differenxes in distances
			



#print(dicts_chain1)
#print(len(dicts_chain1))


for key, value in dicts_chain1.items():
    key_parts = key.split('-')
    chain_1 = key_parts[0][0]
    chain_2 = key_parts[1][0]
    residue_1 = key_parts[0][2:]
    residue_2 = key_parts[1][2:]
    atom = value[1].strip()
    distance = value[0]

    output = f"{residue_1},{atom},{chain_1},{residue_2},{value[2]},{chain_2},{distance}"
    print(output)
