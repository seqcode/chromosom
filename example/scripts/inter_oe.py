from multimds import data_tools as dt
import sys
import numpy as np
from multimds.hic_oe import oe

def interMatFromBed(path, struct1, struct2):
	mat = np.zeros((len(struct1.getPoints()), len(struct2.getPoints())))
	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			index1 = struct1.get_rel_index(int(line[4]))
			index2 = struct2.get_rel_index(int(line[1]))
			if index1 is not None and index2 is not None:
				mat[index1, index2] += float(line[6])
		infile.close()
	return mat

def threshold(mat, value):
	"""Cuts off values above threshold"""
	n = len(mat)
	thresholded = np.zeros_like(mat)
	for i in range(n):
		for j in range(len(mat[0])):
			thresholded[i,j] = min((mat[i,j], value))
	return thresholded

prefix = sys.argv[1]	#e.g. GM12878
resol = sys.argv[2] #e.g. 100kb
hicdir = sys.argv[3] #e.g. hic_K562_100000

chroms = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
#if prefix == 'K562':
#	chroms = ["1","2","3","4","5","6","7","8","10","11","12","13","14","15","16","17","18","19","20","21","X"] #ignore translocated

structs = [dt.structureFromBed("{}/{}_{}_{}.bed".format(hicdir, prefix, chrom, resol)) for chrom in chroms]	#structures hold non-empty loci ("points") and metadata

num_bins = sum([len(struct.getPoints()) for struct in structs])		#total number of bins across genome

full_mat_contact = np.zeros((num_bins, num_bins))
full_mat_oe = np.zeros((num_bins, num_bins))

#Print coordinates of bins
coordfilename = "{}/{}_{}.coords".format(hicdir, prefix, resol)
cfile = open(coordfilename, "w")
for i in range(len(chroms)):
	coords = structs[i].getGenCoords()
	for x in range(len(coords)): 
		cfile.write(str(chroms[i])+'\t'+str(coords[x])+'\n')

cfile.close()


#Generate O/E matrices
row_offset = 0
for i in range(len(chroms)):
	col_offset = 0		#reset column offset at beginning of each row
	num_bins2 = len(structs[i].getPoints())
	for j in range(i):
		num_bins1 = len(structs[j].getPoints())
		mat = interMatFromBed("{}/{}_{}_{}_{}.bed".format(hicdir, prefix, chroms[j], chroms[i], resol), structs[i], structs[j])
		mat_oe = mat/np.mean(mat)	#observed/expected
		mat_oe = threshold(mat_oe, 10)	#remove outliers
		full_mat_contact[row_offset:row_offset+num_bins2, col_offset:col_offset+num_bins1] = mat	#put data in full matrix
		full_mat_oe[row_offset:row_offset+num_bins2, col_offset:col_offset+num_bins1] = mat_oe	#put data in full matrix
		col_offset += num_bins1		#update offset
	mat = dt.matFromBed("{}/{}_{}_{}.bed".format(hicdir, prefix, chroms[i], resol), structure=structs[i]) 	#intrachromosomal
	mat_oe = oe(mat, structs[i])
	full_mat_contact[row_offset:row_offset+num_bins2, col_offset:col_offset+num_bins2] = mat	#put data in full matrix
	full_mat_oe[row_offset:row_offset+num_bins2, col_offset:col_offset+num_bins2] = mat_oe	#put data in full matrix
	row_offset += num_bins2		#update offset

#make symmetric
for i in range(num_bins):
	for j in range(i):
		full_mat_contact[j,i] = full_mat_contact[i,j]
		full_mat_oe[j,i] = full_mat_oe[i,j]

np.savetxt("{}/{}_{}_interchromosomal.contact.tsv".format(hicdir, prefix, resol), full_mat_contact, delimiter="\t")
np.savetxt("{}/{}_{}_interchromosomal.oe.tsv".format(hicdir, prefix, resol), full_mat_oe, delimiter="\t")

