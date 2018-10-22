#!/usr/bin/python

###########################################################################################################

import sys

###########################################################################################################
#
#												Description
#
# concat_RM.py concatenates similar successive hits from RM output tab file 
#
# It considers only hits closer than 500 bp (can be modified in the concatHSP function)
#
# It returns only the final hits larger than 300 bp (can be modified in the length_filter function)
#
#												Requirements
#
# The RM output tab file needs to be sorted on the chromosome id and then on the start position of the hit
#
# Sequences from the database used for the RM search should have an id format similar to 'Gypsy_xx'
#
#												Usage
#
# concat_RM.py RM_output_tab_file > output_concat_file
#
###########################################################################################################

def parse_args():      
	file = sys.argv[1:]
	return file

def readList(f):
	hsp = {}
	infile=open(f, 'r')
	lines = infile.readlines()
	for l in lines :
		li=l.split()
		lii=str(li[10]).split('_')

			
		if (hsp.has_key(li[4])):
			hsp[li[4]].append([str(li[5]), str(li[6]), str(lii[0])])
		else :
			hsp[li[4]] = [[str(li[5]), str(li[6]), str(lii[0])]]

	return hsp


def concatHSP (dico):
	hsp2 = {}
	for i in dico.keys():
		x=1
		v=0
		if len(dico[i]) > 1 :
			while x<len(dico[i]) :
				if (hsp2.has_key(i)):
					tag1 = str(hsp2[i][v][2])
					tag2 = str(dico[i][x][2])
					start1 = int(hsp2[i][v][0])
					end1 = int(hsp2[i][v][1])
					start2 = int(dico[i][x][0])
					end2 = int(dico[i][x][1])
					if (tag1==tag2 and start2<=end1+500) :
						if end2>=end1 :
							new_end = int(end2)
							hsp2[i].pop(v)
							if (hsp2.has_key(i)):
								hsp2[i].append([int(start1), int(new_end), str(tag1)])
							else :
								hsp2[i] = [[int(start1), int(new_end), str(tag1)]]
						else :
							nothing=0
					else :
						hsp2[i].append([int(start2), int(end2), str(tag2)])
						v=v+1
						if start2<=end1 :
							print >> sys.stderr, 'Alert: two overlaping hits belonging to two different superfamilies on', str(i), 'at', str(start1), 'and', str(start2)
						else :
							nothing=0
				else :
					tag1 = str(dico[i][x-1][2])
					tag2 = str(dico[i][x][2])
					start1 = int(dico[i][x-1][0])
					end1 = int(dico[i][x-1][1])
					start2 = int(dico[i][x][0])
					end2 = int(dico[i][x][1])
					if (tag1==tag2 and start2<=end1+500) :
						if end2>=end1 :
							new_end = int(end2)
						else :
							new_end = int(end1)
						hsp2[i] = [[int(start1), int(new_end), str(tag1)]]
						v=0							
					else :
						hsp2[i] = [[int(start1), int(end1), str(tag1)]]
						hsp2[i].append([int(start2), int(end2), str(tag2)])
						if start2<=end1 :
							print >> sys.stderr, 'Alert: two overlaping hits belonging to two different superfamilies on', str(i), 'at', str(start1), 'and', str(start2)
						else :
							nothing=0
						v=1
				x=x+1
		else :
			tag1 = str(dico[i][x-1][2])
			start1 = int(dico[i][x-1][0])
			end1 = int(dico[i][x-1][1])
			hsp2[i] = [[int(start1), int(end1), str(tag1)]]

	return hsp2
	
def length_filter (dico3) :
	hsp3 = {}
	for i in dico3.keys():
		x=0
		if len(dico3[i]) > 1 :
			while x<len(dico3[i]) :
				tag2 = str(dico3[i][x][2])
				start2 = int(dico3[i][x][0])
				end2 = int(dico3[i][x][1])
				length = end2-start2
				if length > 300 :
					if (hsp3.has_key(i)):
						hsp3[i].append([int(start2), int(end2), str(tag2)])
					else :
							hsp3[i] = [[int(start2), int(end2), str(tag2)]]
				else :
					nothing=0
				x=x+1
		else :
			tag2 = str(dico3[i][x][2])
			start2 = int(dico3[i][x][0])
			end2 = int(dico3[i][x][1])
			length = end2-start2
			if length > 300 :
				hsp3[i] = [[int(start2), int(end2), str(tag2)]]
			else :
				nothing=0

	return hsp3

					
def printout (dico2) :
	for i in dico2.keys():
		x=0
		if len(dico2[i]) > 1 :
			while x<len(dico2[i]) :
				start = int(dico2[i][x][0])
				end = int(dico2[i][x][1])
				length = end-start
				Len = str(length)+'bp'
				print str(i), str(dico2[i][x][0]), str(dico2[i][x][1]), str(dico2[i][x][2]), str(Len)
				x=x+1
		else:
			start = int(dico2[i][x][0])
			end = int(dico2[i][x][1])
			length = end-start
			Len = str(length)+'bp'
			print str(i), str(dico2[i][x][0]), str(dico2[i][x][1]), str(dico2[i][x][2]), str(Len)

	return


####################################################
# Main
####################################################
  
def main():
	args = parse_args()
	bl = readList(args[0])
	d1 = concatHSP(bl)
	d2 = length_filter(d1)
	d3 = concatHSP(d2)
	printout (d3)
  
main()			

