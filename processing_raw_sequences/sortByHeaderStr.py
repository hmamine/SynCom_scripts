#!/usr/bin/python
#author: Benli Chai
#sort reads by a unique substring in the header - requires Biopython 1.51 or later
from Bio import SeqIO
import itertools, sys

#Setup variables (could parse command line args instead)

format = 'fasta'
def sortByTag(iter):
	handles = {}#list of opened output file handles
	handle = open('notag.fa', 'w')
	handles['notag'] = handle
	for read in iter:
		find = 0 #switch notag and tagged reads
		for tag in tagList:
			if read.id.split(';')[0].split('_')[0] == tag:
				find = 1
				if not tag in handles.keys():#new sample new output handler
					name = tag  + '.fa'
					handle = open(name, 'w')
					handles[tag] = handle

					SeqIO.write(read, handle, format)
				else:
					handle = handles[tag]
					SeqIO.write(read, handle, format)
				break#find the barcode then break out

		if find == 0:#no tag was found for this pair of reads
			handle = handles['notag']
			SeqIO.write(read, handle, format)
	
	for handle in handles.values():#close all the output file handles
		handle.close()

if len(sys.argv) != 3:
	print 'sortByTag.py <ReadFile.fasta> <tagfile> '
	sys.exit()
#Main#
#sequence input handles
ReadFile = SeqIO.parse(open(sys.argv[1], "rU"), format)

#get sample and tags
tagFile = open(sys.argv[2], 'r').readlines()

hash = {}
tagList = []
for line in tagFile:
	tag = line.strip()
	if not tag in tagList:
		tagList.append(tag)
	else:
		print 'Duplicated tag', line.strip()
		sys.exit()
sortByTag(ReadFile) 
