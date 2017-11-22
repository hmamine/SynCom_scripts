#!/usr/bin/python 
#author: Benli Chai
#generate count table from alignment files
import sys, string, os
if len(sys.argv) < 5:
	print 'Usage: getSum.py sampleIDlist otuIDlist cutoff pairwiseKNNfiles'
	sys.exit()


fS = open(sys.argv[1], 'r').readlines()#sample ID file
fR = open(sys.argv[2], 'r').readlines()#OTU ID file
IDENT = float(sys.argv[3])
names = sys.argv[4:]#pairwise-knn files
samples = []
otus = []
for line in fS:
	samples.append(line.strip())

for line in fR:
	otus.append(line.strip())

hash = {}
total = 0
passCount = 0
for input in names:
	f = open(input, 'r')
	while 1:
		try:
			line = f.next()
			if line.strip()[0] == '@':
				total += 1
				cols = line.strip().split('\t')
				rID = cols[0]
				sampleID = rID.split('_')[0].strip('@')
				otuID = cols[10].split()[0]
				ident = float(cols[4])
				if ident >= IDENT:
					passCount += 1
					if not sampleID in hash.keys():
						hash[sampleID] = {}
						hash[sampleID][otuID] = 1
					elif not otuID in hash[sampleID].keys():
						hash[sampleID][otuID] = 1
					else:
						hash[sampleID][otuID] += 1
		except StopIteration:
                	break
print '\t' + string.join(otus, '\t')

for sampleID in samples:
	out = [sampleID]
	for otuID in otus:
		try:
			count = hash[sampleID][otuID]
		except KeyError:
			count = 0
		out.append('%s'%count)
	print string.join(out, '\t')
print 'Total:', total
print '%s sequences >= %s identity to reference'%(passCount, IDENT)
		

