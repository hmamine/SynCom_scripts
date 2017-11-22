#!/usr/bin/python
#rename read sequence IDs by sample IDs using the barcode mapping file 
#author: Benli Chai

import os, sys, string

def NameBySample(f):
	OUT = open(ReadOut, 'w')
	notag = open('NotMapped.fastq', 'w')

	while 1:
		try:
			Read = []#lines in one read
			if format == 'fasta':
				for i in [1, 2]:
					Read.append(f.next().strip())
			elif format == 'fastq':
				for i in [1, 2, 3, 4]:#collect four lines of each record
					Read.append(f.next().strip())
			head = Read[0]
			barcode = head.split(':')[-1]
			try:
				sample = Map[barcode] 
				if not sample in sampleRead.keys():
					sampleRead[sample] = 1
				else:
					sampleRead[sample] += 1
				count = sampleRead[sample] # the count of read from this sample so far
				if format == 'fasta':
					head = '>' + sample + '_Read%s'%count + ' ' + head.strip('>')
				elif format == 'fastq':
					head = '@' + sample + '_Read%s'%count + ' ' + head.strip('@')
				Read[0] = head #replace the header line with modified head with sample name and count as the I 
				OUT.write(string.join(Read, '\n') + '\n')	
			except KeyError:
				notag.write(string.join(Read, '\n') + '\n')	
				 

		except StopIteration:
			break
	OUT.close()
	notag.close()
		
#Main#
#sequence input handles
if not len(sys.argv) == 4:
	print 'NameBySample.py pandaseqAssembled.fastq barcodeMap format[fasta|fastq]'
	sys.exit()

ReadIn = open(sys.argv[1], "r")#pandaseq assembled reads
mapFile = open(sys.argv[2], "r").readlines()#barcode-sample map
format = sys.argv[3]

Map = {}
sampleRead = {}#record the count of reads from ech sample
for line in mapFile:
	line = line.strip()
	if line[0] == '#':
		continue
	sample, barcode = line.split('\t')[:2]
	Map[barcode] = sample	


print Map
ReadOut = os.path.basename(sys.argv[1]).split('.%s'%format)[0] + 'SampleName.%s'%format


NameBySample(ReadIn)
