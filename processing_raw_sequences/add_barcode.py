#!/usr/bin/python
#adding barcode sequences to forward and reverse read sequences
#author: Benli Chai
 
import sys, string, os

def addTag(f1, f2, f3):
	R1_OUT = open(R1_out, 'w')
	R2_OUT = open(R2_out, 'w')
	while 1:
		try:
			R1 = []#lines in R1
			R2 = []#lines in R2
			R3 = []#lines in R3
			for i in [1, 2, 3, 4]:#collect four lines of each record
				R1.append(f1.next().strip())
				R2.append(f2.next().strip())
				R3.append(f3.next().strip())
			assert R1[0].split()[0] == R2[0].split()[0] == R3[0].split()[0]
			tag = R3[1].strip()
			R1[0] = R1[0].strip('0') + tag

			R2[0] = R2[0].strip('0') + tag
		 
			R1_OUT.write(string.join(R1, '\n') + '\n')	
			R2_OUT.write(string.join(R2, '\n') + '\n')	

		except StopIteration:
			break
	R1_OUT.close()
	R2_OUT.close()
		
#Main#
#sequence input handles
if not len(sys.argv) == 4:
	print 'addTag.py forwordRead.fastq reverseRead.fastq indexRead.fastq'
	sys.exit()

R1_in = open(sys.argv[1], "r")#forward read
R2_in = open(sys.argv[2], "r")#reverse read
R3_in = open(sys.argv[3], "r")#index read

R1_out = os.path.basename(sys.argv[1]).split('.fastq')[0] + '_tagged.fastq'
R2_out = os.path.basename(sys.argv[2]).split('.fastq')[0] + '_tagged.fastq'


addTag(R1_in, R2_in, R3_in)
