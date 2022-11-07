#!/usr/bin/env python

import sys
#Opens an infile specified by the user. Should be a list of the sequence names you want to extract, each name on a new line
HITS = open(sys.argv[1], 'r')

#Opens a fasta file with the fasta sequences to extract from
DBSEQS = open(sys.argv[2], 'r')

#Opens an output text file as specified by user
OUT = open(sys.argv[3], 'w')

status=0
matches={}

linenum=0
for match in HITS:
	linenum+=1
	match=match.rstrip()
	cols=match.split('\t')
	if linenum>=1:
		matches[cols[0]]=''
	

for line in DBSEQS:
	line=line.rstrip()
	if line.startswith('>'):

#	For returning sequence IN the HitsList
		try:
			matches[line[1:]]
			OUT.write(str(line) + '\n')
			status = 1
		except KeyError:
			status=0

#   For returning sequences NOT in the HitsList
#  		if matches.has_key(line[1:]):
# 			status=0
#  		else:
#  			OUT.write(str(line) + '\n')
#  			status = 1
	else:
		if status ==1:
			OUT.write(str(line) + '\n')