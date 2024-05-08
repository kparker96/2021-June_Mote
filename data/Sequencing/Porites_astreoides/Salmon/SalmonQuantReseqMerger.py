#!/usr/bin/env python

import sys
#sys.argv[1] file with list of your paired files to combine

columntoextract = 4

def make_dict1(file):
	fin = open(file, 'r')
	dict={}
	count=0
	for line in fin:
		count+=1
		line=line.rstrip()
		cols=line.split('\t') #for tab-delimited text files
		if count > 1:
			dict[cols[0]]=float(cols[columntoextract])
	fin.close()
	print('Read in first file')
	return dict
	
InfileList=open(sys.argv[1], 'r')
FileCount=0
for line in InfileList:
	FileCount+=1
	line=line.rstrip()
	Files=line.split('\t')
	FirstFileDict=make_dict1(Files[0])
	SecondFile=open(Files[1], 'r')
	SecondFileLine=0
	Outfile=open('%s_merged.txt'%(Files[0][:-3]), 'w')
	for Record in SecondFile:
		SecondFileLine+=1
		Record=Record.rstrip()
		Items=Record.split('\t')
		if SecondFileLine==1:
			Outfile.write('ContigName\tNumReads\r')
		if SecondFileLine>1:
			Newcount=float(Items[columntoextract])+FirstFileDict[Items[0]]
			Outfile.write('%s\t%.3f\r'%(Items[0],Newcount))
	SecondFile.close()
	Outfile.close()
	print('Merged %s and %s into %s' %(Files[0], Files[1],'%s_merged.txt'%(Files[0][:-3])))
		

