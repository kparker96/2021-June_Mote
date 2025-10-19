#!/usr/bin/env python

Usage = """
DEcalc.py - version 1.0
Calculate the number of Differentially Expressed Genes for each contrast 
and the up or down regulation of DEGs.
Usage:
	DEcalc.py *.txt
"""

# import modules 
import sys # lets us access system-specific parameters and functions

sig = float(0.10)

if len(sys.argv)<1:
	print(Usage) #print usage statements to the screen
else:
	FileList = sys.argv[1:]
	FileNumber = 0
	CombinedOutFileName = "2025-06-03_DECalc_deduped-data.tab"
	CombinedOutFile = open(CombinedOutFileName, 'w')
	for InFileName in FileList:
		InFile = open(InFileName, 'r') # open the infile
		Contrast = InFileName[:-4]
		Headerline = "Contrast,NumDEGsUp,NumDEGsDown,TotalNumDEGs"
		if FileNumber == 0:
			CombinedOutFile.write(Headerline + '\n')
		LineNumber = 0
		TotalNumDEGs = 0
		NumDEGsUp = 0
		NumDEGsDown = 0
		for Line in InFile:
			LineNumber+=1 # Line Number + 1
			if LineNumber > 1 :
				Line = Line.strip("\n").strip("\r")
				List = Line.split('\t')
				Contig = List[0]
				baseMean = List[1]
				log2FoldChange = float(List[2])
				lfcSE = List[3]
				stat = List[4]
				if List[5] == 'NA' or List[6] == 'NA':
					continue 
				else:
					pvalue = float(List[5])
					padj = float(List[6])
				
				# "Contrast,NumDEGsUp,NumDEGsDown,TotalNumDEGs"
				if padj < sig:
					TotalNumDEGs += 1
					if log2FoldChange > 0:
						NumDEGsUp += 1
					elif log2FoldChange < 0:
						NumDEGsDown += 1

		CombinedOutFile.write('%s\t%d\t%d\t%d\n' %(Contrast,NumDEGsUp,NumDEGsDown,TotalNumDEGs))
		InFile.close()
		FileNumber+= 1
		print("finished processing file: %s"%(InFileName)) 