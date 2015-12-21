#!/usr/bin/env python
from __future__ import division,print_function
from collections import defaultdict
import random, sys, argparse

def getRefSeq(fastaFile):
	"""
	Returns the first sequence in a fasta file.
	"""
	print("Getting ref sequence...")
	with open(fastaFile) as inFile:
		seq=''
		for line in inFile:
			line=line.strip()
			if line.startswith('>'):
				if seq:
					print("Found refSeq of length {}".format(len(seq)))
					break
			else:
				seq+=line.upper()
	return seq


def trimEnds(alt,ref):
	"""
	Takes two strings and trims off any starting and ending sequence
	that is the same in both.  Returns both trimmed sequences and
	the coordinates of the first and last different characters.
	"""
	trimmedAlt,trimmedRef='',''
	minLen=min([len(alt),len(ref)])
	firstDif=0
	lastDif=0
	for altChar,refChar in zip(alt,ref):
		if altChar==refChar:
			firstDif+=1
		else:
			break
	for altChar,refChar in zip(alt[::-1],ref[::-1]):
		if firstDif+lastDif==minLen:
			break
		if altChar==refChar:
			lastDif+=1
		else:
			break
	trimmedAlt=alt[firstDif:len(alt)-lastDif]
	trimmedRef=ref[firstDif:len(ref)-lastDif]

	return trimmedAlt,trimmedRef,firstDif,lastDif



def editVCF(inFile,outFile,refSeq,shift,start,stop):
	"""
	Reads in a vcf file and creates a new, modified vcf file.
	Shifts forward the position of all entries by <shift> bp 
	(restarting at the beginning if it exceeds the length of the reference)
	"""
	print("Editing vcf file...")

	with open(outFile,'w') as outFile:
		with open(inFile) as inFile:
			for line in inFile:
				if not line.startswith('#'):
					line=line.split('\t')

					#Get pos,ref,and alt
					pos=int(line[1])
					ref=line[3].upper()
					altList=line[4].upper().split(',')


					#Change from 1-based to 0-based
					pos=pos-1


					#Circular shift by shift
					if pos+shift+len(ref)-1>stop:
						overshoot=pos+shift+len(ref)-1-stop
						pos=start+overshoot-1
					else:
						pos=pos+shift
					newRef=refSeq[pos:pos+len(ref)]


					SNPList=['A','C','T','G']
					takenSNPList=[]

					#Classify all alts
					for index, alt in enumerate(altList):

						#Trim ends of alt and ref
						trimmedAlt,trimmedRef,firstDif,lastDif=trimEnds(alt,ref)
						newTrimmedRef=newRef[firstDif:len(ref)-lastDif]


						#Classify variant

						#SNP
						if len(trimmedAlt)==len(trimmedRef)==1:
							takenSNPList.append(newTrimmedRef)
							newSNP=random.choice([base for base in SNPList if base not in takenSNPList])
							takenSNPList.append(newSNP)
							newAlt=newRef[:firstDif]+newSNP+newRef[len(newRef)-lastDif:]

						#MNPs, Insertions, and Deletions
						else:
							newAlt=newRef[:firstDif]+trimmedAlt+newRef[len(newRef)-lastDif:]


						altList[index]=newAlt

					#Change back to 1-based
					pos=pos+1

					#Update vcf line
					line[1]=str(pos)
					line[3]=newRef
					line[4]=','.join(altList)
					line='\t'.join(line)

				outFile.write(line)

def parseArgs():
	parser = argparse.ArgumentParser(description="""Takes an uncompressed vcf file and circularly shifts all the 
		variants right by <shift> number of nucleotides.  Requires a fasta file containing the reference sequence.
		If you wish to perform a circular shift within a substring of the ref, you must use the start and/or stop
		arguments.""")

	#Required arguments
	parser.add_argument('inFile',help="""The vcf file to be shifted.""")
	parser.add_argument('--ref',required=True,help="""A fasta file containing the reference sequence of the vcf file.  Required.""")
	parser.add_argument('--shift',required=True,type=int,help="""The number of bases to shift.  Required.""")

	#Optional arguments
	parser.add_argument('--start',type=int,default=0,help="""The first index of the region in the ref.""")
	parser.add_argument('--stop',type=int,help="""The last index of the region in the ref.""")
	parser.add_argument('--outFile',help="""The output vcfFile.  Default is <inFile_root>_shifted.vcf""")

	args = parser.parse_args()
	return args

def main():
	#Parse args for input vcf and ref info
	args=parseArgs()
	inFile=args.inFile
	if not args.outFile:
		outFile='.'.join(inFile.split('.')[:-1])+'_shifted.vcf'

	#Get ref file
	refFile=args.ref
	refSeq=getRefSeq(refFile)

	#Get shift, and start/stop indices
	shift=args.shift
	start=args.start
	stop=args.stop if args.stop else len(refSeq)-1

	#Obtain the reference sequence to ensure we don't overshoot the pilot region when shifting
	#and also to make sure that the reference sequence at the shifted position isn't the same
	#as the variant sequence.
	editVCF(inFile,outFile,refSeq,shift,start,stop)


if __name__ == "__main__":
	main()