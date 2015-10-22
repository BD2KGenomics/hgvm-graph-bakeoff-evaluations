#!/usr/bin/env python
from __future__ import print_function,division
import urllib2, json, argparse, sys, glob
from collections import OrderedDict, defaultdict
from itertools import izip

def getAlleles(url):
	"""
	Sends a series of get-allele requests to the server url specified by the user, starting
	with an allele id of 0 and incrementing the id until the get request fails.  Returns
	all alleles in the dictionary alleleDict, in which allele names are the key and
	allele path items (stored as a list of dicts) are values.

	Currently assumes that all alt IDs are greater than the ref ID, and only starts
	storing alleles once the ref allele is reached.  This is so that ancestral
	alleles (which somtimes precede the ref) are skipped.
	"""
	print("Getting alleles from server...")
	alleleDict=OrderedDict()
	success=True
	alleleNumber=-1
	refReached=False
	while success==True:
		alleleNumber+=1
		try:
			text=urllib2.urlopen(url+'/alleles/'+str(alleleNumber)).read()
			json_=json.loads(text)
			id_=int(json_['id'])
			name=json_['name'].split('.')[-1]
			print(name)
			if name in ['ref','GRCh38:0','GRCh38_2c5:0']:
				refReached=True
				name='ref'
			# elif not name.startswith('hu'):
			# 	pass
			# else:
			# 	break
			if refReached: #and not name.startswith('hu'):
				variantSetID=json_['variantSetId']
				path=json_['path']
				segments=path['segments']
				allelePathItemList=[]
				for segment in segments:
					allelePathItem={}
					length=int(segment['length'])
					start=segment['start']
					strand=start['strand']
					base=start['base']
					pos=int(base['position'])
					seqId=int(base['sequenceId'])
					allelePathItem['seq']=seqId
					allelePathItem['strand']=strand
					allelePathItem['pos']=pos
					allelePathItem['length']=length
					allelePathItemList.append(allelePathItem)
				alleleDict[name]=allelePathItemList
		except urllib2.HTTPError:
			if alleleNumber>0:
				success=False
	return alleleDict

def maf2Indices(inputFile):
	"""
	Extracts from the maf file input by the user a dictionary containing
	allele names as keys. Values are lists of integers with lengths equal to
	the lengths of the allele sequences; integers correspond to the 1-based index of
	the reference base that the allele base is aligned to (or 0 if the
	base is not aligned.)  If the alt base aligns to the '-' strand of the ref,
	the integer will be negative.

	ToDo: This program currently assumes that there is only one
	block per MAF file, which is not generally true of the format,
	but happens to be true of the GRC maf files.
	"""
	print("Converting maf file into alignment indices...")
	altDict=defaultdict(list)
	cursorDict={}
	infoList=[]
	sequenceList=[]
	seqID=-1
	for line in inputFile:
		line=line.strip().split()
		if line and line[0]=='s':
			seqID+=1
			name,start,length,strand,sourceLength,sequence=line[1:]
			start=int(start)
			sourceLength=int(sourceLength)
			infoList.append([name,strand,seqID])
			sequenceList.append(sequence)
			if strand=='+':
				cursorDict[seqID]=start
			elif strand=='-':
				cursorDict[seqID]=sourceLength-start-1
			else:
				raise Exception("strand is not '+' or '-'.")
			if name not in altDict and name!='ref':
				altDict[name]=[0]*sourceLength
	refIndex=0
	for letterList in izip(*sequenceList):
		if letterList[0]!='-':refIndex+=1
		for info,altLetter in izip(infoList[1:],letterList[1:]):
			name,strand,seqID=info
			altIndex=cursorDict[seqID]
			if letterList[0]=='-':
				if altLetter=='-':
					pass
				else:
					if strand=='-':
						cursorDict[seqID]-=1
					else:
						cursorDict[seqID]+=1
			else:
				if altLetter=='-':
					pass
				else:
					if strand=='-':
						altDict[name][altIndex]=-refIndex
						cursorDict[seqID]-=1
					else:
						altDict[name][altIndex]=refIndex
						cursorDict[seqID]+=1
	return altDict

def graph2Indices(alleleDict):
	"""
	Extracts from the graph server alleles input by the user a dictionary 
	containing allele names as keys. Values are lists of integers with lengths 
	equal to the lengths of the allele sequences; integers correspond to the 1-based 
	index of the reference base that the allele base is aligned to (or 0 if the
	base is not aligned.)  If the alt base aligns to the '-' strand of the ref,
	the integer will be negative.
	"""
	print("Converting graph alleles into alignment indices...")
	altDict=defaultdict(list)
	refSegMap=defaultdict(dict)
	refIndex=1
	#Record which sequence bases map to which ref bases
	for pathItem in alleleDict['ref']:
		seqID=pathItem['seq']
		strand=pathItem['strand']
		start=pathItem['pos']
		length=pathItem['length']
		if strand=='POS_STRAND':
			for seqIndex in xrange(start+1,start+length+1):
				refSegMap[seqID][seqIndex]=refIndex
				refIndex+=1
		else:
			for seqIndex in xrange(-start-1,-start+length-1):
				refSegMap[seqID][seqIndex]=refIndex
				refIndex+=1
	for name in alleleDict:
		if name!='ref':
			for pathItem in alleleDict[name]:
				seqID=pathItem['seq']
				strand=pathItem['strand']
				start=pathItem['pos']
				length=pathItem['length']
				if seqID in refSegMap:
					if strand=='POS_STRAND':
						for seqIndex in xrange(start+1,start+length+1):
							if seqIndex in refSegMap[seqID]:
								altDict[name].append(refSegMap[seqID][seqIndex])
							elif -seqIndex in refSegMap[seqID]:
								altDict[name].append(-refSegMap[seqID][-seqIndex])
							else:
								altDict[name].append(0)
					else:
						for seqIndex in xrange(-start-1,-start+length-1):
							if seqIndex in refSegMap[seqID]:
								altDict[name].append(refSegMap[seqID][seqIndex])
							elif -seqIndex in refSegMap[seqID]:
								altDict[name].append(-refSegMap[seqID][-seqIndex])
							else:
								altDict[name].append(0)
				else:
					altDict[name]+=[0]*length
	return altDict

def mergeRefItems(refPathList):
	"""
	Before computing the overlap between alt alleles and a reference
	allele, it helps to eliminate overlapping regions in the reference
	(like from duplicated segments.)

	Thus, this function takes the reference allele (a list of dicts), and
	returns a dict where sequence IDs (integers) are keys, and values are
	lists of tuples corresponding to base-index-ranges of each sequence that 
	the reference allele spans.
	"""
	print("Merging segments in reference allele...")
	unmergedRefDict=defaultdict(list)
	refDict=defaultdict(list)
	for pathItem in refPathList:
		seq=pathItem['seq']
		length=pathItem['length']
		strand=pathItem['strand']
		start=pathItem['pos']
		if strand=='POS_STRAND':
			end=start+length-1
		else:
			end=start-length+1
		unmergedRefDict[seq].append(sorted([start,end]))
	for seq,rangeList in unmergedRefDict.iteritems():
		for begin,end in sorted(rangeList):
		    if refDict[seq] and refDict[seq][-1][1] >= begin - 1:
		        refDict[seq][-1] = (refDict[seq][-1][0], end)
		    else:
		        refDict[seq].append((begin, end))
	return refDict

def getRefOverlap(allelePathItemList,refDict):
	"""
	Takes an alt allele (a list of dicts), and a ref dict 
	(a dict containing base-index-ranges for each Sequence...
	see function mergeRefItems.)

	Computes a percent overlap of the alt allele with the reference,
	i.e. number of bases in alt overlapping ref, divided by total 
	alt bases.
	"""
	overlapLength=0
	totalLength=0
	for pathItem in allelePathItemList:
		totalLength+=pathItem['length']
	for pathItem in allelePathItemList:
		segSeq=pathItem['seq']
		segStrand=pathItem['strand']
		segStart=pathItem['pos']
		segLength=pathItem['length']
		segEnd=segStart+segLength-1
		if segStrand=='NEG_STRAND':
			segEnd=segStart-segLength+1
		segStart,segEnd=sorted([segStart,segEnd])
		for refStart,refEnd in refDict[segSeq]:
			if not (segEnd<refStart or segStart>refEnd):
				if segStart>=refStart and segEnd<=refEnd:
					overlapLength+=segLength
				elif segStart<=refStart and segEnd>=refEnd:
					overlapLength+=refEnd-refStart+1
				elif segStart<=refStart and segEnd>=refStart:
					overlapLength+=segEnd-refStart+1
				elif segStart<=refEnd and segEnd>=refEnd:
					overlapLength+=refEnd-segStart+1
	refOverlapFraction=overlapLength/totalLength
	# print("Total length is {}".format(totalLength))
	# print("Total overlap length is {}".format(overlapLength))
	return refOverlapFraction

def getGenesFromBed(bedFile):
	"""
	Takes a directory of directories, each containing a single bed file, and returns
	a dict where keys are indices and values are lists of genes.
	"""
	geneMap=defaultdict(set)
	with open(bedFile) as inputFile:
		for line in inputFile:
			line=line.strip().split()
			start=int(line[1])
			end=int(line[2])
			gene=line[3]
			for index in range(start,end+1):
				geneMap[index].add(gene)
	return geneMap

def parseArgs():
	parser = argparse.ArgumentParser(description="""Performs the specified evaluation on a specified graph server.
		Requires a url to be supplied from the user.""")
	parser.add_argument('url',type=str,help="""A string containing the url of the graph server to be evaluated.
		  Must include a version number at the end of the url, e.g. 'v0.6.g' """)
	parser.add_argument('--align2ref', action='store_const', const=True,
		help="""Compares each allele returned by the server to the reference allele, returning a percent overlap.
		Assumes that the reference allele is named "ref", or "ref.ref".""")
	parser.add_argument('--maf',type=file, help="""Compares the graph-alignments of each allele to the reference, 
		to the corresponding graph-alignments in a separate MAF file.  Requires the name of the maf file.  
		Assumes the reference allele is named 'ref' or 'ref.ref'.  Also currently assumes the MAF only contains
		a single block.""")
	parser.add_argument('--gene',help="""For each alignment-column in a graph, and given
		a directory containing a bed file for each path in the graph, computes the number of ortholog and paralog alignments.""")
	args = parser.parse_args()
	return args

def main():
	args=parseArgs()
	if args.align2ref:
		alleleDict=getAlleles(args.url)
		refAllele=alleleDict['ref']
		refDict=mergeRefItems(refAllele)
		print("Computing overlap with reference allele...")
		for id_,allele in sorted(alleleDict.items()):
			refOverlapFraction=getRefOverlap(allele,refDict)
			print("{}\t{}".format(id_,refOverlapFraction))
	elif args.maf:
		alleleDict=getAlleles(args.url)
		graphAltDict=graph2Indices(alleleDict)
		mafAltDict=maf2Indices(args.maf)
		assert set(mafAltDict.keys())==set(graphAltDict.keys())
		for alt in mafAltDict:
			assert len(mafAltDict[alt])==len(graphAltDict[alt])
		print("Computing precision and recall...\n")
		totalMatchCount=0
		totalNumMafAlignedBases=0
		totalNumGraphAlignedBases=0
		for alt in mafAltDict:
			mafAlt=mafAltDict[alt]
			graphAlt=graphAltDict[alt]
			matchCount=sum(map(lambda n:1 if n[0]!=0 and n[0]==n[1] else 0,izip(mafAlt,graphAlt)))
			numMafAlignedBases=len([i for i in mafAlt if i!=0])
			numGraphAlignedBases=len([i for i in graphAlt if i!=0])
			totalMatchCount+=matchCount
			totalNumMafAlignedBases+=numMafAlignedBases
			totalNumGraphAlignedBases+=numGraphAlignedBases
			try: precision=matchCount/numGraphAlignedBases
			except ZeroDivisionError:
				precision=0
				print("ZeroDivisionError when computing precision for {}.".format(alt))
			try: recall=matchCount/numMafAlignedBases
			except ZeroDivisionError:
				recall=0
				print("ZeroDivisionError when computing recall for {}.".format(alt))
			print("Allele: {}\nPrecision: {:.3f}\nRecall: {:.3f}\n".format(alt,precision,recall))
		try: avePrecision=totalMatchCount/totalNumGraphAlignedBases
		except ZeroDivisionError:
			avePrecision=0
			print("ZeroDivisionError when computing average precision.")
		try: aveRecall=totalMatchCount/totalNumMafAlignedBases
		except ZeroDivisionError:
			aveRecall=0
			print("ZeroDivisionError when computing average recall.")
		print("Overall\nPrecision: {:.3f}\nRecall: {:.3f}".format(avePrecision,aveRecall))
	elif args.gene:
		alleleDict=getAlleles(args.url)
		graphAltDict=graph2Indices(alleleDict)
		refGeneDict=getGenesFromBed(args.gene+'/'+'ref'+'/genes.bed')
		for allele,indexList in graphAltDict.iteritems():
			altGeneDict=getGenesFromBed(args.gene+'/'+allele+'/genes.bed')
			# print(allele+":")
			#For each base in the alt
			#	check if there are any genes in the alt, and any genes in the index of the ref it's aligned to
			#	If there are genes in both
			#		If they're the same genes (or perhaps more), then count it as an ortholog alignment
			#		If they're different genes, then count it as a paralog alignment
			paralogCount=0
			orthologCount=0
			for altIndex,refIndex in enumerate(indexList):
				if refIndex<0:
					refIndex=-refIndex
				refIndex-=1
				if altIndex in altGeneDict and refIndex in refGeneDict:
					altGenes=altGeneDict[altIndex]
					refGenes=refGeneDict[refIndex]
					if altGenes&refGenes:
						orthologCount+=1
					else:
						paralogCount+=1
			# print("Number bases with orthologous ref mappings:"+str(orthologCount))
			# print("Number bases with paralogous ref mappings:"+str(paralogCount))
			print(allele+'\t'+str(orthologCount)+'\t'+str(paralogCount))

	else:
		print("Please specify an evaluation option.")


if __name__ == '__main__':
	main()

