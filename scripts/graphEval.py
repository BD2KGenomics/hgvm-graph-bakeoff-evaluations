#!/usr/bin/env python
from __future__ import print_function,division, unicode_literals
import urllib, urllib2,requests, json, argparse, sys, glob, string, codecs, os, socket, time, contextlib, tarfile
from collections import OrderedDict, defaultdict
from itertools import izip


def getReferenceID(url):
	"""
	Returns a reference set ID from the server at the given URL.
	"""
	url+='/referencesets/search'
	req={
      "accessions": [],
      "assemblyId": '',
      "md5checksums": [],
      "pageSize": 100, 
      "pageToken": '0'
	}
	header={'Content-Type':'application/json'}
	res=requests.post(url,data=json.dumps(req, encoding='utf-8'),headers=header)
	thePage=res.text
	thePage=json.loads(thePage)
	#If the server returned no reference sets, just use '0'
	try:
		referenceSetId=thePage['referenceSets'][0]['id']
	except IndexError:
		referenceSetId='0'
	return referenceSetId


def getSequenceDict(url,referenceSetId=0,wantBases=False):
	sequenceDict={}
	baseDict={}
	url+='/sequences/search'
	req={
      "pageSize": 100,
      "pageToken": '0',
      "referenceSetId": referenceSetId, 
      "variantSetId": "null",
      "listBases":wantBases}
	nextPageToken=True
	header={'Content-Type':'application/json'}
	while nextPageToken:
		res=requests.post(url,data=json.dumps(req),headers=header)
		thePage=json.loads(res.text)
		nextPageToken=thePage['nextPageToken']
		req['pageToken']=nextPageToken
		sequenceList=thePage['sequences']
		for sequence in sequenceList:
			id_=sequence['id']
			length=int(sequence['length'])
			sequenceDict[id_]=length
			if wantBases:
				bases=sequence['bases']
				baseDict[id_]=bases
	if wantBases:
		return sequenceDict, baseDict
	else:
			return sequenceDict


def getJoinDict(url,referenceSetId=0):
	"""
	Returns a dict containing all the joins in the specified server.
	
		Key: seq
		Value: {(pos,strand):[(seq,pos,strand),(seq,pos,strand),...]}
	"""
	def defaultList():
		return defaultdict(list)
	joinDict=defaultdict(defaultList)
	url+='/joins/search'
	req={
      "length": None, 
      "pageSize": 100, 
      "pageToken": '0', 
      "referenceSetId": str(referenceSetId), 
      "sequenceId": '', 
      "start": None, 
      "strand": None, 
      "variantSetId": ''
	}
	nextPageToken=True
	header={'Content-Type':'application/json'}
	while nextPageToken:
		res=requests.post(url,data=json.dumps(req),headers=header)
		thePage=json.loads(res.text)
		nextPageToken=thePage['nextPageToken']
		req['pageToken']=nextPageToken
		joins=thePage['joins']
		for join in joins:
			side1=join['side1']
			side2=join['side2']
			base1=side1['base']
			base2=side2['base']
			seq1=base1['sequenceId']
			seq2=base2['sequenceId']
			pos1=int(base1['position'])
			pos2=int(base2['position'])
			strand1=side1['strand']
			strand2=side2['strand']
			joinDict[seq1][(pos1,strand1)].append((seq2,pos2,strand2))
			joinDict[seq2][(pos2,strand2)].append((seq1,pos1,strand1))
	return joinDict

def getAlleleIDDict(url):
	"""
	Uses the alleles/search endpoint of the ga4gh server to get a dict
	of ID:name allele pairs.
	"""
	alleleIDDict={}
	url+='/alleles/search'
	req={	
		"pageSize": 100,
		"pageToken": '0', 
		"start": 0, 
		"end": 10, 
		"sequenceId": "", 
		"variantSetIds": []
		}
	nextPageToken=True
	header={'Content-Type':'application/json'}
	while nextPageToken:
		res=requests.post(url,data=json.dumps(req),headers=header)
		thePage=json.loads(res.text)
		nextPageToken=thePage['nextPageToken']
		req['pageToken']=nextPageToken
		alleles=thePage['alleles']
		for allele in alleles:
			id_=allele['id']
			name=allele['name']
			alleleIDDict[id_]=name
	return alleleIDDict

def getAlleleDict(url,alleleIDDict):
	"""
	Sends a series of get-allele requests to the server url specified by the user.  
	Requests the allele IDs contained in the keys of alleleIDDict. 
	Returns all alleles in the dictionary alleleDict, in which allele names are the key and
	allele path items (stored as a list of dicts) are values.
	"""
	print("Getting alleles from server...")
	alleleDict=OrderedDict()
	def getAlleles(url,tries):
		if tries==5:
			return 'Tries'
		try:
			return urllib2.urlopen(url).read()
		except urllib2.HTTPError:
			return 'Fail'
		except socket.error:
			print('Server timed out.  Trying again after 5 second break.')
			time.sleep(5)
			return getAlleles(url,tries+1)

	for id_,name in alleleIDDict.items():
			text=getAlleles(url+'/alleles/'+str(id_),0)
			if text=='Tries':
				print('Server timed out 5 times.  Skipping server.')
				return False
			elif text=='Fail':
				print('Error retrieving allele {} from server.  Skipping server.'.format(id_))
				return False

			json_=json.loads(text)
			id_=int(json_['id'])
			name=json_['name']
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
	return alleleDict

def maf2Indices(inFile):
	"""
	Extracts from a maf file a dictionary containing
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
	with open(inFile) as inFile:
		for line in inFile:
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
		        refDict[seq][-1] = (refDict[seq][-1][0], max(end,refDict[seq][-1][1]))
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
	with open(bedFile) as inFile:
		for line in inFile:
			line=line.strip().split()
			start=int(line[1])
			end=int(line[2])
			gene=line[3]
			for index in range(start,end+1):
				geneMap[index].add(gene)
	return geneMap

def countJoins(joinDict):
	"""
	Counts the joins in joinDict, which is structured as follows:	
		Key: seq
		Value: {(pos,strand):[(seq,pos,strand),(seq,pos,strand),...]}
	"""
	joinCount=0
	for seq in joinDict:
		for tuple_ in joinDict[seq]:
			for otherSide in joinDict[seq][tuple_]:
				joinCount+=1
	joinCount/=2
	return joinCount


def countSegments(sequenceDict,joinDict):
	"""
	Count all the segments in the graph by splitting up
	sequences at points where joins attach to them.  Return the
	count as an integer.
	"""
	segmentCount=len(sequenceDict)
	#A set of (int,int) tuples marking the (seq,index) of locations
	#of split-points of sequences.  Each of these will increase the 
	#total number of segments by one.
	#For example: (1,0) means there is a split in sequence 1 AFTER
	#the first (index 0) base.
	#Indices in the second position should be at most length-1.
	nonEndJoinSet=set()
	for seq in joinDict:
		if seq!='0':
			for pos,strand in joinDict[seq]:
				length=sequenceDict[seq]
				if not ((pos==0 and strand=='POS_STRAND') or (pos==length-1 and strand=='NEG_STRAND') or length==1):
					if strand=='NEG_STRAND':
						nonEndJoinSet.add((seq,pos))
					else:
						nonEndJoinSet.add((seq,pos-1))
	segmentCount+=len(nonEndJoinSet)
	return segmentCount

@contextlib.contextmanager
def smartOpen(filename=None):
	if filename and filename != '-':
		fh = open(filename, 'w')
	else:
		fh = sys.stdout

	try:
		yield fh
	finally:
		if fh is not sys.stdout:
			fh.close()


def parseArgs():
	parser = argparse.ArgumentParser(description="""Performs the specified evaluation on a specified graph server.
		Requires a url to be supplied from the user.""")
	evaluation=parser.add_mutually_exclusive_group()
	target=parser.add_mutually_exclusive_group()
	target.add_argument('--list',help="""A tab-separated file, with a header line, containing a list of all
		servers to be evaluated.  First four columns must be: region, url, algorithm, source.""")
	target.add_argument('--url',help="""The url of the graph to be evaluated.  Make sure url starts with "http://" and
		ends with the "/" preceding the server version.  This and the --list argument are mutually exclusive.""")
	evaluation.add_argument('--align2ref', action='store_true',
		help="""Compares each allele returned by the server to the reference allele, returning a percent overlap.
		Assumes that the reference allele is named "ref", or "ref.ref".""")
	evaluation.add_argument('--maf', action='store_true', help="""Compares the graph-alignments of each allele to the reference, 
		to the corresponding graph-alignments in a separate MAF file.  Requires the name of the maf file.  
		Assumes the reference allele is named 'ref' or 'ref.ref'.  Also currently assumes the MAF only contains
		a single block.""")
	evaluation.add_argument('--gene',action='store_true',help="""For each alignment-column in a graph, and given
		a directory containing a bed file for each path in the graph, computes the number of ortholog and paralog alignments.""")
	parser.add_argument('--stats',action='store_true',help="""Compute general stats about the graph, such as number of positions and number of seqs/joins.""")
	parser.add_argument('--out',help="""The name of an output file.""")
	args = parser.parse_args()
	return args



def main():
	args=parseArgs()
	if args.list:
		serverDict=defaultdict(list)
		with open(args.list) as inFile:
			for line in inFile:
				if line and not line.startswith('#'):
					line=line.strip().split('\t')
					try:
						region,url,algo,source=line[:4]
					except:
						print(line)
					if region!='region':
						serverDict[region].append([url,algo,source])
	elif args.url:
		serverDict={None:[[args.url,None,None]]}
	else:
		sys.exit('Please specify a file listing servers to be evaluated using the --list or --url options.')
	
	if args.align2ref:
		with smartOpen(args.out) as outFile:
			for region in serverDict:
				if region!='cenx':
					if args.list:
						outFile.write('$'+region+'\n')
					for url,algo,source in serverDict[region]:
						url+='v0.6.g'
						print("Processing "+url)
						#Get a dict of id:name pairs
						alleleIDDict=getAlleleIDDict(url)
						if 'ref' not in alleleIDDict.values() and 'GRCh38_2c5:0' not in alleleIDDict.values() and 'GRCh38_247:0' not in alleleIDDict.values():
							print("Skipping... graph does not contain an allele called either 'ref' or 'GRCh38_2c5:0' or 'GRCh38_247:0'.")
							continue
						else:
							alleleDict=getAlleleDict(url,alleleIDDict)
							if not alleleDict:
								continue
							if args.list:
								outFile.write('@'+algo+' '+source+'\n')
							#Special case for curoverse, since both brca1 and brca2 are in the same graph
							if algo=='curoverse':
								if region=='brca1':
									total=0
									for allelePathItemList in alleleDict.values():
										total+=sum([pathItem['length'] for pathItem in allelePathItemList])
									print(total)
									break
									refAllele=alleleDict['GRCh38_2c5:0']
									alleleDict={'GRCh38_2c5:0':alleleDict['GRCh38_2c5:0'],
												'GI262359905_rc:0':alleleDict['GI262359905_rc:0'],
												'GI528476558:0':alleleDict['GI528476558:0']}
								elif region=='brca2':
									total=0
									for allelePathItemList in alleleDict.values():
										total+=sum([pathItem['length'] for pathItem in allelePathItemList])
									print(total)
									break
									refAllele=alleleDict['GRCh38_247:0']
									alleleDict={'GRCh38_247:0':alleleDict['GRCh38_247:0'],
												'GI388428999:0':alleleDict['GI388428999:0'],
												'GI528476586:0':alleleDict['GI528476586:0']}
								else:
									print("Curoverse servers in non-brca regions not supported.  Skipping.")
									continue
							else:
								refAllele=alleleDict['ref']
							#Make sure the server only has one reference allele
							if len([name for name in alleleDict.keys() if name in ['ref','GRCh38_2c5:0','GRCh38_247:0']])>1:
								print("Skipping because server has more than one reference allele.")
								continue
							#Make sure the server has more than just a reference allele
							if len([name for name in alleleDict.keys() if name not in ['ref','GRCh38_2c5:0','GRCh38_247:0']])==0:
								print("Skipping because server only has reference allele.")
								continue
							refDict=mergeRefItems(refAllele)
							print("Computing overlap with reference allele...")
							totalOverlap=0
							for name,allele in sorted(alleleDict.items()):
								if not name in ['ref','GRCh38_2c5:0','GRCh38_247:0']:
									refOverlapFraction=getRefOverlap(allele,refDict)
									totalOverlap+=refOverlapFraction
									outFile.write("{}\t{}\n".format(name,refOverlapFraction))
							averageOverlap=totalOverlap/(len(alleleDict)-1)
							outFile.write('average\t'+str(averageOverlap)+'\n')


	elif args.maf:
		#Check for the existence of maf files in the data directory
		mafDir='../data/alignments/'
		if not os.path.isdir(mafDir):
			if not os.path.exists("../data/alignments.tar.gz"):
				sys.exit("This evaluation needs 'alignments.tar.gz' in the data directory.")
			else:
				os.chdir('../data')
				tf=tarfile.open('alignments.tar.gz')
				tf.extractall()
				tf.close()
				os.chdir('../scripts')
		#Get maf alt dicts for all regions before opening any servers
		#Store in mafMasterDict
		lrc_kirMaf=mafDir+"LRC_KIR_GRCAlignment.maf"
		smaMaf=mafDir+"SMA_GRCAlignment.maf"
		mhcMaf=mafDir+"MHC_GRCAlignment.maf"
		mafMasterDict={
			'lrc_kir':maf2Indices(lrc_kirMaf),
			'mhc':maf2Indices(mhcMaf),
			'sma':maf2Indices(smaMaf)
		}
		with smartOpen(args.out) as outFile:
			for region in ['sma','mhc','lrc_kir']:
				outFile.write('$'+region+'\n')
				for url,algo,source in serverDict[region]:
					url+='v0.6.g'
					print("Processing "+url)
					
					#Get a dict of id:name pairs
					alleleIDDict=getAlleleIDDict(url)
					if 'ref' not in alleleIDDict.values() and 'GRCh38_2c5:0' not in alleleIDDict.values():
						print("Skipping... graph does not contain an allele called either 'ref' or 'GRCh38_2c5:0'.")
						continue
					else:
						alleleDict=getAlleleDict(url,alleleIDDict)

						#Check if the server has any alleles (and if not, skip to the next one)
						if not alleleDict:
							continue

						#Make sure the server only has one reference allele
						if len([name for name in alleleDict.keys() if name in ['ref','GRCh38_2c5:0','GRCh38_247:0']])>1:
							print("Skipping because server has more than one reference allele.")
							continue
						#Make sure the server has more than just a reference allele
						if len([name for name in alleleDict.keys() if name not in ['ref','GRCh38_2c5:0','GRCh38_247:0']])==0:
							print("Skipping because server only has reference allele.")
							continue
						graphAltDict=graph2Indices(alleleDict)
						mafAltDict=mafMasterDict[region]
						if set(mafAltDict.keys())!=set(graphAltDict.keys()):
							print("Skipping because maf alleles:\n{}\naren't the same as graph alleles:\n{}".format(mafAltDict.keys(),graphAltDict.keys()))
							continue
						stop=False
						for alt in mafAltDict:
							if len(mafAltDict[alt])!=len(graphAltDict[alt]):
								print("""Skipping because corresponding {} alleles in maf and graph don't have the same length.
									Maf alt has length {} while graph alt has length {}.""".format(alt,len(mafAltDict[alt]),len(graphAltDict[alt])))
								stop=True
						if stop:
							continue
						print("Computing precision and recall...\n")
						totalMatchCount=0
						totalNumMafAlignedBases=0
						totalNumGraphAlignedBases=0
						for alt in mafAltDict:
							if alt!='ref':
								mafAlt=mafAltDict[alt]
								graphAlt=graphAltDict[alt]
								matchCount=sum(map(lambda n:1 if n[0]!=0 and n[0]==n[1] else 0,izip(mafAlt,graphAlt)))
								numMafAlignedBases=len([i for i in mafAlt if i!=0])
								numGraphAlignedBases=len([i for i in graphAlt if i!=0])
								totalMatchCount+=matchCount
								totalNumMafAlignedBases+=numMafAlignedBases
								totalNumGraphAlignedBases+=numGraphAlignedBases
						try: avePrecision=totalMatchCount/totalNumGraphAlignedBases
						except ZeroDivisionError:
							avePrecision=0
							print("Skipping due to ZeroDivisionError when computing average precision.")
							continue
						try: aveRecall=totalMatchCount/totalNumMafAlignedBases
						except ZeroDivisionError:
							aveRecall=0
							print("Skipping due to ZeroDivisionError when computing average recall.")
							continue
						outFile.write('@'+algo+'\t'+source+'\n')
						outFile.write("average\t{}\t{}\n".format(avePrecision,aveRecall))
	

	elif args.gene:
		with smartOpen(args.out) as outFile:
			for region in ['sma','mhc','lrc_kir']:
				outFile.write('$'+region+'\n')
				bedDir="../data/genes/"+region.upper()+'/'
				if not os.path.isdir(bedDir):
					if not os.path.exists("../data/genes.tar.gz"):
						sys.exit("This evaluation needs 'genes.tar.gz' in the data directory.")
					else:
						os.chdir('../data')
						tf=tarfile.open('genes.tar.gz')
						tf.extractall()
						tf.close()
						os.chdir('../scripts')

				print("Reading reference bed file for {}...".format(region))
				refGeneDict=getGenesFromBed(bedDir+'ref/genes.bed')
				for url,algo,source in serverDict[region]:
					url+='v0.6.g'
					print("Processing "+url)
					
					#Get a dict of id:name pairs
					alleleIDDict=getAlleleIDDict(url)
					if 'ref' not in alleleIDDict.values() and 'GRCh38_2c5:0' not in alleleIDDict.values():
						print("Skipping... graph does not contain an allele called either 'ref' or 'GRCh38_2c5:0'.")
						continue
					else:
						alleleDict=getAlleleDict(url,alleleIDDict)

						#Check if the server has any alleles (and if not, skip to the next one)
						if not alleleDict:
							continue

						#Make sure the server only has one reference allele
						if len([name for name in alleleDict.keys() if name in ['ref','GRCh38_2c5:0','GRCh38_247:0']])>1:
							print("Skipping because server has more than one reference allele.")
							continue
						#Make sure the server has more than just a reference allele
						if len([name for name in alleleDict.keys() if name not in ['ref','GRCh38_2c5:0','GRCh38_247:0']])==0:
							print("Skipping because server only has reference allele.")
							continue
						outFile.write('@'+algo+'\t'+source+'\n')
						graphAltDict=graph2Indices(alleleDict)
						totalOrthologCount=0
						totalParalogCount=0
						for allele,indexList in graphAltDict.iteritems():
							altGeneDict=getGenesFromBed(bedDir+allele+'/genes.bed')
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
							totalOrthologCount+=orthologCount
							totalParalogCount+=paralogCount
							# print("Number bases with orthologous ref mappings:"+str(orthologCount))
							# print("Number bases with paralogous ref mappings:"+str(paralogCount))
							outFile.write(allele+'\t'+str(orthologCount)+'\t'+str(paralogCount)+'\n')
						outFile.write("total\t"+str(totalOrthologCount)+'\t'+str(totalParalogCount)+'\n')


	elif args.stats:
		with smartOpen(args.out) as outFile:
			for region in serverDict:
				outFile.write('$'+region+'\n')
				for url,algo,source in serverDict[region]:
					url+='v0.6.g'
					print("Processing "+url)
					print('Getting referenceSetId...')
					referenceSetId=getReferenceID(url)
					print('Getting sequences...')
					sequenceDict=getSequenceDict(url,'0')
					print('Getting joins...')
					joinDict=getJoinDict(url,'0')
					print('Computing stats...')
					# Report general features of the graph:
						# Number of sequences
					print("There are {} sequences in the graph".format(len(sequenceDict)))
						# Number of joins
					joinCount=countJoins(joinDict)
					print("There are {} joins in the graph".format(joinCount))
						# Number of positions
					positionSum=sum(sequenceDict.values())
					print("There are {} positions in the graph.".format(positionSum))
					# 	Number of segments (split up by joins)
					segmentCount=countSegments(sequenceDict,joinDict)
					print("There are {} segments in the graph.".format(segmentCount))
					# 	Number of joins (in split-up graph)
					fullJoinCount=joinCount+(segmentCount-len(sequenceDict))
					print("There are {} joins in the segmentized graph.".format(fullJoinCount))
					outFile.write(algo+'\t'+str(positionSum)+'\n')

	else:
		print("Please specify an evaluation option.")

if __name__ == '__main__':
	main()

