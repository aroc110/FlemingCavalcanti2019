def parse_fasta(fasta_file):
	"""Opens fasta assembly file and makes a dictionary for assembly regions with region name as the key 
	for the sequence itself"""
	from Bio import SeqIO
	fasta_dic ={}
	
	with open(fasta_file, 'r') as genome: #parse FASTA file
		for record in SeqIO.parse(genome, "fasta"):
			if record.id not in fasta_dic:
				fasta_dic[record.id] = str(record.seq).upper() #make dictionary for assembly regions
	
	return fasta_dic
	
def process_gff(gff_file,marker,gfftype):
	"""Reads gff, produces a dataframe from the information,
	filters out uneccessary columns and rows, combines cDS 
	of the same genes, and reduces attributes to a simplified 
	gene name"""
	import pandas as pd
	import re
	
	df = pd.read_table(gff_file, sep='\t',names=["seqid", "source", "type","start","end","score","strand","phase","attributes"], dtype={'seqid': str, 'source': str, 'type': str, 'start': str, 'end': str, 'score': str, 'strand': str, 'phase': str, 'attributes': str})

	start_end = {}
	
	df = df.loc[df['type'] == gfftype] #filter to only CDS or exon depending on the gff convention used
	df = df[pd.notnull(df['end'])] #eliminate null rows
	df = df.drop(['source','score','phase'], axis=1) #drop uneccesary columns
	
	regex = re.compile(marker+'([^;,]*)')
	
	for index, row in df.iterrows():
		found = regex.findall(str(row['attributes']))
		if found:
			row['attributes'] = found[0]
		else:
			#print(gff_file,"\tRemoved: Not found  ",marker,row['attributes'])
			continue
		
		if row['attributes'] not in start_end: #produce dictionary with seqid as the key
			start_end[row['attributes']] = [row['seqid'],row['type'],row['start'],row['end'],row['strand'],row['attributes']]
		if int(row['end']) > int(start_end[row['attributes']][3]): #replaces the 'end' of the gene with a larger 'end' location from the next exon
			start_end[row['attributes']][3] = row['end']
		if int(row['start']) < int(start_end[row['attributes']][2]): #replaces the 'start' of the gene with a smaller 'start' location from the next exon
			start_end[row['attributes']][2] = row['start']
   
	df_new = pd.DataFrame(start_end) #dictionary to dataframe
	df_new = df_new.transpose() #rotate dataframe
	df_new.columns = ['seqid','type','start','end','strand','attributes']
	
	return df_new
	
		
def get_UTRs(org,df,fasta_dic,genCode,N,utrFile):
	"""Extract gene sequences + N nt downstream based on locations
	from GFF files (now in dataframe format), and write files
	containing all UTRS and sequences for which there were no
	stop codons."""
	import pandas as pd
	
	if genCode == 1:
		stops = ['TAA','TAG','TGA']
	elif genCode == 6:
		stops = ['TGA']
	elif genCode == 10:
		stops = ['TAA','TAG']
	else:
		print("Problem with genetic code for ", org)
	
	UTRs = []
	stopList = []
	good = 0
	bad = 0
	total = 0
	badStart = 0
	badStop = 0
	badSmall = 0
	badN = 0
	missingDir = 0
	
	TAG = 0
	TAA = 0
	TGA = 0
	
	for index, row in df.iterrows(): #iterates over rows of the distilled table

		total += 1
		location = str(row['seqid'])
		start = int(row['start'])
		end = int(row['end'])
		strand = str(row['strand'])
		name = str(row['attributes'])
		start = start - 1

		if location in fasta_dic:
			sequence = fasta_dic[location]

			small = 0

			if strand == "+":
				UTR = sequence[start:end+N] #UTR assigned to 100nt after gene
				if end+N > len(sequence):
					small = 1
					UTR = sequence[end-3:]
			elif strand == "-":
				UTR = sequence[start-N:end]
				UTR = reverse_complement(UTR)
				if start-N < 0:
					small = 1
					UTR = sequence[:start+3]
					UTR = reverse_complement(UTR)
			else:
				missingDir += 1
				continue
			
			if small:
				badSmall += 1
			elif UTR[-N:].count('N') > N // 10:
				badN += 1
			elif UTR[-N-3:-N] in stops and UTR[0:3]=='ATG':
				good +=1
				UTRs.append([name,UTR[-N-3:-N],UTR[-N:]])
				stopList.append(UTR[-N-3:-N])
			else: 
				if UTR[-N-3:-N] not in stops:
					badStop += 1
				elif UTR[0:3]!='ATG':
					badStart += 1
	
	utrFrame = pd.DataFrame(UTRs, columns = ['Gene','Stop','UTR'])
	utrFrame.to_csv(utrFile,sep='\t')
	

	# Send info with small report
	TAG = stopList.count('TAG')
	TAA = stopList.count('TAA')
	TGA = stopList.count('TGA')
	total_stops = TGA+TAA+TAG
	
	if total_stops != good:
		print("Problem number of stops doesn't match genes: ", org)
	
	bad = badN + badSmall + badStart + badStop
	info = [org, total, good, bad, badN, badStop, badStart, badSmall, missingDir, TGA, TAA, TAG]
	
	return (utrFrame,info)
	
def reverse_complement(s):
	"""returns reverse complement of DNA"""
	
	seq = s[::-1]
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	bases = list(seq) 
	bases = [complement.get(base,"N") for base in bases] 
	return ''.join(bases)

def main():
	import os
	import pandas as pd

	info = pd.read_excel('GenomeInfoNew.xlsx')

	for number, organism in info.iterrows():
		if (os.path.exists('NewGenomes/'+organism['GFF File']) and
			os.path.exists('NewGenomes/'+organism['Assembly File'])):
			#print(organism['Organism'],': OK')
			continue
		else:
			print(organism['Organism'],': Not OK')
		#print(organism['Organism'],organism['GFF File'],organism['Assembly File'])
		
	from Bio import SeqIO
	import re
	import functionsUTR as utr
	utrSize = 50

	allInfo = []
	for number, organism in info.iterrows():
		org = organism['Organism']
		gffFile = 'NewGenomes/' + organism['GFF File']
		fasFile = 'NewGenomes/' + organism['Assembly File']
		marker = organism['GFF Marker']
		gffType = organism['Feature']
		utrFile = 'utrsNew/' + organism['UTR File']
		stops = organism['GeneticCode']
		
		fasta_dic = utr.parse_fasta(fasFile)
		gff = utr.process_gff(gffFile,marker,gffType)
		utrs, details = utr.get_UTRs(org,gff,fasta_dic,stops,utrSize,
						   utrFile)
		print(details)
		allInfo.append(details)
	return(0)
	

#main()