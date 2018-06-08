# makes a genbank file from a .gff file and a .fasta file and a .faa file. change output file name below.
# usage: make_genbank_file.py [fasta file] [gff file] [faa file]

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

#import JGI files
fasta_filehandle = sys.argv[1]
gff_filehandle = sys.argv[2]
aa_filehandle = sys.argv[3]

#make a list of translated proteins
proteinlist = []
protein_file = open(aa_filehandle).read()
protein_fastas = protein_file.split('>')
for protein_fasta in protein_fastas:
	proteinlist.append(protein_fasta)
	
#open and parse gff file
gfffile = open(gff_filehandle).read()
gff_lines = gfffile.split('\n')

#open and parse faa file and make a dictionary
Dfaa = {}
faafile = open(aa_filehandle).read()
faas = faafile.split('>')
for faa in faas:
	faa_lines = faa.split('\n')
	firstline = faa_lines[0]
	firstline_split = firstline.split(' ')
	faa_id = firstline_split[0]
	faa_seq = faa.strip(firstline)
	faa_seq = faa_seq.replace('\n', '')
	#print faa_seq
	Dfaa[faa_id] = faa_seq


	
#make an outfile-- **EDIT THIS, I WAS TOO LAZY TO USE INPUTS**
outfile = open('AD-155-C09.gbk', 'w')

#parse by fastas
fastafile = open(fasta_filehandle).read()
fastas = fastafile.split('>')
for fasta in fastas[1:]:
	
	lines = fasta.split('\n')
	titleline = lines[0]
	title_split = titleline.split(' ')
	title = title_split[0]
	sequence = fasta.strip(titleline)
	#print sequence
	sequence = sequence.replace('\n', '') #remove carriage returns from sequence
	#print fasta
	
	#create a sequence
	sequence_string = sequence
	sequence_object = Seq(sequence_string, IUPAC.unambiguous_dna)

	#create a record
	record = SeqRecord(sequence_object, id='unknown', name=title, description= title)
	
	#add a source feature for the record
	source_end = len(sequence)
	feature = SeqFeature(FeatureLocation(start=0, end=source_end), type='source', qualifiers={'db_xref':'taxon:155862', 'organism':'Methanothermococcus archaeon SAG AD-155-C09'})
	record.features.append(feature)

	# parse gff file second time, where each line is a feature
	for gff_line in gff_lines[1:-1]:
		D = {} #this is the qualifiers dictionary
		cols = gff_line.split('\t')
		contig = cols[0]
		feat = cols[2]
		start_coord = int(cols[3])-1
		stop_coord = int(cols[4])
		str = cols[6]
		if str == '+':
			gffstrand = 1
		elif str == '-':
			gffstrand = -1
		else:
			 gffstrand == 0
		
		annote = cols[8]
		annote_split = annote.split(';')
		id = annote_split[0].strip('ID=')
		loc_tag = annote_split[1].strip('locus_tag=')
		
		#assign keys in qualifiers dictionary
		D['db_xref'] = 'taxon:155862'
		D['locus_tag'] = loc_tag
		
		if feat == 'CDS':
			product = annote_split[2].strip('product=')
			D['product'] = product
			D['translation'] = Dfaa[id]

		
		#if the line matches the contig we're parsing, write it out	
		if contig == title:
			#print contig
	
			#now add annotation
			feature = SeqFeature(FeatureLocation(start=start_coord, end=stop_coord), type=feat, strand=gffstrand, qualifiers=D)
			record.features.append(feature)

	#save as genbank file
	SeqIO.write(record, outfile, 'genbank')

outfile.close()