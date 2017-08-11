#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Constructs lists of contigs and SNVs from input files"""
from Contig import *
from SNV import *
from SAAV import *
from ORF import *
import os

#example SAG fasta file name:
#AD-155-C09-many_assemblers_simple.CISA.ctg.fa
#example contig name .tsv file:
#C09_contig_names.tsv
#example names_map file:
#C09.names_map
#
#Note: The file names MUST be formatted in the above format

def make_contigs(SAG_fasta_file, contig_name_tsv, names_map, SAG_name):
	
	# uses SAG fasta file to load contig sequences
	# and mapped name/tsv files to update names
	# for contigs that were simplified using anvi-script
	
	SAG_fa_file = open(SAG_fasta_file)
	contig_names = open(contig_name_tsv)
	names_map_file = open(names_map)

	print('filenames checked out')

	contig_list = [] #output list for contigs

	SAG_fa_lines = SAG_fa_file.readlines()
	contig_name_lines = contig_names.readlines()

	for i in range(0, len(SAG_fa_lines), 2): #make contigs from fasta file

		ID = SAG_fa_lines[i][1 : (len(SAG_fa_lines[i]) - 1)]
		sequence = SAG_fa_lines[i + 1]
		SAG = SAG_name
		contig = Contig(ID, sequence, SAG)
		contig_list.append(contig)

	#update contig names from simplified versions
	#anvi-script -> SAG.fasta names
	for line in contig_name_lines:
		line = line.split()
		old_ID = line[0]
		new_ID = line[1]
		for i in range(len(contig_list)):
			if contig_list[i].update_ID(old_ID, new_ID):
				break

	#SAG.fasta names -> JGI contig names
	names_map_lines = names_map_file.readlines()
	for line in names_map_lines:
		line = line.split()
		old_ID = line[0]
		new_ID = line[1]
		for i in range(len(contig_list)):
			if contig_list[i].update_ID(old_ID, new_ID):
				break
	
	SAG_fa_file.close()
	contig_names.close()
	names_map_file.close()
	
	return contig_list


#makes list of SNVs from variability profile
def make_SNVs(variability_file):
	
	SNV_list = []
	
	variability_file = open(variability_file)
	variability_lines = variability_file.readlines()
	
	for line in variability_lines:
		line = line.split()
		if not (line[0] == 'entry_id'):
			
			snv = SNV(line)
			SNV_list.append(snv)
			
	return SNV_list

#makes list of SAAVs from variability profile
def make_SAAVs(variability_file):
	
	SAAV_list = []
	
	variability_file = open(variability_file)
	variability_lines = variability_file.readlines()
	
	for line in variability_lines:
		line = line.split()
		if not (line[0] == 'entry_id'):

			saav = SAAV(line)
			SAAV_list.append(saav)
	
	variability_file.close()
	return SAAV_list  


# makes ORF file from JGI gff and ko file
def make_ORFs(gff_file, ko_file, alias_file, cluster_file):
	
	ORF_list = []
	
	ko_file = open(ko_file)
	gff_file = open(gff_file)
	alias_file = open(alias_file)
	c_file = open(cluster_file)
	gff_lines = gff_file.readlines()
	ko_lines = ko_file.readlines()
	alias_lines = alias_file.readlines()
	c_lines = c_file.readlines()
	
	for line in gff_lines:
		line = line.split()
		if len(line) > 3:
			orf = ORF(line)
			ORF_list.append(orf)
	
	for orf in ORF_list:
		for line in ko_lines:
			line = line.split()
			if line[0] == orf.ID:
				orf.KO = line[9]
		
		for line in alias_lines:
			line = line.split()
			if orf.locus_tag == line[1]:
				
				orf.alias = line[0]
				
		for line in c_lines:
			line = line.split()
			if orf.alias == line[2]:
				orf.cluster = line[1]
	
	
	ko_file.close()
	gff_file.close()
	alias_file.close()
	c_file.close()
	
	return ORF_list	   

#iterates through SNVs and SAAVs, makes parent contig names identical
def make_names_consistent( Contigs, SNVs=None, SAAVs=None):
	
	for contig in Contigs:
		
		if SNVs != None:
			for snv in SNVs:
				if snv.update_contig_name(contig.first_ID, contig.ID):
					pass		 
		
		if SAAVs != None:
			for saav in SAAVs:
				if saav.update_contig_name(contig.first_ID, contig.ID):
					pass
						
	return [Contigs, SNVs, SAAVs]

#simple function for retrieving contig from contig list
def get_contig(ID, contigs):
	for contig in contigs:
		if ID == contig.ID:
			return contig
		
#function to determine if an ORF contains and SNV
#and adds SNV to ORF's list of SNVs
def populate_ORF_list(ORFs, SNVs, contigs):
	
	for ORF in ORFs:
		parent_contig_name = ORF.contig
		parent_contig = get_contig(parent_contig_name, contigs)
		potential_SNVs = parent_contig.SNV_list
		potential_SAAVs = parent_contig.SAAV_list
		for SNV in potential_SNVs:
			if SNV.pos_in_contig >= ORF.start and SNV.pos_in_contig <= ORF.stop:
				ORF.number_of_SNVs += 1
				ORF.SNVs.append(SNV)
				
		if potential_SAAVs != None:
			for SAAV in potential_SAAVs:
				if (SAAV.gene_length + 2 >= ORF.length) and (SAAV.gene_length - 2 <= ORF.length):
					ORF.number_of_SAAVs += 1
					ORF.SAAVs.append(SAAV)
				
	return ORFs
				

	
		
if __name__ == '__main__':
	
	c = make_contigs('AD-155-C09-many_assemblers_simple.CISA.ctg.fa', 'C09_contig_names.tsv','C09.names_map')
	s = make_SNVs('C09_variability_new.txt')
	make_SNV_contig_names_consistent(s, c)
	
	
	
	
	
	
	
	
	
	
	