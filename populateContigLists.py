#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Constructs lists of contigs and SNVs from input files"""
from Contig import *
from SNV import *
from ORF import *
from SAAV import *
from makeSequenceObjects import *
import os


# populates lists of SAAVs, SNVs, and ORFs in each contig
def populate_contig_lists(contigs, SNVs, SAAVs, ORFs):
     
    count = 0
    num_contigs = len(contigs)
    for contig in contigs:
        contig_SNV_list = []
        contig_ORF_list = []
        contig_SAAV_list = []
        
        for snv in SNVs:
            if snv.contig_name == contig.ID:
                contig_SNV_list.append(snv)
        
        for saav in SAAVs:
            if saav.contig_name == contig.ID:
                contig_SAAV_list.append(saav)
        
        
        for orf in ORFs:
            if orf.contig == contig.ID:
                contig_ORF_list.append(orf)
        
        contig.ORF_list = contig_ORF_list
        contig.SNV_list = contig_SNV_list
        contig.SAAV_list = contig_SAAV_list

    return contigs
        
if __name__ == '__main__':
    c = make_contigs('AD-155-C09-many_assemblers_simple.CISA.ctg.fa', 'C09_contig_names.tsv','C09.names_map')
    s = make_SNVs('C09_variability_new.txt')
    list1 = make_SNV_contig_names_consistent(s, c)
    o = make_ORFs('2703719243.gff','2703719243.ko.tab.txt')
    x = make_contig_ORF_and_SNV_lists(list1[1], list1[0], o)

        
        
        