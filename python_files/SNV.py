#Author: Michael Hoffert
#Date: 6/15/2017
#PI: Rika Anderson
#Python version 2.7.12


""" This is an SNV class for holding data parsed from Anvi'o variability
profiles. It stores all data associated with the SNV profile."""

"""hoffertm@baross ~ $ anvi-profile --version

Anvi'o version ...............................: 2.3.1-master

Profile DB version ...........................: 20

Contigs DB version ...........................: 8

Pan DB version ...............................: 5

Samples information DB version ...............: 2

Genome data storage version ..................: 1

Auxiliary data storage version ...............: 3

Anvi'server users data storage version .......: 1"""


class SNV:
    
    def __init__(self, linelist): #26 entry list from anvio SNV profile
        self.entry_id = int(linelist[0])
        self.unique_pos_id = int(linelist[1])
        self.gene_length_id = int(linelist[2])
        self.sample_id = linelist[3]
        self.pos = int(linelist[4])
        self.pos_in_contig = int(linelist[5])
        self.corresponding_gene_call = int(linelist[6])
        self.in_partial_gene_call = int(linelist[7])
        self.in_complete_gene_call = int(linelist[8])
        self.base_pos_in_codon = int(linelist[9])
        self.codon_order_in_gene = int(linelist[10])
        self.coverage = int(linelist[11])
        self.cov_outlier_in_split = int(linelist[12])
        self.cov_outlier_in_contig = int(linelist[13])
        self.departure_from_reference = float(linelist[14])
        self.competing_nts = linelist[15]
        self.reference = linelist[16]
        self.A = int(linelist[17])
        self.T = int(linelist[18])
        self.C = int(linelist[19])
        self.G = int(linelist[20])
        self.N = int(linelist[21])
        self.consensus = linelist[22]
        self.departure_from_consensus = float(linelist[23])
        self.n2n1ratio = float(linelist[24])
        self.contig_name = linelist[25]
        self.SAG = self.sample_id[:3]
        self.metagenome = self.sample_id[4:9]
        
    def __str__(self):
        return 'SNV at pos %s from contig %s' % (self.pos, self.contig_name)
    
    def __repr__(self):
        return 'SNV at pos %s from contig %s' % (self.pos, self.contig_name)
    
    def update_contig_name(self, old_name, new_name):
        
        if old_name == self.contig_name:
            self.contig_name = new_name
            return True
        else:
            return False
    
        
        

        
        

        
        
        