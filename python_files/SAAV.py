#Author: Michael Hoffert
#Date: 6/23/2017
#PI: Rika Anderson
#Python version 2.7.12

'''Makes SAAV object using data from anvi'o AA variability file
   Class instance variables correspond to columns in anvi'o output data
'''

class SAAV:
    
    def __init__(self, variability_list):
            
        self.entry_id = variability_list[0]
        self.unique_pos_identifier = variability_list[1]
        self.gene_length = int(variability_list[2])
        self.sample_id = variability_list[3]
        self.corresponding_gene_call = variability_list[4]
        self.codon_order_in_gene = variability_list[5]
        self.reference = variability_list[6]
        self.departure_from_reference = float(variability_list[7])
        self.coverage = int(variability_list[8])
        self.AA_list = variability_list[9:30]
        if variability_list[30][0] in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            #print 'FUCKING FOUND ONE'
            #print variability_list
            self.BLOSUM62 = 'NO VALUE'
            self.BLOSUM90 = 'NO VALUE'
            self.competing_aas = variability_list[30]
            self.consensus = variability_list[31]
            self.departure_from_consensus = variability_list[32]
            self.n2n1ratio = float(variability_list[33])
            self.contig_name = variability_list[34]
            
        else:
            self.BLOSUM62 = variability_list[30]
            self.BLOSUM90 = variability_list[31]
            self.competing_aas = variability_list[32]
            self.consensus = variability_list[33]
            self.departure_from_consensus = variability_list[34]
            self.n2n1ratio = float(variability_list[35])
            self.contig_name = variability_list[36]
            
    def update_contig_name(self, old_name, new_name):
        
        if old_name == self.contig_name:
            self.contig_name = new_name
            return True
        else:
            return False
