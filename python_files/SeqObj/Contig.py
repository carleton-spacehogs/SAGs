#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Class for Contigs and Contig associated data"""

class Contig:
    
    def __init__(self, ID, sequence, SAG):
        self.first_ID = ID
        self.ID = ID
        self.sequence = sequence
        self.SAG = SAG
        self.SNV_list = None
        self.ORF_list = None
        self.SAAV_list = None
        self.length = len(self.sequence)
        self.num_of_SNVs = 0
        
    def update_ID(self, old_ID, new_ID):
        
        if self.ID == old_ID:
            self.ID = new_ID
            return True
        else:
            return False
    
    def __repr__(self):
        return self.ID
    
    def __str__(self):
        return self.ID
    
    
        
        
