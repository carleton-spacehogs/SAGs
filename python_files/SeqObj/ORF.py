#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Class for ORFs"""

class ORF:
	
	def __init__(self, attribute_list):
		x = attribute_list[8].index(';')
		
		if not ((attribute_list[8][-1]) in '0123456789'):
			next_x = attribute_list[8].index(';', x + 1)
		else:
			next_x = -1
			
		self.ID = attribute_list[8][3:x]
		
		if not (next_x == -1):
			self.locus_tag = attribute_list[8][(x + 11) : next_x]
		else:
			self.locus_tag = attribute_list[8][(x + 11):]
		
		if ';' in self.locus_tag:
			x = self.locus_tag.index(';')
			self.locus_tag = self.locus_tag[:x]
		self.alias = None
		self.cluster = None
		self.contig = attribute_list[0]
		self.start = int(attribute_list[3])
		self.stop = int(attribute_list[4])
		self.strand = attribute_list[6]
		self.frame = attribute_list[7]
		self.data=attribute_list[8:]
		self.KO = None
		self.number_of_SNVs = 0
		self.number_of_SAAVs = 0
		self.SAAVs = []
		self.flag = True
		self.SNVs = []
		self.KO_A = None
		self.KO_B = None
		self.KO_C = None
		self.KO_data = None
		self.length = self.stop - self.start
		self.mean_n2n1 = 0
		self.filtered_SNVs = 0
		self.pa = 0
		
		
		
		
		