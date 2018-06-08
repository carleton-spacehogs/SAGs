#module to construct gene pa figure
from fuzzywuzzy import fuzz
def read_file(tsv_file,search_terms):
	
	opd = dict.fromkeys(search_terms,0)
	opd.keys() 
	tsv_file = open(tsv_file)
	highest = 0
	lines = tsv_file.readlines()
	for term in search_terms:
		for line in lines:
			line = line.split('\t')
			#print(line[10], term)
			if fuzz.ratio(line[10],term) > 65:
				#print('a high ratio found!',fuzz.ratio(line[10],term))
				#print(line[10],term)
				opd[term] += 1
				
	return opd
				
	
	
	
	tsv_file.close()

def merge_and_write():
	opf = open('output_test.txt','w')
	dicts = []
	for SAG in ['2703719243','2703719244','2703719245','2703719246','2703719247']:
		dicts.append(read_file('./../SAG_data_files/'+SAG+'.ko.tab.txt',search_terms))
	for term in search_terms:
		writestr = term + '\t'
		for dict in dicts:
			#print(dict[term])
			writestr += str(dict[term]) + '\t'
		opf.write(writestr + '\n')
		
	opf.close()
	

search_terms = ['siderophore biosynthesis',
'TonB-dependent Fe acquisition',
'iron transport protein (FeoB_CN)',
'iron transport protein (FeoA)',
'siderophore biosynthesis (lucA_lucC)',
'ech hydrogenase (echABCE)',
'energy-converting hydrogenase A (ehaBCEGNOP)',
'energy-converting hydrogenase B (ehbABFIJKLNO)',
'coenzyme F420 hydrogenase (frhAGDG)',
'hydrogenase (hyaABC, hybO, hybC)',
'[NiFe] hydrogenase (hydA2A3B2B3)',
'methane/ammonia monooxygenase (pmoA-amoA)',
'methyl-coenzyme M reductase (mcrABCDG)',
'nitrous-oxide reductase (nosZ)',
'nitrite reductase (NO-forming)/hydroxylamine reductase (nirKS)',
'nitric oxide reductase (norB)',
'nitrite reductase (cytochrome c-552) - (nrfA)',
'cytochrome c oxidase (coxABCD,AC,ctaF)',
'cytochrome c oxidase cbb3-type (ccoNOPQ_NO)',
'carbon-monoxide dehydrogenase (cooCFS, acsA)',
'ATP-citrate lyase (aclAB)',
'fumarate reductase (frdAB) subunitAB',
'nitrogenase molybdenum-iron protein NifH, NifD, NifK',
'Periplasmic nitrate reductase (napAGHBFLD)',
'Nitrate reductase (narB)',
'Nitrite and sulfite reductase/Ferrodoxin nitrite reductase (nirA)',
'Cytochrome cd1 nitrite reductase (NO-forming) (nirs)',
'Sulfate permease',
'Sulfate adenylyltransferase  (cysD/cysN)',
'Sulfate adenylyltransferase (met3)',
'Adenylylsulfate kinase (apsK/cysC)',
'Phosphoadenosine phosphosulfate reductase/PAPS reductase',
'NAD(P)H: polysulfide oxidoreductase',
'Serine O-acetyltransferase (cysE)',
'Cysteine synthase (cysK)',
'Polysulfide reductase',
'Sulfur oxidation system (sox)',
'Sulfide:quinone oxidoreductase (sqr)',
'Sulfite:cytochrome c oxidoreductase (sorAB)']

if __name__ == '__main__':
	merge_and_write()