#file for making exclusive gene lists and info files
import os
import itertools

SAG_names =	 ['Methanothermococcus_archaeon_SAG_AD_155_C09',
			  'Methanothermococcus_archaeon_SAG_AD_155_K20',
			  'Methanothermococcus_archaeon_SAG_AD_155_M21',
			  'Methanothermococcus_archaeon_SAG_AD_155_N22',
			  'Methanothermococcus_archaeon_SAG_AD_155_E23']

def make_name_file(pa_str):
	str = ''
	for i in range(len(pa_str)):
		if pa_str[i] == '1':
			str += SAG_names[i] + '\n'
	file = open(pa_str + '_names.txt','w')
	file.write(str)
	file.close()
	
def run_names_to_exclusive(pa_str):
	cline = "cat " + pa_str + "_names.txt" + " | db_findClustersByOrganismList.py -a -s all_I_2.0_c_0.4_m_maxbit > " + pa_str + "_exclusive.txt"
	
	run_program = os.popen(cline)
	status = run_program.close()
	
def make_info_file(pa_str):
	exclusive = open(pa_str + '_exclusive.txt')
	for line in exclusive.readlines():
		line = line.split()
		cluster = line[1]
		
		cline = 'makeTabDelimitedRow.py "all_I_2.0_c_0.4_m_maxbit" ' + str(cluster) + ' | db_getClusterGeneInformation.py >> ' + pa_str + 'exclusive_info.txt'
	
		run_program = os.popen(cline)
		status = run_program.close()

from fuzzywuzzy import fuzz
def make_tsv_from_file(pa_str):
	file = open(pa_str + 'exclusive_info.txt')
	
	rl = {}
	lines = file.readlines()
	nlines = []
	for line in lines:
		nline = line.split('\t')
		nline = nline[9][:nline[9].index('_Ga')]
		nlines.append(nline)
		print(nline)
	rl = {}
	for nline in nlines:
		if not (nline in rl.keys()):
			rl[nline] = 1
		for line in lines:
			if nline in line:
				
				rl[nline] += 1
				break
	
	

	file.close()
	out = open(pa_str + '.tsv', 'w')
	out.truncate()
	out.write('ORF annotation' + '\t' + 'count' + '\n')
	for key in rl.keys():
		out.write(key+'\t'+str(rl[key])+'\n')
	out.close()
		
	

if __name__ == '__main__':
	for pa_str in ["".join(seq) for seq in itertools.product("01", repeat=5)][1:]:
		#make_name_file(pa_str)
		#run_names_to_exclusive(pa_str)
		#make_info_file(pa_str)
		
		try:
			make_tsv_from_file(pa_str)
		except FileNotFoundError:
			print("sorry. shit's fucked")
		