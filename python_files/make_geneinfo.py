import glob
import os


def make_geneinfo_files():
	path = r'./'


	for cluster in range(1, 1269):

		cline = 'makeTabDelimitedRow.py "all_I_2.0_c_0.4_m_maxbit"' + ' "' + str(cluster) + \
		'" | db_getGenesInClusters.py |	 db_getAlignmentBetweenGenes.py -g 3  -l > ' + str(cluster) + '_sequence_alignment'
		
		run_program = os.popen(cline)
		status = run_program.close()
		
		geneinfoline = 'makeTabDelimitedRow.py "all_I_2.0_c_0.4_m_maxbit"' + ' "' + str(cluster) + \
		'" | db_getClusterGeneInformation.py | getClusterFastas.py -n run_cluster_fastas'
		
		
		run_program = os.popen(geneinfoline)
		status = run_program.close()
	
		
def make_FastTree():
	
	for cluster in range(1, 1269):
		cline = 'cat ' + str(cluster) + '_sequence_alignment | FastTree_wrapper.py -m WAG - p fasttree > ' + str(cluster) + '_fasttree'
	
		run_program = os.popen(cline)
		status = run_program.close()
		
def fasta_to_phylip(filelist):
	for file in filelist:
		index = file.index('_')
		sqa_file = open(file)
		lines = sqa_file.readlines()
		rows = 0
		current = 0
		start = 0
		stop = 0
		startlist = []
		
		while current != len(lines):
			
			if lines[current][0] == '>':
				start = current
				current += 1
				startlist.append(start)
				
			else:
				current += 1
				
		startlist.append(len(lines))
		phylip = open(file[:index] + '_aligned.phi', 'w')
		writelist = []
		
		for i in range(len(startlist) - 1):
		
			name = lines[startlist[i]][1:-1] + '  '
			sequence = lines[startlist[i] + 1 : startlist[i + 1]]
			seq_str = ''
			for s in sequence:
				seq_str += s[:-1]
			rows += 1
			columns = len(seq_str)
			print(name, seq_str)

			
			writelist.append(name + seq_str + '\n')
		writelist.append(str(rows) + ' ' + str(columns) + '\n')
		phylip.write(writelist[-1])
		for s in writelist[:-1]:
			phylip.write(s)
		phylip.close()
		sqa_file.close()


			
def make_ctl_file(cluster_in):
	template = open('codeml.ctl','r')
	yn_00_template = open('yn00.ctl','r')
	writefile = open('codeml' + str(cluster_in) + '.ctl', 'w')
	yn00_writefile = open('yn00' + str(cluster_in) + '.ctl', 'w')
	template_lines = template.readlines()
	for line in range(3):
		start = template_lines[line].index('=')
		stop = template_lines[line].index('_')
		template_lines[line] = template_lines[line][: (start + 1)] + ' ' + str(cluster_in) + template_lines[line][stop:]
	#print(template_lines)
	
	for line in template_lines:
		writefile.write(line)
	template.close()
	writefile.close()
	
	yn00_lines = yn_00_template.readlines()
	
	for line in range(2):
		start = yn00_lines[line].index('=')
		stop = yn00_lines[line].index('_')
		yn00_lines[line] = yn00_lines[line][: (start + 1)] + ' ' + str(cluster_in) + yn00_lines[line][stop:]
		
	for line in yn00_lines:
		yn00_writefile.write(line)
	
	
	
	
	
	yn00_writefile.close()
			
		
def call_paml():
	for cluster in range(1, 1269):
		make_ctl_file(cluster)
		cline = 'codeml ' + 'codeml' + str(cluster) + '.ctl'
		run_program = os.popen(cline)
		status = run_program.close()
		
		cline = 'yn00 ' + 'yn00' + str(cluster) + '.ctl'
		run_program = os.popen(cline)
		status = run_program.close()
		

		
		
				
		
		
if __name__ == '__main__':

	#make_geneinfo_files()
	
	# filelist = []
# 	for i in range(1, 1839):
# 		filelist.append(str(i) + '_sequence_alignment')
# 	fasta_to_phylip(filelist)
#	make_FastTree()
	call_paml()
