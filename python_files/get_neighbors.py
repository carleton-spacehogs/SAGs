import sys
import os

def get_genes_and_neighbors(cluster, pa_list):

	command = 'makeTabDelimitedRow.py "all_I_2.0_c_0.4_m_maxbit"' + ' "' + str(cluster) + \
	'" | db_getGenesInClusters.py > outfile.txt'
	run_program = os.popen(command)
	status = run_program.close()
	
	outfile = open('outfile.txt')
	
	lines = outfile.readlines()
	
	for line in lines:
		line = line.split()
		
		orf_name = line[2]
		neighbor_command = 'echo "' + orf_name + '" | db_getGeneNeighborhoods.py  >> ' + str(pa_list) + 'neighbors.txt'
		
		run_program = os.popen(neighbor_command)
		status = run_program.close()
	

	outfile.close()
	
def get_all_neighbors_from_pa(pa_list):
	cluster_list = []
	pa_file = open(pa_list + '_exclusive.txt')
	lines = pa_file.readlines()
	for line in lines:
		line = line.split()
		cluster_list.append(int(line[1]))
		
	for c in cluster_list:
		get_genes_and_neighbors(c, pa_list)
		
if __name__ == '__main__':
	get_all_neighbors_from_pa(sys.argv[1])
		