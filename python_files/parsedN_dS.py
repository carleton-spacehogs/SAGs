#function for parsing dN/dS outputs
import pandas as pd
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
def parse_yn00_output():

	NUM_COLORS = 32 
	cm = plt.get_cmap('gist_rainbow')
	colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
	
	dict1 = {'cluster' : [],
			 'colors'  : [],
			 'dN' : [],
			 'dS' : []}

	for i in range(1269):
		filename = str(i) + + '_yn_output'
		dN = 0
		dS = 0
		cluster = int(filename[:filename.index('_')])
		output = open(filename)
		lines = output.readlines()
	
		for i in range(len(lines)):
			lines[i] = lines[i].split()
		
			if len(lines[i]) > 0:
				if lines[i][0] == '(B)':
					start = i
				if lines[i][0] == '(C)':
					stop = i
		for i in range((start + 8), (stop - 2)):

			dN += float(lines[i][7])
			dS += float(lines[i][10])
		dN, dS = dN / ((stop - 2) - (start + 8)), dS / ((stop - 2) - (start + 8))
		if dS != 0:
			dNdS = dN / dS
		else:
			dNdS = 'NAN'
		pa = get_cluster_pa(cluster)
		
		dict1['cluster'].append(cluster)
		dict1['dN'].append(dN)
		dict1['dS'].append(dS)
		dict1['colors'].append(colors[int(pa, 2)])
			
		output.close()
	return dict1
		

def get_cluster_pa(cluster):
	cluster_pa_file = open('cluster_output_noArea.txt')
	cluster_lines = cluster_pa_file.readlines()
	
	line = cluster_lines[cluster]
	line = line.split()
	
	rts = line[-5] + line[-4] + line[-3] + line[-2] + line[-1]
	rts = rts.replace(',','')[1:-1]
	
	return rts
		


	cluster_pa_file.close()
if __name__ == '__main__':
	mydict = parse_yn_output()
	df = pd.DataFrame.from_dict(data=graph_dict)
	df.plot.scatter(x='dN',y='dS',c='colors',alpha=0.75)
	
	plt.savefig('dN_dS_scatter.png')
	
	
	
	
	
	
	