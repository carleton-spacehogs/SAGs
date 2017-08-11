#mapping script
import os


import pandas as pd
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np


SAG_LIST = ['C09','E23','M21','N22','K20']


def map_SAGs():
	path = r'./'
	
	for ref_SAG in SAG_LIST:
		df_list = []
		colors = []
		base_list = []
		coverage = []
		for read_SAG in SAG_LIST:
			if ref_SAG != read_SAG:
				reads = './SAG_reads/SAG_reads/AD-155-'+read_SAG+'/merged/MidCaymanRise_DCO_Methanothermococcus_SAGs_'+read_SAG+'_merged_MERGED_filtered.fa'
				ref = 	'../hoffertm/summer_2017/mappings/'+ref_SAG+'_mappings/AD-155-'+ref_SAG+'-many_assemblers_simple.CISA.ctg.fa'
		
# 				command1 = 'bowtie2-build '+ref+' '+ref[46:-3]+'.btindex'
# 				run_program = os.popen(command1)
# 				status = run_program.close()
		
				# command2 = 'bowtie2 -x '+ ref[46:-3] + '.btindex' + ' -p 80 -f -U '+reads+' -S '+read_SAG + '_to_' + ref_SAG + '_mapped.sam'
# 				run_program = os.popen(command2)
# 				status = run_program.close()
# 				
# 				command3 = 'samtools view -bS ' + read_SAG + '_to_' + ref_SAG + '_mapped.sam' + ' > ' + read_SAG + '_to_' + ref_SAG + '_mapped.bam'
# 				run_program = os.popen(command3)
# 				status = run_program.close()
# 				
# 				command4 = 'samtools sort -@ 80 ' + read_SAG + '_to_' + ref_SAG + '_mapped.bam' + ' -o ' + read_SAG + '_to_' + ref_SAG + '_sorted' + '_mapped.bam'
# 				run_program = os.popen(command4)
# 				status = run_program.close()
# 				
# 				command5 = 'samtools faidx ' + ref
# 				run_program = os.popen(command5)
# 				status = run_program.close()
# 				
# 				command6 = 'samtools index ' + read_SAG + '_to_' + ref_SAG + '_sorted' + '_mapped.bam'
# 				run_program = os.popen(command6)
# 				status = run_program.close()
# 				
# 				command7 = 'bedtools genomecov -d -ibam ' + read_SAG + '_to_' + ref_SAG + '_sorted' + '_mapped.bam -g ' + 'AD-155-' + ref_SAG + '_genome.txt > ' + read_SAG + '_to_' + ref_SAG + '.genomecov ' + '&'
# 				run_program = os.popen(command7)
#  				status = run_program.close()

				x = make_genomecov_graph(read_SAG + '_to_' + ref_SAG + '.genomecov')
				df, color, dict = x[0], x[1], x[2]
				base_list += dict['base']
				coverage += dict['coverage']
				
				colors += color
				df_list.append(df)
		dN_dS_heatmap(base_list, coverage, ref_SAG)
		df = pd.concat(df_list)
		#ax = df.plot.scatter(x='base',y='coverage',color=colors,s=5, alpha=0.005)
	
		#plt.savefig(ref_SAG + '_scatter_edit_sums.png')
		
def dN_dS_heatmap(base_list, coverage_list, SAG):
	
	x = base_list
	y = coverage_list
	bins_x = len(base_list) / 1000
	
	heatmap, xedges, yedges = np.histogram2d(x, y, bins=(bins_x,1))
	extent = [0, 1000, 0, 100]

	plt.clf()
	plt.imshow(heatmap.T, extent=extent, origin='lower')
	plt.savefig(SAG + '_coverage_heatmap.png')
	plt.show()
					
def make_genomecov_graph(filename):
	
	mydict = {'base' : [],
			  'coverage' : [],
			  'colors' : []}
	
	file = open(filename)	
	lines = file.readlines()
	count = 1
	for line in lines:
		line = line.split()
		
		
		if int(line[2]) == 0:
			mydict['coverage'].append(0)
			mydict['base'].append(count)
			mydict['colors'].append('r')
		else:
			for i in range(2,100):
				mydict['base'].append(count)
				mydict['coverage'].append(i)
				mydict['colors'].append('b')
			
		
		count += 1
		if count % 10000 == 0:
			print(count)
		
	
	
	df = pd.DataFrame.from_dict(data=mydict)
	return [df, mydict['colors'], mydict]
	# ax = df.plot.scatter(x='base',y='coverage',color=mydict['colors'],s=5)
# 	
# 	plt.savefig(filename + '_scatter_edit.png')
	
	#plt.show()
	
	
	
	
	
	
	
	
			
				
if __name__ == '__main__':
	map_SAGs()