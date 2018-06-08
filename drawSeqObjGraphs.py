#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

'''Constructs lists of contigs and SNVs from input files'''
from Contig import *
from SNV import *
from SAAV import *
from ORF import *

import makeSequenceObjects
import organizeKO
import populateContigLists
import mergeORFpaData
#import contrastText

import pandas as pd
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn import mixture
from scipy.spatial import ConvexHull
import os
import sys
import time
import pickle


plt.rcParams["font.family"] = "serif"

RUNTRUE = True

SAG_lengths = { #lengths of SAGS in bp
	'E23':1097435,
	'C09':1753024,
	'N22':1028235,
	'M21':1105061,
	'K20':628966
	}
	
metagenome_reads = { #number of reads in each metagenome, bp
	'FS841':53842644,
	'FS842':64300743,
	'FS844':82138820,
	'FS848':77333359,
	'FS849':26061468,
	'FS851':75912423,
	'FS852':37753638,
	'FS854':140664299,
	'FS856':104123722,
	'FS866':24088594,
	'FS872':40397883,
	'FS874':65834175,
	'FS877':152063932,
	'FS879':24723930,
	'FS881':167339253
				   }
sample_sites = [
			'FS841',
			'FS842',
			'FS844',
			'FS848',
			'FS849',
			'FS851',
			'FS852',
			'FS854',
			'FS856',
			'FS866',
			'FS872',
			'FS874',
			'FS877',
			'FS879',
			'FS881']

temps = {
		'FS841':114,
		'FS842':86,
		'FS844':50,
		'FS848':47,
		'FS849':109,
		'FS851':106,
		'FS852':44,
		'FS854':18,
		'FS856':108,
		'FS866':130,
		'FS872':30,
		'FS874':140,
		'FS877':131,
		'FS879':29,
		'FS881':114}
				   
mapping_coverage = { #coverage of each metagenome mapping
	'C09_FS841':0.85074292,'K20_FS841':0.521733092,'N22_FS841': 0.776108175,'M21_FS841':0.005984661,'E23_FS841':0.578942662,
	'C09_FS842':2.687491353,'K20_FS842':2.687537325,'N22_FS842':13.70974501,'M21_FS842':1.714508436,'E23_FS842':4.885482009,
	'C09_FS844':0.087377417,'K20_FS844':0.034499176,'N22_FS844':0.026152219,'M21_FS844':0.108001639,'E23_FS844':0.024126383,
	'C09_FS848':108.2655521,'K20_FS848':127.1934995,'N22_FS848':34.96652093,'M21_FS848':116.9982544,'E23_FS848':128.0653937,
	'C09_FS849':0.266736873,'K20_FS849':0.359053378,'N22_FS849':0.626042108,'M21_FS849':0.050728756,'E23_FS849':0.23392349,
	'C09_FS851':1.705616837,'K20_FS851':3.396945504,'N22_FS851':4.520771031,'M21_FS851':0.067383183,'E23_FS851':2.324827237,
	'C09_FS852':3.966187332,'K20_FS852':8.987104562,'N22_FS852':11.8071716,'M21_FS852':0.080472645,'E23_FS852':5.298544029,
	'C09_FS854':5.1191866,'K20_FS854':5.79856688,'N22_FS854':4.743035247,'M21_FS854':5.729131702,'E23_FS854':1.072089439,
	'C09_FS856':1.524237066,'K20_FS856':1.634765377,'N22_FS856':1.000254865,'M21_FS856':1.783870947,'E23_FS856':0.210610362,
	'C09_FS866':14.39162359,'K20_FS866':19.1174054,'N22_FS866':69.23744747,'M21_FS866':12.33293687,'E23_FS866':288.0378849,
	'C09_FS872':0.11233089,'K20_FS872':0.009309417,'N22_FS872':0.013612907,'M21_FS872':0.010954871,'E23_FS872':0.020586681,
	'C09_FS874':0.197947084,'K20_FS874':0.010388171,'N22_FS874':0.071649855,'M21_FS874':0.01218419,'E23_FS874':0.060701726,
	'C09_FS877':79.29027457,'K20_FS877':95.06176238,'N22_FS877':18.12472016,'M21_FS877':91.7939652,'E23_FS877':51.7186239,
	'C09_FS879':3.342027416,'K20_FS879':4.233221034,'N22_FS879':7.233275936,'M21_FS879':3.452026452,'E23_FS879':1.957391507,
	'C09_FS881':13.2611894127,'K20_FS881':41.8011974258,'N22_FS881':42.26812025,'M21_FS881':8.45602434432,'E23_FS881':434.018177602
	}


#function to generate Contig, SNV, SAAV, and ORF objects
def make_sequence_objects(SAG):
	
	print('Running for SAG',SAG)
	
	#path = '/Users/hoffertm/Documents/Summer_2017/KO_SNV/'
	path = '/Users/hoffertm/Documents/Summer_2017/Git_DC/SAGs/SAG_data_files/'
	
	SAG_fasta = 'AD-155-' + SAG + '-many_assemblers_simple.CISA.ctg.fa'
	names_tsv = SAG + '_contig_names.tsv'
	names_map = SAG + '.names_map'
	SNV_var = SAG + '_variability_new.txt'
	AA_var = SAG + '_variability_AA.txt'
	clusters = 'all_I_2.0_c_0.4_m_maxbit'
	
	JGI_dict = {'E23': ['2703719244.gff', '2703719244.ko.tab.txt'],
				'C09': ['2703719243.gff', '2703719243.ko.tab.txt'],
				'M21': ['2703719246.gff', '2703719246.ko.tab.txt'],
				'K20': ['2703719245.gff', '2703719245.ko.tab.txt'],
				'N22': ['2703719247.gff', '2703719247.ko.tab.txt']
			   }
	
	Contigs = makeSequenceObjects.make_contigs(path + SAG_fasta, \
											   path + names_tsv, \
											   path + names_map, \
											  SAG)
	SNVs = makeSequenceObjects.make_SNVs(path + SNV_var)
	SAAVs = makeSequenceObjects.make_SAAVs(path + AA_var)
	ORFs = makeSequenceObjects.make_ORFs(path + JGI_dict[SAG][0], path + JGI_dict[SAG][1],path+'aliases',path+clusters)
	opl = makeSequenceObjects.make_names_consistent(Contigs, SNVs, SAAVs)
	Contigs = populateContigLists.populate_contig_lists(opl[0], opl[1], opl[2], ORFs)
	KO_dict = organizeKO.parse_KO(path + 'ko00001.keg')
	ORFs = makeSequenceObjects.populate_ORF_list(ORFs, SNVs, Contigs)
	ORFs = organizeKO.assign_KO_to_ORFs(ORFs, KO_dict)
	
	return [SAG, Contigs, SNVs, SAAVs, ORFs]
	

def pickle_SeqObj(SeqObjList):
	
	SAG = SeqObjList[0]
	opf = open(SAG + '_SeqObj.pickle','wb')
	pickle.dump(SeqObjList, opf)
	print('The SAG data was successfully pickled.')
	opf.close()


def unpickle_SeqObj(SAG):

	openFile = open(SAG + '_SeqObj.pickle', 'rb')
	SeqObjList = pickle.load(openFile)
	print('The SAG data was successfully unpickled.')
	openFile.close()
	
	return SeqObjList
		
	
def make_DataFrame(graphtype, SeqObjList):
	
	SAG = SeqObjList[0]
	if SAG == 'C09':
		color = 'green'
	elif SAG == 'E23':
		color = 'blue'
	elif SAG == 'M21':
		color = 'orange'
	elif SAG == 'K20':
		color = 'red'
	else:
		color = 'purple'
	
	if graphtype == 'bubble':
		
		graph_dict = {'names' : [],
					  'n2n1' : [],
					  'SNVs' : [],
					  'SAAVs' : [],
					  'coverage' : []
					 }
		
		for c in SeqObjList[1]:

			if c.SAG == SAG:
				SNVs = c.SNV_list
				for snv in SNVs:
					if not (snv.sample_id in graph_dict['names']):
						graph_dict['names'].append(snv.sample_id)
						graph_dict['coverage'].append(mapping_coverage[snv.sample_id[:9]] * 10.0)
						graph_dict['n2n1'].append([snv.n2n1ratio, 1])
						graph_dict['SNVs'].append(1)
					else:
						sample_index = graph_dict['names'].index(snv.sample_id)
						graph_dict['n2n1'][sample_index][0] += snv.n2n1ratio
						graph_dict['n2n1'][sample_index][1] += 1
						graph_dict['SNVs'][sample_index] += 1
		
		
		length = SAG_lengths[SAG] / 1000.0
		for i in range(len(graph_dict['n2n1'])):
			
			graph_dict['n2n1'][i] = graph_dict['n2n1'][i][0] / graph_dict['n2n1'][i][1]
			sample_name = graph_dict['names'][i][4:9]
			graph_dict['SNVs'][i] = graph_dict['SNVs'][i] / length
			   
		NUM_COLORS = 15
		pc_num_colors = 91
		vd_num_colors = 112
		
		sizes = graph_dict['coverage']
		
		cm = plt.get_cmap('gist_rainbow')
		von_damm = plt.get_cmap('autumn')
		piccard = plt.get_cmap('winter')
		
		pc_colors = [piccard(1.*i/pc_num_colors) for i in range(pc_num_colors)]
		vd_colors = [von_damm(1.*i/vd_num_colors) for i in range(vd_num_colors)]
		pc_colors.reverse()
		vd_colors.reverse()
		normal_colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
		
		colors = []
		for item in graph_dict['names']:
			if '5' in item:
				temp = temps[item[4:9]]
				print('temp',temp-18)
				c = pc_colors[temp - 18]
				colors.append(c)
			else:
				temp = temps[item[4:9]]
				print('temp',temp-29)
				c = vd_colors[temp - 29]
				colors.append(c)
		
		
		df = pd.DataFrame.from_dict(data=graph_dict)
		
		return [df, colors, sizes]

		
		
	elif graphtype == 'histo':
		graph_dict = {'n2n1': []}
		
		for SNV in SeqObjList[2]:
			if SNV.n2n1ratio > 0.0:
				graph_dict['n2n1'].append(SNV.n2n1ratio)
		df = pd.DataFrame.from_dict(data=graph_dict)
		
	elif graphtype == 'scatter_ORFs_SAAV':
	
		color_dict = {'A' : [],
						  'B' : [],
						  'C' : []
						  }
		for c in SeqObjList[1]: 
			if c.SAG == SAG:
				SAAVs = c.SAAV_list
				ORFs = c.ORF_list
			
				for orf in ORFs:
				
					A = orf.KO_A
					B = orf.KO_B
					C = orf.KO_C
				
					if not (A in color_dict['A']):
						color_dict['A'].append(A)
				
					if not (A in color_dict['B']):
						color_dict['B'].append(B)
				
					if not (A in color_dict['C']):
						color_dict['C'].append(C)
				
					NUM_COLORS_A = len(color_dict['A'])
					NUM_COLORS_B = len(color_dict['B'])
					NUM_COLORS_C = len(color_dict['C'])
					cm = plt.get_cmap('gist_ncar')
					colors_A = [cm(1.*i/NUM_COLORS_A) for i in range(NUM_COLORS_A)]
					colors_B = [cm(1.*i/NUM_COLORS_B) for i in range(NUM_COLORS_B)]
					colors_C = [cm(1.*i/NUM_COLORS_C) for i in range(NUM_COLORS_C)]
				
					color_dict['colors_A'] = colors_A
					color_dict['colors_B'] = colors_B
					color_dict['colors_C'] = colors_C
	
		graph_dict = {'names' : [],
					  'n2n1' : [],
					  'SAAVs' : [],
					  'length' : [],
					  'colors' : []
					 }
		print('CAUTION: \nA robust method for determining ORF/SAAV placement has not been implemented')
		print('A grain of salt is required for these results')			 
		
		for c in SeqObjList[1]:
			ORFs = c.ORF_list
			if c.SAG == SAG:
				for orf in ORFs:
					A = orf.KO_A
					B = orf.KO_B
					C = orf.KO_C
					
					mean_n2n1 = 0
					count = 0
					if len(orf.SAAVs) > 0 and orf.KO_A != None:
						
						for saav in orf.SAAVs:
							
							if saav.departure_from_reference > 0.05 and \
							(mapping_coverage[saav.sample_id[:9]] >= 10.0):
								mean_n2n1 += saav.n2n1ratio
								count += 1
						
						if count > 0:
							mean_n2n1 = mean_n2n1 / float(count)
							graph_dict['names'].append(orf.ID)
							graph_dict['n2n1'].append(mean_n2n1)
							graph_dict['SAAVs'].append(count)
							graph_dict['length'].append(orf.length)
							
#							ind = color_dict['A'].index(A)
#							color = color_dict['colors_A'][ind]
#							graph_dict['colors'].append(color)
						
#							ind = color_dict['B'].index(B)
#							color = color_dict['colors_B'][ind]
#							graph_dict['colors'].append(color)
#						
							ind = color_dict['C'].index(C)
							color = color_dict['colors_C'][ind]
							graph_dict['colors'].append(color)
					

		for i in range(len(graph_dict['n2n1'])):
			graph_dict['SAAVs'][i] = graph_dict['SAAVs'][i] / float(graph_dict['length'][i])
		df = pd.DataFrame.from_dict(data=graph_dict)
				
	
		
	elif graphtype == 'scatter_ORFs_SNV':
	
		graph_dict = {'names' : [],
					  'n2n1' : [],
					  'SNVs' : [],
					  'length' : [],
					  'colors' : [],
					 }
					  
		#NUM_COLORS = 7
		#colors = ['r','b','g','y','m','c','k']
		
		
		
		for c in SeqObjList[1]:
			color_dict = {'A' : [],
						  'B' : [],
						  'C' : []
						  }
			
			if c.SAG == SAG:
				SNVs = c.SNV_list
				ORFs = c.ORF_list
				
				for orf in ORFs:
					
					A = orf.KO_A
					B = orf.KO_B
					C = orf.KO_C
					
					if not (A in color_dict['A']):
						color_dict['A'].append(A)
					
					if not (A in color_dict['B']):
						color_dict['B'].append(B)
					
					if not (A in color_dict['C']):
						color_dict['C'].append(C)
					
					NUM_COLORS_A = len(color_dict['A'])
					NUM_COLORS_B = len(color_dict['B'])
					NUM_COLORS_C = len(color_dict['C'])
					cm = plt.get_cmap('gist_ncar')
					colors_A = [cm(1.*i/NUM_COLORS_A) for i in range(NUM_COLORS_A)]
					colors_B = [cm(1.*i/NUM_COLORS_B) for i in range(NUM_COLORS_B)]
					colors_C = [cm(1.*i/NUM_COLORS_C) for i in range(NUM_COLORS_C)]
					
					color_dict['colors_A'] = colors_A
					color_dict['colors_B'] = colors_B
					color_dict['colors_C'] = colors_C
					
				for orf in ORFs:
					A = orf.KO_A
					B = orf.KO_B
					C = orf.KO_C
					
					mean_n2n1 = 0
					count = 0
					
					if len(orf.SNVs) > 0 and orf.KO_A != None:
						
						for snv in orf.SNVs:
							mean_n2n1 += snv.n2n1ratio
							count += 1
						
						mean_n2n1 = mean_n2n1 / float(count)
						graph_dict['names'].append(orf.ID)
						graph_dict['n2n1'].append(mean_n2n1)
						graph_dict['SNVs'].append(count)
						graph_dict['length'].append(orf.length)
#						ind = color_dict['A'].index(A)
#						color = color_dict['colors_A'][ind]
#						graph_dict['colors'].append(color)
						
#						ind = color_dict['B'].index(B)
#						color = color_dict['colors_B'][ind]
#						graph_dict['colors'].append(color)
						
						ind = color_dict['C'].index(C)
						color = color_dict['colors_C'][ind]
						graph_dict['colors'].append(color)
					
						
						
		

		for i in range(len(graph_dict['n2n1'])):
			graph_dict['SNVs'][i] = graph_dict['SNVs'][i] / float(graph_dict['length'][i])
		df = pd.DataFrame.from_dict(data=graph_dict)
		
	elif graphtype == 'scatter_ORF_ITEP_cluster':
	
		NUM_COLORS = 10 
		cm = plt.get_cmap('gist_rainbow')
		colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
		
	
		graph_dict = {'names' : [],
					  'n2n1' : [],
					  'SNVs' : [],
					  'length' : [],
					  'colors' : [],
					  'clusters' : []
					  }
		
		#clusters = {279:colors[0],124:colors[1],139:colors[2],4:colors[3],9:colors[4],154:colors[5]}
		
		ORFs = SeqObjList[4]
		
		for orf in ORFs:
			mean_n2n1 = 0
			count = 0
			cluster = orf.cluster
			if cluster != None:
				if len(orf.SNVs) > 0 and int(cluster) <= NUM_COLORS - 1:
				
					for snv in orf.SNVs:
						if mapping_coverage[snv.sample_id[:9]] > 10 :
							mean_n2n1 += snv.n2n1ratio
							count += 1
				
					if count > 0:
						mean_n2n1 = mean_n2n1 / float(count)
					graph_dict['names'].append(orf.ID)
					graph_dict['n2n1'].append(mean_n2n1)
					graph_dict['SNVs'].append(count)
					graph_dict['length'].append(orf.length)
					graph_dict['colors'].append(colors[int(cluster) - 1])
					graph_dict['clusters'].append(int(cluster))
			

		
		for i in range(len(graph_dict['n2n1'])):
			graph_dict['SNVs'][i] = graph_dict['SNVs'][i] / float(graph_dict['length'][i])
		df = pd.DataFrame.from_dict(data=graph_dict)
		

	return [df, graph_dict['colors']]

def get_cluster_dist(cluster):
	
	sum = 0
	count = 0
	for orf1 in cluster:
		for orf2 in cluster:
			d = np.sqrt((orf1[0] - orf2[0]) ** 2 + (orf1[1] - orf2[1]) ** 2)
			sum += d
			count += 1
	return sum / float(count)
		
def parse_ORFs_for_clusters(ORFs):
	clusters = {}
	for orf in ORFs:
		mean_n2n1 = 0
		count = 0
		cluster = orf.cluster
		#print(cluster)
		if cluster != None:
			if int(cluster) <= 500:
			
				if len(orf.SNVs) > 0:
					for snv in orf.SNVs:
						if mapping_coverage[snv.sample_id[:9]] > 10:
							mean_n2n1 += snv.n2n1ratio
							count += 1
			
					mean_n2n1 = mean_n2n1 / float(count)
		
				if not (cluster in clusters.keys()):
					clusters[cluster] = [(count, mean_n2n1)]
				else:
					clusters[cluster].append((count, mean_n2n1))
		
	graph_dict = {
				  'clusters' : [],
				  'distances' : []
	
				  }
	for key in clusters.keys():
	
		#clusters[key] = get_cluster_dist(clusters[key])
		print('cluster',key,'  number of ORFs:',len(clusters[key]),'  dist:',[get_cluster_dist(clusters[key])])
	
		if get_cluster_dist(clusters[key]) < 60:
			print('The following cluster has tight grouping and is composed of the following proteins: ')
			# for item in clusters[key]:
				#print(item[2])
	
		print('\n\n')
		graph_dict['clusters'].append(key)
		graph_dict['distances'].append(get_cluster_dist(clusters[key]))
		#clusters[key] = get_cluster_dist(clusters[key])
	
	
	
	
	df = pd.DataFrame.from_dict(data=graph_dict)
	bins = np.linspace(0, 300, 35)
	data=df.hist(column='distances',bins=bins,grid=False)
	
	return clusters
	
'''#########################################################################################################'''
'''#########################################################################################################'''

NUM_CLUSTERS = 1838
#all clusters with at least 3 ORFs

'''#########################################################################################################'''
'''#########################################################################################################'''
def make_ORF_cluster_dataframes(ORFs):
	graph_dict = {'names' : [],
				  'n2n1' : [],
				  'SNVs' : [],
				  'length' : [],
				  'clusters' : []
				  }
				  
	for orf in ORFs:

		mean_n2n1 = 0
		count = 0
		cluster = orf.cluster
		if cluster != None:
			if len(orf.SNVs) > 0:
		
				for snv in orf.SNVs:
					if mapping_coverage[snv.sample_id[:9]] > 10 \
					and snv.departure_from_reference > 0.05:
						mean_n2n1 += snv.n2n1ratio
						count += 1
		
				if count > 0:
					mean_n2n1 = mean_n2n1 / float(count)
				graph_dict['names'].append(orf.ID)
				graph_dict['n2n1'].append(mean_n2n1)
				graph_dict['SNVs'].append(count)
				graph_dict['length'].append(orf.length)
				graph_dict['clusters'].append(int(cluster))
				orf.mean_n2n1 = mean_n2n1
				orf.filtered_SNVs = count / float(orf.length)
		
		

	
	for i in range(len(graph_dict['n2n1'])):
		graph_dict['SNVs'][i] = graph_dict['SNVs'][i] / float(graph_dict['length'][i])
	
	return_dict = {}
	cluster_count = 1
	
	while cluster_count <= NUM_CLUSTERS:
	
		dfd = {'n2n1' : [],
			   'SNVs' : []
			   }
		
		for i in range(len(graph_dict['n2n1'])):
			
			if graph_dict['clusters'][i] == cluster_count:
				dfd['n2n1'].append(graph_dict['n2n1'][i])
				dfd['SNVs'].append(graph_dict['SNVs'][i])
				
		dfdf = pd.DataFrame.from_dict(data=dfd)
		return_dict[cluster_count] = dfdf
		
		cluster_count += 1
		
	# area_dict = {}
#	for key in return_dict.keys():
#		hull = ConvexHull(return_dict[key][['SNVs','n2n1']])
#		#plt.scatter(df['SNVs'], df['n2n1'], 'o')
#		# for simplex in hull.simplices:
# #			plt.plot(return_dict[key]['SNVs'].iloc[simplex], return_dict[key]['n2n1'].iloc[simplex], 'k-')
#			
#		area = hull.area
#		area_dict[key] = area
#		
#	return area_dict
	
def make_cluster_table(orflist, pa_dict, area_dict=None):
	cluster_output = open('cluster_output_noArea.txt', 'w')
	if area_dict != None:
		cluster_output.write('Cluster ID\tArea\tAverage SNVs\tAverage n2n1\tpa\n')
	else:
		cluster_output.write('Cluster ID\tAverage SNVs\tAverage n2n1\tpa\n')
	for c in range(1, NUM_CLUSTERS):
		ops = ''
		orfs_in_cluster = []
		for orf in orflist:
			if orf.cluster != None:
				if int(orf.cluster) == c:
					orfs_in_cluster.append(orf)
		mean_n2n1 = 0
		SNVs = 0
		
		for orf in orfs_in_cluster:
			mean_n2n1 += orf.mean_n2n1
			SNVs += orf.filtered_SNVs
	
		mean_n2n1 = mean_n2n1 / len(orfs_in_cluster)
		SNVs = SNVs / len(orfs_in_cluster)
		if area_dict != None:
			cluster_output.write(str(c) + '\t' + str(area_dict[c]) + '\t' + str(SNVs)  + '\t' + str(mean_n2n1)	+ '\t' + str(pa_dict[c]) + '\n')	
		else:
			cluster_output.write(str(c) + '\t' + str(SNVs)	+ '\t' + str(mean_n2n1)	 + '\t' + str(pa_dict[c]) + '\n')	
		
	cluster_output.close()

def summarize_cluster_file_noArea(cluster_file, average):
	
	cluster_file = open(cluster_file)
	if average == True:
		opf = open('cluster_summary_ave_noArea.txt', 'w')
	else:
		opf = open('cluster_summary_sum_noArea.txt', 'w')
	opf.write('SNVs\tn2n1\tSet\n')
	lines = cluster_file.readlines()
	opd = {}
	first = lines[0]
	for line in lines:
		if line != first:
			line = line.split('\t')
			if not(line[3] in opd.keys()):
			
				opd[line[3]] = [1, float(line[1]), float(line[2])]
				opd[line[3]][0] = 1
			else:
				opd[line[3]][0] += 1
				opd[line[3]][1] += float(line[1])
				opd[line[3]][2] += float(line[2])

	if average == True:
		for key in opd.keys():
			count = opd[key][0]
			opd[key][1] = float(opd[key][1]) / float(count)
			opd[key][2] = float(opd[key][2]) / float(count)
			
			ops = ''
			for ch in key:
				if ch == '[' or ch == ']' or ch == ',':
					x = 0
				elif ch == '0':
					ops += '0'
				elif ch == '1':
					ops += '1'
				else:
					x = 1
			
			opf.write(str(opd[key][1]) + '\t' + str(opd[key][2]) + '\t' + ops + '\n')
	else:
		for key in opd.keys():
			count = opd[key][0]
			opd[key][1] = float(opd[key][1])
			opd[key][2] = float(opd[key][2])
			
			ops = ''
			for ch in key:
				if ch == '[' or ch == ']' or ch == ',':
					x = 0
				elif ch == '0':
					ops += '0'
				elif ch == '1':
					ops += '1'
				else:
					x = 1
			
			opf.write(str(opd[key][1]) + '\t' + str(opd[key][2]) + '\t' + ops + '\n')		
	
	cluster_file.close()
	opf.close()


def summarize_cluster_file(cluster_file, average):
	
	cluster_file = open(cluster_file)
	if average == True:
		opf = open('cluster_summary_ave.txt', 'w')
	else:
		opf = open('cluster_summary_sum.txt', 'w')
	opf.write('Area\tSNVs\tn2n1\tSet\n')
	lines = cluster_file.readlines()
	opd = {}
	first = lines[0]
	for line in lines:
		if line != first:
			line = line.split('\t')
			if not(line[4] in opd.keys()):
				opd[line[4]] = [1, float(line[1]), float(line[2]), float(line[3])]
				opd[line[4]][0] = 1
			else:
				opd[line[4]][0] += 1
				opd[line[4]][1] += float(line[1])
				opd[line[4]][2] += float(line[2])
				opd[line[4]][3] += float(line[3])
	print(opd)	
	if average == True:
		for key in opd.keys():
			count = opd[key][0]
			opd[key][1] = float(opd[key][1]) / float(count)
			opd[key][2] = float(opd[key][2]) / float(count)
			opd[key][3] = float(opd[key][3]) / float(count)
			
			ops = ''
			for ch in key:
				if ch == '[' or ch == ']' or ch == ',':
					x = 0
				elif ch == '0':
					ops += '0'
				elif ch == '1':
					ops += '1'
				else:
					x = 1
			
			opf.write(str(opd[key][1]) + '\t' + str(opd[key][2]) + '\t' + str(opd[key][3]) + '\t' + ops + '\n')
	else:
		for key in opd.keys():
			count = opd[key][0]
			opd[key][1] = float(opd[key][1])
			opd[key][2] = float(opd[key][2])
			opd[key][3] = float(opd[key][3])
			
			ops = ''
			for ch in key:
				if ch == '[' or ch == ']' or ch == ',':
					x = 0
				elif ch == '0':
					ops += '0'
				elif ch == '1':
					ops += '1'
				else:
					x = 1
			
			opf.write(str(opd[key][1]) + '\t' + str(opd[key][2]) + '\t' + str(opd[key][3]) + '\t' + ops + '\n')		
	
	cluster_file.close()
	opf.close()
	
'''#########################################################################################################'''
'''#########################################################################################################'''

# def cluster(method, DataFrame, N):
# 
#	Num_clusters = N
#	cm = plt.get_cmap('gist_ncar')
#	colors = [str(mp.colors.to_hex(cm(1.*i/Num_clusters))) for i in range(Num_clusters)]
#	
#	if method == 'kmeans':
# 
#		kmeans = KMeans(n_clusters = Num_clusters)
#		X = DataFrame.as_matrix(columns=['SAAVs','n2n1'])
#		kmeans.fit(X)
#	
#		centroids = kmeans.cluster_centers_
#		labels = kmeans.labels_
#		plt.scatter(centroids[:, 0],centroids[:, 1], marker = "x", s=150)
#	
#		for i in range(len(X)):
#			#print(colors[labels[i]])
#			#print("coordinate:",X[i], "label:", labels[i])
#			plt.scatter(X[i][0], X[i][1], c=colors[labels[i]], alpha=0.5)
#			
#	elif method == 'dbscan':
#	
#		X = DataFrame.as_matrix(columns=['SNVs','n2n1'])
#		dbscan = DBSCAN(eps=0.1).fit(X)
# 
#		
#		#centroids = dbscan.cluster_centers_
#		labels = dbscan.labels_
#		#plt.scatter(centroids[:, 0],centroids[:, 1], marker = "x", s=150)
#		for i in range(len(X)):
#			#print(colors[labels[i]])
#			print("coordinate:",X[i], "label:", labels[i])
#			plt.scatter(X[i][0], X[i][1], markersize = 10, alpha=1,marker=markers[labels[i]])
#	
#	else:
#		pass
#		
#	plt.show()
#	#return centroids

def get_SNVs_by_sample(ORFs):
	sample_dict = {}
	count = 0
	for orf in ORFs:
		SNVs = orf.SNVs
		
		for snv in SNVs:
			SAG = snv.sample_id[:3]
			sample = snv.sample_id[4:9]
			
		if not (SAG in sample_dict.keys()):
			sample_dict[SAG] = {}
			sample_dict[SAG][sample] = [orf.cluster]
		else:
			if not (sample in sample_dict[SAG].keys()):
				sample_dict[SAG][sample] = [orf.cluster]
			else:
				if not (orf.cluster in sample_dict[SAG][sample]):
					sample_dict[SAG][sample].append(orf.cluster)
					print(SAG,'	  ',sample,'   ',orf.data)
				else:
					#print('Another added!',count)
					count += 1
	x = 0				
	for key in sample_dict.keys():
		
		for key1 in sample_dict[key].keys():
			
			x += len(sample_dict[key][key1])
			print('SAG: ',key,'sample: ',key1,'num clusters: ',len(sample_dict[key][key1]))
	
	testlist = []
	for SAG in sample_dict.keys():
		
		for sample in sample_dict[SAG].keys():
			if sample == 'FS866':
				testlist.append(sample_dict[SAG][sample])
	counter = 0
	samelist = []		
	for e1 in testlist:
		for e2 in testlist:
			if e1 != e2:
				x = list(set(e1).intersection(e2))
				print(len(x), len(e1), len(e2), float(len(x))/len(e2))
							
					
		
	#print(x)

def read_cluster_list_from_file(filename):
	file = open(filename)
	lines = file.readlines()
	return_list = []
	for line in lines:
		line = line.split()
		return_list.append(int(line[1]))
	
	file.close()
	return return_list
	
		

all_list = ['all','E23','N22','M21','K20','C09']

def get_SNVs_by_cluster(all_list, orflist):
	
	opf = open('exclusive_SNV_counts.txt', 'w')
	opf.write('File describing average numbers of SNVs per ORF for clusters exclusive to elements in SAG column\n')
	opf.write('SAG\tSNV count\n')
	for file in all_list:
		cluster_list = read_cluster_list_from_file(file + '_exclusive.txt')
		total = 0
		count = 0
		for c in cluster_list:
			for orf in orflist:
				if orf.cluster != None:
					if int(orf.cluster) == c:
						
						total += len(orf.SNVs)
						count += 1
					
		opf.write(file + '\t' + str(total / float(count)) + '\n')
						
		
def get_pa_text(cluster_output, orflist):
	
	alist = []
	
	for a in [0,1]:
		for b in [0,1]:
			for c in [0,1]:
				for d in [0,1]:
					for e in [0,1]:
						alist.append('[' + str(a) + ', ' + str(b) + ', ' + str(c) + ', ' + str(d) + ', ' + str(e) + ']')
	
	cluster_output1 = open(cluster_output, 'r')
	lines = cluster_output1.readlines()
	cluster_data_dict = {}
	for pa in alist:
		cluster_data_dict[pa] = []
	
		for orf in orflist:
			for line in lines[1:]:
				line = line.split('\t')
				cluster = line[0]
				if str(orf.cluster) == cluster:
					if line[-1][:-1] == pa:
						cluster_data_dict[pa].append(orf.data)
				
	print(cluster_data_dict)
	
	cluster_output1.close()
	
	for key in cluster_data_dict.keys():
		cluster_data_dict[key] = contrastText.get_unique(cluster_data_dict[key])
		print(key,'\n',cluster_data_dict[key],'\n\n\n\n')
				
		
def make_violin_plot_files(SeqObjList, SAG, pa_dict):
	#print(SAG) #function to make datafiles for violin plot script
	
	#switch between these file types by editing lines 1020-1030
	#opf = open(SAG + '_violin_SAAVs_locustagtext.txt', 'w')
	#opf = open(SAG + '_violin_SAAVs.txt', 'w')
	opf = open(SAG + '_violin_SAAVs_no_KO_data.txt', 'w')
	opf.write('orf sag/sample\tSNVs\n')
	summary_out = open(SAG + '_all_samples_KO_A_q.tsv', 'w')
	summary_out.truncate()
	
	mtf = open('output_mtf.txt', 'w')
	mtf.write('samp/SAG\tave\tKO_B\thi\tlo\n')
	my_Temp_Dict = {}
	print('Number of ORFs analyzed',len(SeqObjList[4]))
	yescount = 0
	for site in sample_sites:
		
		orf_names = []
		orf_SNV_densities = []
		#print('orfs',len(SeqObjList[4]))
		
		
		
		for orf in SeqObjList[4]:
			#print(SAG)
			ave_SNVs = 0
			ave_SAAVs = 0
			length = (orf.length / 1000.0)
			
			
			for SNV in orf.SNVs:

				if site in SNV.metagenome:

					ave_SNVs += 1

			for SAAV in orf.SAAVs:
				
				if site in SAAV.metagenome:

					ave_SAAVs += 1
					
					
					
			############################################################
			if ave_SAAVs == 0:
				if ave_SNVs == 0:
					yescount += 1
					print(yescount)
					ave = 0
				ave = 0
			elif ave_SAAVs > ave_SNVs:
				#print('An anomaly was identified,',ave_SAAVs, ave_SNVs)
				ave = 0
			else:
				ave = ave_SAAVs / ave_SNVs
				
				
				
			
			#print(ave)
			
			try:
				if (orf.alias) in my_Temp_Dict.keys():
					if ave > 0.5:
						my_Temp_Dict[orf.alias][-2] += 1
					else:
						my_Temp_Dict[orf.alias][-1] += 1
				
				else:
					if ave > 0.5:
						my_Temp_Dict[orf.alias] = [str(orf.alias),str(ave),str(orf.KO_B),1,0]
					
					else:
						my_Temp_Dict[orf.alias] = [str(orf.alias),str(ave),str(orf.KO_B),0,1]
			except TypeError:
				''
				
			try: 
				opf.write(SAG + '_' + site + '\t' + str(ave) + '\n')# + '\t' + orf.locus_tag + '\t' + orf.KO + '\t' + orf.KO_A + '\t' + orf.KO_B + '\t' + orf.KO_C + '\n')
			except TypeError:
				if orf.KO == None:
					
					opf.write(SAG + '_' + site + '\t' + str(ave) + '\t' + orf.locus_tag + '\t' + 'no annotation' + '\n')
				else:
					opf.write(SAG + '_' + site + '\t' + str(ave) + '\t' + orf.locus_tag + '\t' + orf.KO + '\n')
				
				
			
			orf_names.append(orf)
			orf_SNV_densities.append(ave_SNVs)
				
	
		orf_SNV_densities_np = np.array(orf_SNV_densities)
		upper = np.percentile(orf_SNV_densities_np,75.0)
		mean = np.mean(orf_SNV_densities_np)
		lower  = np.percentile(orf_SNV_densities_np,25.0)
		print('UPPER:',upper,'MEAN:',mean,'LOWER:',lower)
		hi = []
		lo = []
		mid = []
		for value in range(len(orf_SNV_densities)):
			if orf_SNV_densities[value] > upper:
				hi.append(orf_names[value])
			elif orf_SNV_densities[value] < (upper) and orf_SNV_densities[value] > (lower):
				mid.append(orf_names[value])
			elif orf_SNV_densities[value] < (lower):
				lo.append(orf_names[value])
			else:
				'Nothing'
		freq_dict_hi = {}
		freq_dict_mid = {}
		freq_dict_lo = {}
		
		for i in hi:
			try:
				if (i.KO_A in freq_dict_hi.keys()) and i.KO_A != None:
					freq_dict_hi[i.KO_A] += 1
				else:
					freq_dict_hi[i.KO_A] = 1
					
# 				print(i.KO_A,i.KO,i.KO_data,pa_dict[int(i.cluster)].count(1))
			except TypeError:
				print('no pa data available')
		
		for i in lo:
			try:
				if (i.KO_A in freq_dict_lo.keys()):
					freq_dict_lo[i.KO_A] += 1
				else:
					freq_dict_lo[i.KO_A] = 1
					
# 				print(i.KO_A,i.KO,i.KO_data,pa_dict[int(i.cluster)].count(1))
			except TypeError:
				print('no pa data available')
				
		for i in mid:
			try:
				if (i.KO_A in freq_dict_mid.keys()):
					freq_dict_mid[i.KO_A] += 1
				else:
					freq_dict_mid[i.KO_A] = 1
					
# 				print(i.KO_A,i.KO,i.KO_data,pa_dict[int(i.cluster)].count(1))
			except TypeError:
				print('no pa data available')
				
		summary_out.write(SAG + '_' + site + '\n\n')
		summary_out.write('\tHIGH SNV DENSITY\n')
		for key in freq_dict_hi:
			if key == None:
				summary_out.write('\t\t' + 'None' + ' : ' + str(freq_dict_hi[key]) + '\n')
			else:
				summary_out.write('\t\t' + key + ' : ' + str(freq_dict_hi[key]) + '\n')
		summary_out.write('\n')
		summary_out.write('\tMID SNV DENSITY\n')
		for key in freq_dict_mid:
			if key == None:
				summary_out.write('\t\t' + 'None' + ' : ' + str(freq_dict_mid[key]) + '\n')
			else:
				summary_out.write('\t\t' + key + ' : ' + str(freq_dict_mid[key]) + '\n')
		summary_out.write('\n')
		summary_out.write('\tLOW SNV DENSITY\n')
		for key in freq_dict_lo:
			if key == None:
				summary_out.write('\t\t' + 'None' + ' : ' + str(freq_dict_lo[key]) + '\n')
			else:
				summary_out.write('\t\t' + key + ' : ' + str(freq_dict_lo[key]) + '\n')
		summary_out.write('\n\n')	
		
	summary_out.close()
	opf.close()
	
	for key in my_Temp_Dict.keys():
		#print(my_Temp_Dict[key])
		str1 = str(my_Temp_Dict[key][0]) + ',' + str(my_Temp_Dict[key][1]) + ',' + str(my_Temp_Dict[key][2]) + ',' + str(my_Temp_Dict[key][3]) + ',' + str(my_Temp_Dict[key][4]) + '\n'
		#print([str1])
		mtf.write(str1)
			
			
			#print(length, ave_SNVs)
			
			
	mtf.close()	
					


					
			
				
		


if __name__ == '__main__':
	SAG_list = ['E23','C09','N22','K20','M21']
	SeqObjDict = {}
	for SAG in SAG_list:
		if os.path.isfile(SAG + '_SeqObj.pickle'):
			print('Unpacking data for SAG',SAG)
			list1 = unpickle_SeqObj(SAG)
			SeqObjDict[SAG] = list1
			
		else:
			list1 = make_sequence_objects(SAG)
			SeqObjDict[SAG] = list1
			pickle_SeqObj(list1)
			print('Packing data for SAG',SAG)
	
	pa_dict = mergeORFpaData.merge_PA_data('PresenceAbsenceTable_All.txt')
	for SAG in SAG_list:
		
		make_violin_plot_files(SeqObjDict[SAG], SAG, pa_dict)
		
		
	# df_list = []
#	color_list = []
#	orflist = []
#	for SAG in SAG_list:
#		
#		orflist += SeqObjDict[SAG][4]
#		
#	for orf in orflist:
#		index = 0
#		if len(orf.data) > 1:
#			index = orf.data[0].index('=', 26)
#		if index > 0:
#			newdata = orf.data[0][index + 1:] + ' '
#			for datapoint in orf.data[1:]:
#				newdata += datapoint + ' '
#				
#		orf.data = newdata
#	
#	alist = get_pa_text('cluster_output_noArea.txt', orflist)
		
	
	
				
				
		
		
	#get_SNVs_by_sample(orflist)
	#get_SNVs_by_cluster(all_list, orflist)
		
	#area_dict = make_ORF_cluster_dataframes(orflist)
	#pa_dict = mergeORFpaData.merge_PA_data('PresenceAbsenceTable_All.txt')
	#make_cluster_table(orflist, pa_dict)
	#summarize_cluster_file_noArea('cluster_output_noArea.txt', True)
	#summarize_cluster_file_noArea('cluster_output_noArea.txt', False)
	# for SAG in SAG_list:
# 		#df = make_DataFrame('bubble', SeqObjDict[SAG])
# 		#df = make_DataFrame('histo', SeqObjDict[SAG])
# 	
# 	
# 		returnList = make_DataFrame('bubble', SeqObjDict[SAG])
# 		df = returnList[0]
# 		colors = returnList[1]
# 		size = returnList[2]
# 		fig, ax = plt.subplots()
# 		df.plot.scatter(x='SNVs',y='n2n1',s=size,c=colors,alpha=0.5,ax=ax)
# 		#plt.savefig('SNVs_vs_n2n1_scatter_temp_'+SAG+'.svg')
# 		
# 	plt.show()
	
		#df_list.append(df)
		#color_list += colors

#	combined = pd.concat(df_list)
#	color_list = ['b']*15+['r']*15+['y']*14+['g']*12+['k']*14
#	combined.plot.scatter(x='SNVs',y='n2n1',c=color_list)
#	plt.savefig('SNVs_vs_n2n1_scatter_all_SAGs_TEST.svg')
# # 
#		
#	
#	
#	
#	
#	#combined.plot.scatter(x='SNVs',y='n2n1')
# # combined.plot.scatter(x='SAAVs',y='n2n1',c=color_list)
# # plt.show()
# # print('RUNNING CLUSTERING ALGORITHM')
# # cluster('kmeans',combined,5)
# # cluster('kmeans',combined,10)
# # cluster('kmeans',combined,15)
# # cluster('kmeans',combined,20)
# # 
#	
#	#combined.plot.scatter(x='SNVs',y='n2n1', c=colors) #for bubble plots
#	
#	#FOR n2n1 HISTOGRAM
# # fig = plt.figure()
# # ax1 = fig.add_subplot(111)
# # 
# # #FOR n2n1 HISTOGRAM
# # colors = ['b']*15+['r']*15+['y']*14+['g']*12+['k']*14
# # bins = np.linspace(0, 1, 500)
# # data=df_list[0].hist(column='n2n1',ax=ax1,bins=bins,grid=False,color='b',alpha=0.4)
# # data=df_list[1].hist(column='n2n1',ax=ax1,bins=bins,grid=False,color='r',alpha=0.4)
# # data=df_list[2].hist(column='n2n1',ax=ax1,bins=bins,grid=False,color='y',alpha=0.4)
# # data=df_list[3].hist(column='n2n1',ax=ax1,bins=bins,grid=False,color='g',alpha=0.4)
# # data=df_list[4].hist(column='n2n1',ax=ax1,bins=bins,grid=False,color='k',alpha=0.4)
#	
#	plt.show()
	
	
