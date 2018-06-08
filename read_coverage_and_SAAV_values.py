import seaborn
import pandas
import matplotlib.pyplot as plt
import matplotlib as mp
import math
SAGs = ['C09','E23','K20','M21','N22']
metagenomes = ['FS841','FS848','FS851','FS852','FS854','FS856','FS879','FS881']
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

def make_SAAV_dict(SAG_SAAV_filename):
	
	rtd = {}
	SAAV_file = open(SAG_SAAV_filename)
	lines = SAAV_file.readlines()
	
	# for line in lines:
# 		line = line.split()
# 		print(line)
	#print(lines)
	
	
	for line in lines[1:]:
		line = line.split('\t')
		if len(line) > 4:
			#print([line])
			if line[2] in rtd.keys():
				if line[0][4:] in metagenomes:
					rtd[line[2]]['ratio'].append(float(line[1]) * mapping_coverage[line[0]])
				
				
			else:
				rtd[line[2]] = {'ratio':[float(line[1]) * mapping_coverage[line[0]]],'KO':line[3],'KO_A':line[4],'KO_B':line[5],'KO_C':line[6][:-1]}
				
	SAAV_file.close()
	
	return rtd
	
def make_coverage_dict(coverage_filename, SAG_dict):
	covfile = open(coverage_filename)
	lines = covfile.readlines()
	for line in lines:
		line = line.split()
		if (line[1] != 'no annotation'):
			if (line[0] in SAG_dict.keys()):
				if 'coverages' in SAG_dict[line[0]].keys():
					SAG_dict[line[0]]['coverages'].append(float(line[3]))
				else:
					SAG_dict[line[0]]['coverages'] = [float(line[3])]
				
	covfile.close()
	return SAG_dict
 
def write_file(dict):
	opf = open('SAAV_coverage_test_new.txt','w')
	opf.write('genome' + '\t' + 'orf' + '\t' + 'KO' + '\t' + 'SAAV/SNV' + '\t' + 'RNA coverage' + '\t' + 'KO A' + '\t' + 'KO B' + '\t' + 'KO C\n')
	genomes = dict.keys()
	for genome in dict.keys():
		for orf in dict[genome].keys():
			#print(dict[genome][orf]['coverages'])
			opf.write(genome + '\t' + orf + '\t' + dict[genome][orf]['KO'] + '\t' + str(sum(dict[genome][orf]['ratio']) / len(dict[genome][orf]['ratio'])) + '\t' +  str(sum(dict[genome][orf]['coverages']) / len(dict[genome][orf]['coverages'])) + '\t' + dict[genome][orf]['KO_A'] + '\t' + dict[genome][orf]['KO_B'] + '\t' + dict[genome][orf]['KO_C'] + '\n')
	
	opf.close()

levels = ['A','B','C']
cutoffs = {'C09':{'x':10,'y':-25},
		   'E23':{'x':35,'y':-25},
		   'K20':{'x':10,'y':-25},
		   'M21':{'x':10,'y':-25},
		   'N22':{'x':4,'y':-25}}
	
def read_summary_file(summary_filename):
	
	openfile = open(summary_filename)
	lines = openfile.readlines()
	
		
	rtd = {}
	index = 0
	line4sub = 0
	for KO_level in levels:
		for sag in SAGs:
			outputfile = open(sag + '_' + KO_level + '_permuttest_new.txt','w')
			rtd = {}
			for line in lines:
				line4sub = 0
				line = line.split('\t')
				if line[0] == sag:
				
					if KO_level == 'A':
						index = 5
					if KO_level == 'B':
						index = 6
					if KO_level == 'C':
						index == 7
					
					if not (line[index] in rtd.keys()):
						rtd[line[index]] = {'1':0,'2':0,'3':0,'4':0}
					if line[4] == '0.0':
						line4sub = math.log(float(1e-15))
					else:
						line4sub = math.log(float(line[4]))
						
					
					
					if float(line[3]) > cutoffs[sag]['x']:
						if float(line4sub) > cutoffs[sag]['y']:
							rtd[line[index]]['2'] += 1
						else:
							rtd[line[index]]['3'] += 1
					else:
						if float(line4sub) > cutoffs[sag]['y']:
							rtd[line[index]]['1'] += 1
						else:
							rtd[line[index]]['4'] += 1
							
			outputfile.write(sag + '\t' + KO_level + '\n')
			outputfile.write('category' + '\t' + '1' + '\t' + '2' + '\t' + '3' + '\t' + '4' + '\n')
			for key in rtd.keys():
				if key[-1] != '\n':
					outputfile.write(key + '\t' + str(rtd[key]['1']) + '\t' + str(rtd[key]['2']) + '\t' + str(rtd[key]['3']) + '\t' + str(rtd[key]['4']) + '\n')
				else:
					outputfile.write(key[:-1] + '\t' + str(rtd[key]['1']) +  '\t' + str(rtd[key]['2']) + '\t' + str(rtd[key]['3']) + '\t' + str(rtd[key]['4']) + '\n')
			outputfile.close()
						
	openfile.close()
if __name__ == '__main__':
	dict = {}
	for sag in SAGs:
		dict[sag] = make_SAAV_dict(sag + '_violin_SAAVs_locustagtext.txt')
		
	for sag in SAGs:
		for m in metagenomes:
			dict[sag] = make_coverage_dict(sag + '_' + m + '_RNA_ORF_coverage_calculated_and_annotated_KO_normalized_new.txt',dict[sag])
			
			
	write_file(dict)
	read_summary_file('SAAV_coverage_test_new.txt')
	opefile = open('SAAV_coverage_test_new.txt')
	lines = opefile.readlines()
	
	colors = []
	color = 0

	color_label_dict = {'A':[],'B':[],'C':[]}
	
	for key in color_label_dict.keys():
		for line in lines[1:]:
			line = line.split();
		
			if key == 'A':
				if not (line[5] in color_label_dict[key]):
					color_label_dict[key].append(line[5])
			elif key == 'B':
				if not (line[6] in color_label_dict[key]):
					color_label_dict[key].append(line[6])
			elif key == 'C':
				if not (line[7] in color_label_dict[key]):
					color_label_dict[key].append(line[7])
					
					
	
	for key in color_label_dict.keys():
		
		NUM_COLORS_A = len(color_label_dict[key])
		cm = plt.get_cmap('gist_ncar')
		colors_A = [cm(1.*i/NUM_COLORS_A) for i in range(NUM_COLORS_A)] 
		
		
		
	
		
		for sag in SAGs:
			specific_colors = []
			x = []
			y =[]
			for line in lines[1:]:
				line = line.split()
				if line[0] == sag:
					if line[4] == '0.0':
						y.append(math.log(1.0e-15))
					else:
						y.append(math.log(float(line[4])))
					
					x.append(float(line[3]))
					
					if key == 'A':
						index = color_label_dict[key].index(line[5])
						specific_colors.append(colors_A[index])
					if key == 'B':
						index = color_label_dict[key].index(line[6])
						specific_colors.append(colors_A[index])
					if key == 'C':
						index = color_label_dict[key].index(line[7])
						specific_colors.append(colors_A[index])
						
						
						
			dict1 = {'x':x,'y':y}
			pd = pandas.DataFrame.from_dict(dict1)
			fart = seaborn.jointplot('x','y',data=pd,kind='reg')
			fart.ax_joint.cla()
			plt.sca(fart.ax_joint)
			plt.scatter(x, y, c=specific_colors)
			#plt.show()	
		#print(color_label_dict[key][i])
		for i in range(len(colors_A)):
			print(color_label_dict[key][i])
			print('color',colors_A[i],'   :    ',mp.colors.rgb2hex(colors_A[i]))
			#color += 1
			#plt.savefig(sag + '_double_scatter_reg_new.png')
			
	opefile.close()
			
			
			
						