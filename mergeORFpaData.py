#function for analyzing gene presence/absence table from ITEP
#Michael Hoffert


def merge_PA_data(PA_file):
	
	pa_file = open(PA_file)
	pa_lines = pa_file.readlines()
	first = pa_lines[0]
	pa_dict = {}
	for line in pa_lines:
		if line != first:
			line=line.split('\t')
			cluster = int(line[1])
			pa = line[3:]
			opl = [0,0,0,0,0]
			for i in range(5):
				#print(pa[i])
				if pa[i] == 'NONE' or pa[i] == 'NONE\n':
					opl[i] = 0
				else:
					opl[i] = 1
			#print(cluster,'\t',opl)
			pa_dict[cluster] = opl
			
	return pa_dict
		
			
		
	


	
	pa_file.close()
	
if __name__ == '__main__':
	merge_PA_data('PresenceAbsenceTable2.txt')