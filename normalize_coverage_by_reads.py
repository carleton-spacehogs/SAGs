import sys

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
def normalize(filename):	   
	outputfile = open(filename[:54] + '_normalized.txt','w')
	openfile = open(filename)
	lines = openfile.readlines()
	for line in lines:
		line = line.split()
		outputfile.write(line[0] + '\t' + line[1] + '\t' +  line[2] + '\t' + str(line[2] / metagenome_reads[filename[4:9]]) + '\n')
	
	
	openfile.close()
	outputfile.close()
	
if __name__ == '__main__':
	print(sys.argv[1])
	normalize(sys.argv[1])
	