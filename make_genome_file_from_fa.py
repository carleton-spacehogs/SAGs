#This script makes an genome file in the format
#contig-name \t contig-size
#from a fasta file.
import sys

def read(fasta):
    
    myfile = open(fasta, 'r')
    lines = myfile.readlines()
    opf = open(fasta[:10]+'_genome.txt', 'w')
    
    for i in range(len(lines)):
                          
        if lines[i][0] == '>':
            ops = lines[i][1:-1] + '\t' + str(len(lines[i+1]) - 1) + '\n'
            opf.write(ops)
            ops = ''
    opf.close()
    myfile.close()
    
if __name__ == '__main__':
    read(sys.argv[1])
            
            