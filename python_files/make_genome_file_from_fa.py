#This script makes an genome file in the format
#contig-name \t contig-size
#from a fasta file.

def read(fasta):
    
    myfile = open(fasta, 'r')
    lines = myfile.readlines()
    opf = open(fasta[:10]+'_genome.txt', 'w')
    
    for i in range(len(lines)):
                          
        if lines[i][0] == '>':
            print lines[i]
            ops = lines[i][1:-1] + '\t' + str(len(lines[i+1]) - 1) + '\n'
            print 'THIS IS OPS',[ops]
            opf.write(ops)
            ops = ''
    opf.close()
    myfile.close()
    
if __name__ == '__main__':
    read('AD-155-E23-many_assemblers_simple.CISA.ctg.fa')
    read('AD-155-C09-many_assemblers_simple.CISA.ctg.fa')
    read('AD-155-N22-many_assemblers_simple.CISA.ctg.fa')
    read('AD-155-M21-many_assemblers_simple.CISA.ctg.fa')
    read('AD-155-K20-many_assemblers_simple.CISA.ctg.fa')
            
            