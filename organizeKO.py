#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Makes dictionary for organizing KEGG Orthologies"""

def parse_KO(keg):
    
    d = {}
    keg = open(keg)
    keg_lines = keg.readlines()
    
    A = ''
    B = ''
    C = ''
    
    for i in range(8, len(keg_lines)):

        if keg_lines[i][0] == 'A':
            start = keg_lines[i].index('<b>')
            stop = keg_lines[i].index('</b>')
            A = keg_lines[i][(start + 3):stop]
            d[A] = {}
            
        elif keg_lines[i][0] == 'B':
            if len(keg_lines[i]) > 3:
                start = keg_lines[i].index('<b>')
                stop = keg_lines[i].index('</b>')
                B = keg_lines[i][(start + 3):stop]
                d[A][B] = {}
                
        elif keg_lines[i][0] == 'C':
            if len(keg_lines[i]) > 3:
                if 'BR' in keg_lines[i]:
                    start = keg_lines[i].index('0')
                    stop = keg_lines[i].index('BR')
                    C = keg_lines[i][(start + 5):stop-2]
                elif 'PATH' in keg_lines[i]:   
                    start = keg_lines[i].index('0')
                    stop = keg_lines[i].index('PATH')
                    C = keg_lines[i][(start + 5):stop-2]
                else:
                    start = keg_lines[i].index('0')
                    C = keg_lines[i][(start + 5):]
            d[A][B][C] = {}       

        elif keg_lines[i][0] == 'D':

            x = keg_lines[i].index('K')
            KO = keg_lines[i][x: x + 6]
            name = keg_lines[i][x + 8: -1]
            d[A][B][C][KO] = name
            
    return d

#Assigns KEGG orthologies to ORFs from organized dict
def assign_KO_to_ORFs(ORFs, d):
    for orf in ORFs:
        KO = orf.KO
        if KO != None:
            KO = orf.KO[3:]
        for A in d.keys():
            for B in d[A].keys():
                for C in d[A][B].keys():
                    for key in d[A][B][C].keys():
                        if KO == key:
                            orf.KO_A = A
                            orf.KO_B = B
                            orf.KO_C = C
                            orf.KO_data = d[A][B][C][key]
                            
                            
    return ORFs
	
        
if __name__ == '__main__':
    parse_KO('ko00001.keg')
    
    