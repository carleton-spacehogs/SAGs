#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Determines frequency of SNVs in particular ORFs/KOs"""
from Contig import *
from SNV import *
from ORF import *
from SAAV import *
from make_SNVs_and_Contigs import *
from make_contig_ORF_and_SNV_lists import *
from organize_KO import *
import os



def sort_ORFs_by_SNV(ORFs):
    for fs in range(len(ORFs)-1,0,-1):
        pom=0
        for location in range(1,fs+1):
            if ORFs[location].number_of_SNVs > ORFs[pom].number_of_SNVs:
                pom = location

        temp = ORFs[fs]
        ORFs[fs] = ORFs[pom]
        ORFs[pom] = temp
    
    return ORFs

def assign_KO_to_ORFs(ORFs, d):
    for orf in ORFs:
        KO = orf.KO
        if KO != None:
            KO = orf.KO[3:]
        print 'THIS IS THE ORF KO',KO
        for A in d.keys():
            for B in d[A].keys():
                for C in d[A][B].keys():
                    for key in d[A][B][C].keys():
                        if KO == key:
                            orf.KO_A = A
                            orf.KO_B = B
                            orf.KO_C = C
                            orf.KO_data = d[A][B][C][key]
                            print "KO assigned with data:",orf.KO_data,orf.KO
                            
                            
    return ORFs
                            
def determine_KO_category_SNV_count(level, ORFs):
    
    d = {}
    
    for orf in ORFs:
        if level == 'A':
            x = orf.KO_A
            SNVs = orf.number_of_SNVs
            if x in d.keys():
                d[x][0] += SNVs
                d[x][1] += 1
            else:
                d[x] = [SNVs,1]
        elif level == 'B':
            x = orf.KO_B
            SNVs = orf.number_of_SNVs
            if x in d.keys():
                d[x][0] += SNVs
                d[x][1] += 1
            else:
                d[x] = [SNVs,1]
        elif level == 'C':
            x = orf.KO_C
            SNVs = orf.number_of_SNVs
            if x in d.keys():
                d[x][0] += SNVs
                d[x][1] += 1
            else:
                d[x] = [SNVs,1]
                
    for key in d.keys():
        d[key][0] = d[key][0] / float(d[key][1])
        print key,':',d[key][0]
        
    return d

def sort_ORFs_by_sample(level, sample, ORFs):
    
    d = {}
    
    for orf in ORFs:
        #print 'ORF',orf
        snvs = orf.SNVs
        #print 'SNVs',snvs
        for snv in snvs:
            #print snv.sample_id
            if snv.sample_id == sample:
                if level == 'A':
                    x = orf.KO_A
                    if x in d.keys():
                        d[x] += 1
                    else:
                        d[x] = 1
                
                elif level == 'B':
                    x = orf.KO_B
                    if x in d.keys():
                        d[x] += 1
                    else:
                        d[x] = 1
                    
                elif level == 'C':
                    x = orf.KO_C
                    if x in d.keys():
                        d[x] += 1
                    else:
                        d[x] = 1
                        
    print 'THIS IS d',[{}]
    
    for key in d.keys():
        print sample,':', key,':',d[key]

def write_SAAVs_for_graph(SAAVs):
    
    d = {}
    
    for saav in SAAVs:
        mapping = saav.sample_id
        if mapping not in d.keys():
            d[mapping] = 1
        else:
            d[mapping] += 1
            
    for key in d.keys():
        print key,'\t',d[key]
        
        
def get_SAAV_vs_SNV_for_ORF(ORFs, SAAVs):
    
    count = 0
    number = len(ORFs)
    num = 0
    
    for orf in ORFs:
        print num / float(number) * 100,'percent done'
        #print 'The first orf has a length of:',orf.length
        #print '... and is located on contig:',orf.contig
        for saav in SAAVs:
            if saav.contig_name == orf.contig:
                #print saav.gene_length, (orf.length + 100), (orf.length - 100)
                if int(saav.gene_length) < (int(orf.length) + 3):
                    if int(saav.gene_length) > (int(orf.length) - 3):
                        #print 'A MATCHING ORF/SAAV pair was found on contig %s' % (orf.contig)
                        #print 'ORF.length:',orf.length,'SAAV gene of origin length',saav.gene_length
                        orf.SAAVs.append(saav)
                        if len(orf.SAAVs) > len(orf.SNVs):
                            orf.flag = False
                        else:
                            orf.flag = True
                        count += 1
        num += 1
                        
    print 'A total of %s ORF/SAAV pairs were *probably* identified' % (count)
    
    return ORFs
    
def write_text_file_for_SAAV_v_SNV(ORFs):
    
    opf = open('outfile_E23.txt', 'w')
    
    for orf in ORFs:
        if orf.flag:
            opf.write(orf.ID + '\t' + str(len(orf.SAAVs)) + '\t' + str(len(orf.SNVs)) + '\t' + str(orf.KO_A) + '\t' + '\n')
    
    
    
    opf.close()
        

        
if __name__ == '__main__':
#    c = make_contigs('AD-155-C09-many_assemblers_simple.CISA.ctg.fa', 'C09_contig_names.tsv','C09.names_map')
#    s = make_SNVs('C09_variability_new.txt')
#    list1 = make_SNV_contig_names_consistent(s, c)
#    o = make_ORFs('2703719243.gff','2703719243.ko.tab.txt')
#    x = make_contig_ORF_and_SNV_lists(list1[1], list1[0], o)
#    d = parse_KO('ko00001.keg')
#    new_o = SNVs_in_ORFs(o, s, x)
#    assign_KO_to_ORFs(new_o, d)
#    
#    
#    d = determine_KO_category_SNV_count('A',new_o)
#    print 'END OF A'
#    d1 = determine_KO_category_SNV_count('B',new_o)
#    print 'END OF B'
#    d2 = determine_KO_category_SNV_count('C',new_o)
#    print 'END OF C'
#    #sorted_orfs = sort_ORFs_by_SNV(new_o)

    SAG_list_M21 = ['AD-155-M21-many_assemblers_simple.CISA.ctg.fa',
                'M21_contig_names.tsv',
                'M21.names_map',
                'M21_variability_new.txt',
                '2703719246.gff',
                '2703719246.ko.tab.txt']
    
    c = make_contigs('AD-155-E23-many_assemblers_simple.CISA.ctg.fa', 'E23_contig_names.tsv','E23.names_map')
    s = make_SNVs('E23_variability_new.txt')
    sa = make_SAAVs('E23_variability_AA.txt')
    list1 = make_SNV_SAAV_contig_names_consistent(s, c, sa)
    o = make_ORFs('2703719244.gff','2703719244.ko.tab.txt')
    x = make_contig_ORF_SAAV_SNV_lists(list1[1], list1[0], o, list1[2])
    d = parse_KO('ko00001.keg')
    new_o = SNVs_in_ORFs(o, s, x)
    newer_o = assign_KO_to_ORFs(new_o, d)
    
#    sort_ORFs_by_sample('A','C09_FS881_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS879_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS877_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS874_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS866_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS856_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS854_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS852_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS851_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS849_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS848_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS844_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS842_MAPPED_SORTED',newer_o)
#    sort_ORFs_by_sample('A','C09_FS841_MAPPED_SORTED',newer_o)
    
    write_SAAVs_for_graph(sa)
    saav_orfs = get_SAAV_vs_SNV_for_ORF(newer_o, sa)
    write_text_file_for_SAAV_v_SNV(saav_orfs)
    
    
    
    
#    d = determine_KO_category_SNV_count('A',new_o)
#    print 'END OF A'
#    d1 = determine_KO_category_SNV_count('B',new_o)
#    print 'END OF B'
#    d2 = determine_KO_category_SNV_count('C',new_o)
#    print 'END OF C'
    #sorted_orfs = sort_ORFs_by_SNV(new_o)
    
        
        