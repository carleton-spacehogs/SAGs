def get_total_shared_genes():
	SAG_LIST = ['C09','K20','M21','N22','E23']
	pa_sum = open('cluster_output_noArea.txt')
	
	lines = pa_sum.readlines()
	clusters = []
	for line in lines[1:]:
		clusters.append((line.split('\t')[0], simplify(line.split('\t')[3])))
	
	
	
	
	
	for i in range(len(SAG_LIST)):
		count_total = 0
		for item in clusters:
					# print([item])
# 					print(item[1],item[1][i], i)
# 					print(item[1],item[1][j], j)
					if item[1][i] ==  '1':
						count_total += 1
		
		for j in range(len(SAG_LIST)):
			count = 0
			if i != j:
				for item in clusters:
					# print([item])
# 					print(item[1],item[1][i], i)
# 					print(item[1],item[1][j], j)
					if item[1][i] ==  '1' and item[1][j] == '1':
						
						count += 1
						
				print(SAG_LIST[i],'shares', count / float(count_total) * 100, '% of its genes with ',SAG_LIST[j])
				
def simplify(ast):
	 
	rts = ast[1] + ast[4] + ast[7] + ast[10] + ast[13]
	
	return rts
			

if __name__ == '__main__':
	get_total_shared_genes()