#module for analyzing python text

def contrast(text1, text2, threshold):
	length = len(text1)
	count = 0
	for ch in text1:
		if ch in text2:
			count += 1
			
	if (count / float(length)) > threshold:
		return True
	else:
		return False

def get_unique(returnlist):
	ls = []
	for element in returnlist:
		if not (element in ls):
			ls.append(element)
			
	return ls
		
testlist = ['UDP-glucose 4-epimerase',
'Nitrogenase molybdenum-iron protein alpha chain',
'Molybdate transport system permease protein',
'Hypothetical protein',
'Hypothetical protein',
'Nitrogen regulatory protein P-II family',
'Nitrogen regulatory protein P-II family',
'Membrane protein related to O-antigen export',
'Uncharacterized membrane protein',
'Nitrogenase molybdenum-iron protein beta chain',
'Nitrogenase molybdenum-iron protein beta chain',
'Molybdate transport system substrate-binding protein']

if __name__ == '__main__':
	returnlist = []
	for i in testlist:
		for j in testlist:
			if i != j:
				if contrast(i, j, 0.85):
					returnlist.append(i)
	print(get_unique(returnlist))
					
	