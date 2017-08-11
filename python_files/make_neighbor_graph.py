import networkx as nx
import matplotlib.pyplot as plt

def make_neighbor_graph(input_text_name):
	G = nx.Graph()
	start_list = []
	stop_list = []
	myfile = open(input_text_name)
	lines = myfile.readlines()
	for line in lines:
		line = line.split()
		start_list.append(line[0])
		stop_list.append(line[1])
		
	for i in range(len(start_list)):
		G.add_node(start_list[i])
		print('start',start_list[i])
		G.add_node(stop_list[i])
		print('stop',stop_list[i])
		G.add_edge(start_list[i], stop_list[i])
		
		
	nx.draw_networkx(G, with_labels=False,node_size=100)
	
	plt.show()
	
def test():
	G = nx.Graph()
	G.add_node(1)
	G.add_node(2)
	G.add_edge(1,2)
	nx.draw_networkx(G)
	plt.show()
if __name__ == '__main__':
	#test()
	#make_neighbor_graph('C09_M21neighbors.txt')
	#make_neighbor_graph('C09_K20_M21_neighborhood.txt')
	make_neighbor_graph('E23_N22neighbors.txt')
		
		