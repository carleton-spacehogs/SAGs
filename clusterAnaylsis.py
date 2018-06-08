






def get_cluster_dist(cluster):
	
	sum = 0
	count = 0
	for orf1 in cluster:
		for orf2 in cluster:
			d = np.sqrt((orf1[0] - orf2[0]) ** 2 + (orf1[1] - orf2[1]) ** 2)
			sum += d
			count += 1
	return sum / float(count)