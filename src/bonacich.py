# bonacich.py

import networkx as nx 
import collections
import matplotlib.pyplot as plt

__author__ = """Aric Hagberg <aric.hagberg@gmail.com>"""

def company_name():
	company_id = {}
	companyId = "../data/companyIDs.txt"
	f = open(companyId, 'r')
	for line in f: 
		if("id" not in line):
			linfo = line.strip().split("\t")
			company_id[int(linfo[0])] = linfo[1]
	return company_id

def get_folded_investors_graph():
	investors_folded = "../data/investorsFolded.txt"
	G = nx.Graph()
	nodes = set()
	edges = set()
	f = open(investors_folded, "r")
	for line in f:
		if("X" not in line):
			n1, n2 = map(int, line.strip().split("\t"))
			nodes.add(n1)
			nodes.add(n2)
			edges.add((n1, n2))

	for n in nodes:
		G.add_node(n)

	for n1, n2 in edges:
		G.add_edge(n1, n2)

	return G


def bonacich_centrality():
	f = open("bonacich.dat", 'w')
	company_id = company_name()
	investors_folded = get_folded_investors_graph()
	centrality = collections.Counter(nx.eigenvector_centrality(investors_folded))
	top_k = 10
	print("Top %d Ranked Bonacich Centrality" %(top_k))
	for c_id, cent_val in centrality.most_common(len(centrality)):
		f.write("%d\t%f\n" % (c_id, cent_val))

	# centrality_seq = sorted(nx.eigenvector_centrality(investors_folded).values(), reverse=True)
	# plt.loglog(centrality_seq, 'b-', marker='o')
	# plt.title("Bonacich Centrality Rank Plot")
	# plt.ylabel("Bonacich Centrality")
	# plt.xlabel("Rank")
	# plt.show()

def degree_distribution():
	investors_folded = get_folded_investors_graph()
	degree_sequence=sorted(nx.degree(investors_folded).values(),reverse=True) # degree sequence
	plt.loglog(degree_sequence,'b-',marker='o')
	plt.title("Degree rank plot")
	plt.ylabel("degree")
	plt.xlabel("rank")
	plt.show()

def cluster_coefficient():
	company_id = company_name()
	investors_folded = get_folded_investors_graph()
	clust_seq = sorted(nx.clustering(investors_folded).values(), reverse=True)
	plt.loglog(clust_seq, 'b-', marker = 'o')
	plt.title("Clustering Coefficient Rank Plot")
	plt.ylabel("Clustering Coefficient")
	plt.xlabel("Rank")
	plt.show()



if __name__ == "__main__":
	# degree_distribution()
	# cluster_coefficient()
	bonacich_centrality()


