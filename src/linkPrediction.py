from snap import *
from util import *
from collections import Counter
import math
import random

#G = loadDirectedGraph("../data/investmentsBipartite.txt")
F = loadUndirectedGraph("../data/investorsFolded.txt")
#Gt = loadDirectedGraph("../data/investmentsBipartiteTime.txt")
#Ft = loadDirectedGraph("../data/investorsFoldedTime.txt")

#ids = createNodeIdentificationMap("companyIDs.txt")

#metric that looks at number of same neighbors
def AdamicAdar(node1Id,node2Id, G, directed = False):
	node1 = G.GetNI(node1Id)
	node2 = G.GetNI(node2Id)
	neighbors1 = [node1.GetNbrNId(x) for x in range(node1.GetDeg())]
	neighbors2 = [node2.GetNbrNId(x) for x in range(node2.GetDeg())]
	neighborsOverlap = [G.GetNI(x) for x in neighbors1 if x in neighbors2]
	metric = sum([1/math.log(x.GetDeg()) for x in neighborsOverlap if x.GetDeg() > 0 ])
	return metric

#When using this function you should make sure that there are no self edges or else
#it will screw up the number of paths 
#I made sure not to include when there is already a path of 1 because that indicates the
#2 nodes are trivially connected.
#should increase beta so it is important to consider longer paths
#having high threshold makes it hard to run weightedKatz
def WeightedKatz(node1Id,node2Id, G, directed = False, threshold = 2, beta = 0.95):
	if(CntSelfEdges(G) > 0):
		print "THERE ARE SELF EDGES IN THE GRAPH. WEIGHTED KATZ METRIC WILL BE FAULTY."
	edgeCounts = Counter()
	currentNodes = Counter({node1Id: 1})
	for i in range(threshold):
		nextNodes = Counter()
		for(currentN, currentNcount) in currentNodes.items():
			nodeIt = G.GetNI(currentN)
			neighbors = [nodeIt.GetNbrNId(x) for x in range(nodeIt.GetDeg())]
			neighbors = Counter(neighbors)
			neighbors = Counter({x: y*currentNcount for x,y in neighbors.items()})
			nextNodes += neighbors

		print "Length next nodes", i , ":", len(nextNodes)
		edgeCounts[i+1] =  nextNodes[node2Id]
		currentNodes = nextNodes
		currentNodes[node2Id] = 0

	metric = 0
	for pathlength, numpaths in edgeCounts.items():
		if(pathlength > 1):
			metric += math.pow(beta,pathlength)*numpaths
	return metric

def generateSamplesAndData(G,filename,numsamples = 1000, directed = False, overWrite = False):
	GraphEdges = [(x.GetSrcNId(), x.GetDstNId()) for x in G.Edges()]
	GraphNodeIds = [x.GetId() for x in G.Nodes()]

	GraphSampleEdges = random.sample(GraphEdges,numsamples)
	GraphSampleNonEdges = []
	for i in range(numsamples):
		node1id, node2id = random.sample(GraphNodeIds,2)
		while node1id == node2id or G.IsEdge(node1id,node2id):
			node1id, node2id = random.sample(GraphNodeIds,2)
		GraphSampleNonEdges.append((node1id,node2id))

	if overWrite:
		txtfile = open(filename, 'w')
		txtfile.write("Node1\tNode2\tAdamicAdar\tWeightedKatz\tIsEdge\n")
	else:
		txtfile = open(filename, 'a')

	i = 1
	for (node1Id,node2Id) in GraphSampleEdges:
		print i
		aa = AdamicAdar(node1Id,node2Id,G)
		weightedKatz = WeightedKatz(node1Id,node2Id,G)
		output = [str(node1Id), str(node2Id), str(aa), str(weightedKatz), str(1)]
		output = '\t'.join(output)
		txtfile.write(output+"\n")
		i = i+1

	for (node1Id,node2Id) in GraphSampleNonEdges:
		print i
		aa = AdamicAdar(node1Id,node2Id,G)
		weightedKatz = WeightedKatz(node1Id,node2Id,G)
		output = [str(node1Id), str(node2Id), str(aa), str(weightedKatz), str(0)]
		output = '\t'.join(output)
		txtfile.write(output + "\n")
		i = i+1
	
	
generateSamplesAndData(F, "../data/FoldedNetworkUndirectedFeatures.tsv",numsamples=1000, overWrite = True)






