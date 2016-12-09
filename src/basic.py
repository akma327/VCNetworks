from snap import *
import numpy as np
import matplotlib.pyplot as plt
import collections

def loadDirectedGraph(inputName):
	G = LoadEdgeList(PNGraph, inputName, 0, 1)
	return G


def loadUndirectedGraph(inputName):
	G = LoadEdgeList(PUNGraph, inputName, 0, 1)
	return G 

def computeDegreeDistribution(G):
	DegToCntV = TIntPrV()
	GetDegCnt(G, DegToCntV)

	x = [DegToCntV[i].GetVal1() for i in range(len(DegToCntV))]
	count = [DegToCntV[i].GetVal2() for i in range(len(DegToCntV))]

	total = sum(count)
	y = [freq * 1.0 / total for freq in count]

	return (x,y) 


def plotLogLogDegreeDistribution(G):
	x1, y1 = computeDegreeDistribution(G)
	plt.loglog(x1,y1, label = 'Graph')

	plt.title('log log plot of degree distribution for three networks')
	plt.xlabel('(log) degree')
	plt.ylabel('(log) degree probability ')
	plt.legend()
	plt.show()

def normalizeVector(v):
	total = sum(v)
	return v / total 

def normalizeDict(d):
	normalized = collections.defaultdict(float)
	total = sum(d.values())
	for k in d:
		normalized[k] = d[k] * 1.0 / total

	return normalized

def expectedExcessDegree(G):
	x, y = computeExcessDegree(G)
	result = 0.0
	for i in range(len(x)):
		result += x[i] * y[i]

	return result

def expectedDegree(G):
	x, y = computeDegreeDistribution(G)
	result = 0.0

	for i in range(len(x)):
		result += x[i] * y[i]

	return result


def computeExcessDegree(G):
	qp = collections.defaultdict(int) # q prime as defined in problem statement
	for node in G.Nodes():
		for neighborID in node.GetOutEdges():
			qp[G.GetNI(neighborID).GetDeg() - 1] += 1

	q = normalizeDict(qp)
	return (q.keys(), q.values())


def plotLogLogExcessDegreeDistribution(G):
	x1, y1 = computeExcessDegree(G)

	plt.loglog(x1,y1, label = 'Graph')

	plt.title('log log plot of excess degree distribution for three networks')
	plt.xlabel('(log) degree')
	plt.ylabel('(log) degree probability ')
	plt.legend()
	plt.show()

def outputExpectedDegrees(G):
	print "Graph has expected degree %f and expected excess degree %f" % (expectedDegree(G), expectedExcessDegree(G))	


def shapeLargestWCC(graph):
	wcc = GetMxWcc(graph)
	numNodes = wcc.GetNodes()
	numEdges = wcc.GetEdges()
	return (numNodes, numEdges)

def getWeaklyConnectedComponents(graph):
	components = TCnComV()
	wccs = GetWccs(graph, components)
	return components.Len()

def computeClusteringCoefficients(G):
	localCoeffs = {}
	for node in G.Nodes():
		k = node.GetDeg() 
		if k < 2:
			localCoeffs[node.GetId()] = 0
		else: 
			neighborEdges = 0
			for neighborID1 in node.GetOutEdges():
				for neighborID2 in node.GetOutEdges():
					if G.IsEdge(neighborID1, neighborID2):
						neighborEdges += 1
			neighborEdges /= 2
			localCoeffs[node.GetId()] = (2.0 * neighborEdges) / (k * (k-1))

	avgClusteringCoeff = sum(localCoeffs.values()) / G.GetNodes()

	return avgClusteringCoeff

def outputAvgClusteringCoefficients(G):
	print "Graph has average clustering coefficient of %f " % (computeClusteringCoefficients(G))	


def computeDegreeCentrality(G):
	degreeCentrality = {}

	for node in G.Nodes():
	    val = GetDegreeCentr(G, node.GetId())
	    degreeCentrality[node.GetId()] = val

	return degreeCentrality

def computeClosenessCentrality(G):
	# closeness centrality
	closenessCentrality = {}
	for node in G.Nodes():
	    val = GetClosenessCentr(G, node.GetId())
	    closenessCentrality[node.GetId()] = val

	return closenessCentrality

def computeBetweennessCentrality(G):
	Nodes = TIntFltH()
	Edges = TIntPrFltH()
	GetBetweennessCentr(G, Nodes, Edges, 1.0)

	betweennessCentrality = {}
	for node in Nodes:
		betweennessCentrality[node] = Nodes[node]

	return betweennessCentrality

''' Main ''' 

# G represents our main directed bipartite graph
# F is our folded investor network 
G = loadDirectedGraph("../data/investmentsBipartite.txt")
F = loadUndirectedGraph("../data/investorsFolded.txt")

outputAvgClusteringCoefficients(G)		
outputAvgClusteringCoefficients(F)	

plotLogLogDegreeDistribution(G)
plotLogLogDegreeDistribution(F)

plotLogLogExcessDegreeDistribution(G)
plotLogLogExcessDegreeDistribution(F)

outputExpectedDegrees(G)
outputExpectedDegrees(F)



