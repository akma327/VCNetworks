from snap import *
import numpy as np
import matplotlib.pyplot as plt
import collections
import operator

def loadDirectedGraph(inputName):
	G = LoadEdgeList(PNGraph, inputName, 0, 1)
	return G


def loadUndirectedGraph(inputName):
	G = LoadEdgeList(PUNGraph, inputName, 0, 1)
	return G 

def createNodeIdentificationMap(inputName):
	f = open(inputName, 'rb')
	headers = f.readline().split('\t')
	names = {}
	for line in f:
		tokens = line.strip().split('\t')
		assert(len(tokens) == 2)
		names[int(tokens[0])] = tokens[1:]

	return names

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

def getMaxInDegreeNodes(G, k = 1, offset = 0, keysOnly = False):
	nodeMap = {}
	for node in G.Nodes():
		nodeMap[node.GetId()] = node.GetInDeg()

	sortedNodes = sorted(nodeMap.items(), key=operator.itemgetter(1), reverse = True)
	topK = sortedNodes[offset:(offset + k)]
	if keysOnly:
		result = []
		for i in range(len(topK)):
			result.append(topK[i][0])
		return result

	return topK

def getMaxOutDegreeNodes(G, k = 1, offset = 0, keysOnly = False):
	nodeMap = {}
	for node in G.Nodes():
		nodeMap[node.GetId()] = node.GetOutDeg()

	sortedNodes = sorted(nodeMap.items(), key=operator.itemgetter(1), reverse = True)
	topK = sortedNodes[offset:(offset + k)]
	if keysOnly:
		result = []
		for i in range(len(topK)):
			result.append(topK[i][0])
		return result

	return topK


def getName(nodeId):
	return ids[nodeId]

def listTopInvestors(G, k = 1, offset = 0):
	topInvestors = getMaxOutDegreeNodes(G, k, offset, True)
	names = []
	for i in range(len(topInvestors)):
		names.append(ids[topInvestors[i]])
	return names

def listTopCompanies(G, k = 1, offset = 0):
	topInvestors = getMaxInDegreeNodes(G, k, offset, True)
	names = []
	for i in range(len(topInvestors)):
		names.append(ids[topInvestors[i]])
	return names

''' Main ''' 
# G represents our main directed bipartite graph
# F is our folded investor network 
G = loadDirectedGraph("../data/investmentsBipartite.txt")
F = loadUndirectedGraph("../data/investorsFolded.txt")
Gt = loadDirectedGraph("../data/investmentsBipartiteTime.txt")
Ft = loadDirectedGraph("../data/investorsFoldedTime.txt")
ids = createNodeIdentificationMap("../data/companyIDs.txt")

def outputConnectedComponent(G, nodeId):
	CnCom = TIntV()
	GetNodeWcc(G, nodeId, CnCom)
	print "Nodes in the same connected component as %s:" % ids[nodeId]
	for node in CnCom:
	    print node
	return CnCom

''' What is the node with the maximum degree? '''
def simulate():

	topNodes = getMaxOutDegreeNodes(G)
	print topNodes


	outputAvgClusteringCoefficients(G)		
	outputAvgClusteringCoefficients(F)	

	plotLogLogDegreeDistribution(G)
	plotLogLogDegreeDistribution(F)

	plotLogLogExcessDegreeDistribution(G)
	plotLogLogExcessDegreeDistribution(F)

	outputExpectedDegrees(G)
	outputExpectedDegrees(F)

	PlotOutDegDistr(G, "G-InDeg", "In degree distribution")
	PlotOutDegDistr(F, "F-InDeg", "In degree distribution")
	PlotOutDegDistr(Gt, "Gt-InDeg", "In degree distribution")
	PlotOutDegDistr(Ft, "Ft-InDeg", "In degree distribution")

	PlotInDegDistr(G, "G-InDeg", "In degree distribution")
	PlotInDegDistr(F, "F-InDeg", "In degree distribution")
	PlotInDegDistr(Gt, "Gt-InDeg", "In degree distribution")
	PlotInDegDistr(Ft, "Ft-InDeg", "In degree distribution")

	PrintInfo(G, "Statistics for G", "G-stats", False)
	PrintInfo(F, "Statistics for F", "F-stats", False)
	PrintInfo(Gt, "Statistics for Gt", "Gt-stats", False)
	PrintInfo(Ft, "Statistics for Ft", "Ft-stats", False)

	maxNode = getMaxOutDegreeNodes(G,1, 0, True)[0]



def retrieveSubgraphTopKInvestors(G, k = 1, offset = 0, degreeCutoff = 0):
	SubGraph_Nodes = TIntV()
	maxNodes = getMaxOutDegreeNodes(G, k, offset, True)
	for node in maxNodes:
		SubGraph_Nodes.Add(node)	
		for neighborId in G.GetNI(node).GetOutEdges():
			if neighborId not in SubGraph_Nodes:
				if G.GetNI(neighborId).GetInDeg() >= degreeCutoff:
					SubGraph_Nodes.Add(neighborId)

	return GetSubGraph(G, SubGraph_Nodes)


def visualizeGraph(G, outFile, ids = ids):
	NIdColorH = TIntStrH()
	labels = TIntStrH()

	for node in G.Nodes():
		companyName = ids[node.GetId()]
		print companyName[0]
		labels[node.GetId()] = companyName[0]
		if (node.GetOutDeg() == 0):
			NIdColorH[node.GetId()] = 'red'
		else:
			NIdColorH[node.GetId()] = 'yellow'

#	DrawGViz(G, gvlNeato, outFile, "Graph Visualization", True, labels, NIdColorH)
	DrawGViz(G, gvlNeato, outFile, "Graph Visualization", labels)

subGraph = retrieveSubgraphTopKInvestors(G, 10, 1, 25)	
visualizeGraph(subGraph, "outfile.png")
