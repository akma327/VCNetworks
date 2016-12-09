library("igraph")
investor_path="/Users/akma327/Desktop/Fall 2016/CS224W/Project/Final/VCNetworks/data/investorsFolded.txt"
graph <- read_graph(investor_path, format="edgelist")
pc_1 <- power_centrality(graph, nodes = V(graph), loops = FALSE, exponent=1)
pc_75 <- power_centrality(graph, nodes = V(graph), loops = FALSE, exponent=0.75)
pc_0_5 <- power_centrality(graph, nodes= V(graph), exponent = 0.5)
pc_25 <- power_centrality(graph, nodes = V(graph), loops = FALSE, exponent=0.25)
pc_0 <- power_centrality(graph, nodes= V(graph), exponent = 0)
pc_neg25 <- power_centrality(graph, nodes = V(graph), loops = FALSE, exponent = 0.25)
pc_neg0_5 <- power_centrality(graph, nodes = V(graph), exponent = -0.5)
pc_neg75 <- power_centrality(graph, nodes = V(graph), exponent = -0.75)
pc_neg1 <- power_centrality(graph, nodes= V(graph), exponent = -1)


nid <- c()
centrality <-c()
count = 0
index = 1
k = 1
for (a in pc_neg75){
  if(a > 0){
    count = count + 1
    nid[k] <- toString(index)
    centrality[k] <- toString(a)
    k = k + 1
  }
  index = index + 1
}

nid_to_centrality <- cbind(nid, centrality)
df <- data.frame(nid_to_centrality)
write.table(nid_to_centrality, file = "/Users/akma327/Desktop/Fall 2016/CS224W/Project/Final/VCNetworks/data/bonacich/neg_75.txt")






