investments_start <- read.csv("../data/investments.csv", header = TRUE, stringsAsFactors = FALSE)
library(plyr)

#Remove investments without investor name and companies that are listed as investors
investments <- investments_start[investments_start$investor_name != "",]
investors <- unique(investments$investor_name)
companies <- unique(investments$company_name)
overlap <- companies[companies %in% investors]
investments <- investments[!(investments$company_name %in% overlap),]

#Formulate SNAP graph data. (bipartite investors -> companies)
#This graph has multiple edges, meaning an investor has invested more than once.
#writes out to investmentsBipartite.txt

companies <- unique(investments$company_name)
investors <- unique(investments$investor_name)
everything <- c(companies, investors)
rangecompanies <- c(1:length(companies))
rangeinvestors <- c(1:length(investors)) + length(companies)

investorEdgeIds <- match(investments$investor_name, everything)
companyEdgeIds <- match(investments$company_name, everything)
bipartiteGraph <- data.frame(investors = investorEdgeIds, companies = companyEdgeIds)
write.table(bipartiteGraph, file = "../data/investmentsBipartite.txt", sep = "\t",row.names = FALSE)

#Formulate the folded graph
#Writes out to investorsFolded.txt

FindCooperatingInvestors = function(x, bipartiteGraph) {
  investorCoopIds <- bipartiteGraph[bipartiteGraph$companies == x,"investors"]
  investorCoopIds <- unique(investorCoopIds)
  if(length(investorCoopIds) > 1){
    investorCoopEdges <- combn(investorCoopIds,2)
    investorCoopEdges <- t(investorCoopEdges)
    investorCoopEdges <- data.frame(investorCoopEdges)
  }else{
    return(NA)
  }
  return(investorCoopEdges)
}

investorEdges <- lapply(rangecompanies,FindCooperatingInvestors, bipartiteGraph)
df <- ldply(investorEdges, data.frame)
investorEdgesDF <- df[!is.na(df$X1), c("X1","X2")]
write.table(investorEdgesDF, file = "../data/investorsFolded.txt", sep = "\t", row.names = FALSE)

#Write txt file that writes company name or investor name to nodeId companyIDs
nameTable <- data.frame(id = c(1:length(everything)), name = everything)
write.table(nameTable, file = "../data/companyIDs.txt", sep = "\t",row.names = FALSE)

#Generate timeseries data of bipartite and folded graph at each time

investments_time <- investments[investments$funded_quarter != "",]
funded_quarters <- sort(unique(investments_time$funded_quarter))
for(quarter in funded_quarters){
 if(substr(quarter,7,8) == "1"){
   investments_history <- investments_time[investments_time$funded_quarter <= quarter,]
   investorEdgeIds <- match(investments_history$investor_name, everything)
   companyEdgeIds <- match(investments_history$company_name, everything)
   bipartiteGraph <- data.frame(investors = investorEdgeIds, companies = companyEdgeIds)
   filename <- paste("../data/timeseriesBipartite/investmentsBipartite",quarter,".txt",sep="")
   write.table(bipartiteGraph, file = filename, sep = "\t",row.names = FALSE)
   
   investorEdges <- lapply(rangecompanies,FindCooperatingInvestors, bipartiteGraph)
   df <- ldply(investorEdges, data.frame)
   if(ncol(df) > 1){
     investorEdgesDF <- df[!is.na(df$X1), c("X1","X2")]
     filename <- paste("../data/timeseriesFolded/investmentsFolded",quarter,".txt",sep="")
     write.table(investorEdgesDF, file = filename, sep = "\t", row.names = FALSE)
   }
 }
}

#Generate folded graph that is unidirectional investor1 -> investor2 if investor 2 invests after investor 1
investorEdgeIds <- match(investments_time$investor_name, everything)
companyEdgeIds <- match(investments_time$company_name, everything)
bipartiteGraph <- data.frame(investors = investorEdgeIds, companies = companyEdgeIds, time = investments_time$funded_month)
write.table(bipartiteGraph, file = "../data/investmentsBipartiteTime.txt", sep = "\t",row.names = FALSE)

FindCooperatingInvestorsTime = function(x, bipartiteGraph){
  investorCoopIds <- bipartiteGraph[bipartiteGraph$companies == x,c("investors","time")]
  investorCoopIds <- investorCoopIds[order(investorCoopIds$time),]
  investorCoopIds <- investorCoopIds$investors
  if(length(investorCoopIds) > 1){
    investorCoopEdges <- combn(investorCoopIds,2)
    investorCoopEdges <- t(investorCoopEdges)
    investorCoopEdges <- data.frame(investorCoopEdges)
  }else{
    return(NA)
  }
  return(investorCoopEdges)
}
investorEdges <- lapply(rangecompanies,FindCooperatingInvestorsTime, bipartiteGraph)
df <- ldply(investorEdges, data.frame)
investorEdgesDF <- df[!is.na(df$X1), c("X1","X2")]
write.table(investorEdgesDF, file = "../data/investorsFoldedTime.txt", sep = "\t", row.names = FALSE)

