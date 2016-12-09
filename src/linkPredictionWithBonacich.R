FoldedDataTrain<- read.table("../data/FoldedNetworkUndirectedFeatures.tsv", header = TRUE, stringsAsFactors = FALSE)
BonacichCentrality <- read.table("../data/bonacich.dat", header = FALSE, stringsAsFactors = FALSE)
colnames(BonacichCentrality) <- c("NID", "BonacichCentrality")

FoldedDataTrain$BonCen1 <- BonacichCentrality[match(FoldedDataTrain$Node1, BonacichCentrality$NID),"BonacichCentrality"]
FoldedDataTrain$BonCen2 <- BonacichCentrality[match(FoldedDataTrain$Node2, BonacichCentrality$NID),"BonacichCentrality"]
FoldedDataTrain$StatusDiff <- abs(FoldedDataTrain$BonCen1 - FoldedDataTrain$BonCen2) 

logisticmodel <- glm(IsEdge ~ AdamicAdar+WeightedKatz+BonCen1+BonCen2+StatusDiff,family=binomial(link='logit'),data=FoldedDataTrain)
write.table(FoldedDataTrain, file = "../data/FoldedNetworkUndirectedFeatures2.tsv", quote = FALSE, sep = "\t",col.names = TRUE, row.names = FALSE)

#Statistics about logistic model
summary(logisticmodel)
anova(model, test = "Chisq")

#Find training error
fitted.results <- predict(logisticmodel,newdata=FoldedDataTrain,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

trainingMisClassificationError <- mean(fitted.results != FoldedDataTrain$IsEdge)
print(paste('Training Error:',trainingMisClassificationError))

#Do k-fold cross validation to find test error
#Randomly shuffle the data
FoldedDataTrain<-FoldedDataTrain[sample(nrow(FoldedDataTrain)),]

#Create 10 equally size folds
numFolds <- 10
folds <- cut(seq(1,nrow(FoldedDataTrain)),breaks=numFolds,labels=FALSE)

#Perform 10 fold cross validation
validationErrors <- c()
for(i in 1:numFolds){
  #Segment your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- FoldedDataTrain[testIndexes, ]
  trainData <- FoldedDataTrain[-testIndexes, ]
  #Use the test and train data partitions to find validation error
  crosslogisticmodel <- glm(IsEdge ~ AdamicAdar+WeightedKatz+BonCen1+BonCen2+StatusDiff,family=binomial(link='logit'),data=trainData)
  fitted.results <- predict(crosslogisticmodel,newdata=testData,type='response')
  fitted.results <- ifelse(fitted.results > 0.5,1,0)
  misClassificationError <- mean(fitted.results != testData$IsEdge)
  validationErrors[i] <- misClassificationError
}
validationError <- mean(validationErrors)
print(paste('Validation Error:',validationError))


