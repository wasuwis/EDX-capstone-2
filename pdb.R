if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(dslabs)) install.packages("dslabs")
if(!require(corrplot)) install.packages("corrplot")
if(!require(randomForest)) install.packages("randomForest")
#Structural Protein Sequences
#https://www.kaggle.com/shahir/protein-data-set

#download pdb_data_no_dups
dl <- tempfile()
download.file("https://raw.githubusercontent.com/wasuwis/EDX-capstone-2/master/pdb_data_no_dups.zip",dl)
pdb_data_no_dups <- read_csv(unzip(dl))

#download pdb_data_seq
dl_2 <- tempfile()
download.file("https://raw.githubusercontent.com/wasuwis/EDX-capstone-2/master/pdb_data_seq.zip",dl_2)
pdb_data_seq <- read_csv(unzip(dl_2))

#exclude chainid column because we dont use it in this project
pdb_data_seq <- pdb_data_seq %>%
  subset(select = -c(chainId)) %>%
  distinct()

#combine pdb_data_no_dubs with pdb__dat_seq
pdb <- left_join(pdb_data_no_dups,pdb_data_seq)

#exclude the following column due to lack of usesage and large amount of NAs
pdb <- pdb %>% subset(select=-c(crystallizationMethod,crystallizationTempK,pdbxDetails,publicationYear,structureId)) %>% 
  drop_na()

pdb <- as.data.frame(pdb)

# clear out phvalue more than 14
pdb <- pdb %>% filter(phValue<=14)


#Most use techniques is X-ray diffraction
pdb %>% group_by(experimentalTechnique) %>%
  tally(sort = TRUE)
pdb %>% group_by(experimentalTechnique) %>%
  tally(sort = TRUE) %>%
  .$n %>% .[1]/nrow(pdb)

#Top 10 numbers of classification
pdb %>% group_by(classification) %>%
  tally(sort = TRUE)

#for 15 first classification accounted for 60% of data that consisted with 3449 classification
for(i in 10:15) {
 
  x <- pdb %>% group_by(classification) %>%
   tally(sort = TRUE) %>% slice(1:i) %>%
   .$n %>% sum() 

 print(x/nrow(pdb))
 
}

#Number of protein in top 15 classification
pdb %>% group_by(classification) %>%
  tally(sort = TRUE) %>%
  slice(1:15) %>% 
  ggplot(aes(classification,n)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  ggtitle("Number of protein in top 15 classification")

#use only 15 classification
name_class <- pdb %>% group_by(classification) %>% tally(sort = TRUE) %>% slice(1:15) %>% .$classification  

#pH interpretation
ph_test <- function(phValue){
  
  if(phValue<7){return("Acid")}
  
  if(phValue>7){return("Base")}
  
  if(phValue==7){return("neutal")}

  }

#graph show number of acid,base and normal protein in database
pdb %>% group_by(phValue) %>%
  mutate(ph=ph_test(phValue)) %>%
  ungroup() %>% group_by(ph) %>%
  summarise(n=n()) %>%
  ggplot(aes(ph,n,fill=ph)) + 
  geom_bar(stat = "identity") +
  ggtitle("Distribution of acid, base and normal protein in database") +
  ylab("number of protein") +
  xlab("Type")

#Distribution of top 15 classification in acid, base and normal pH
pdb%>% filter(classification %in% name_class) %>%
  group_by(phValue) %>%
  mutate(ph=sapply(phValue,ph_test)) %>%
  ungroup() %>%
  ggplot(aes(ph,fill=classification)) + 
  geom_bar() +
  ggtitle("Distribution of top 15 classification in acid, base and neutral pH")



#plot show distibution of ph in top 15 classification
pdb %>% filter(classification %in% name_class) %>%
  ggplot(aes(phValue,fill=classification)) + 
  geom_density() +
  ggtitle("Distibution of ph in top 15 classifications")


#Macromolecule
pdb %>% group_by(macromoleculeType) %>%
  tally(sort = TRUE) %>%
  ggplot(aes(macromoleculeType,n)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Number in different macromolecule type") +
  ylab("Number")

#residue count
#plot show distribution of residue count
pdb %>% ggplot(aes(residueCount)) + 
  geom_density() +
  ggtitle("Distribution of residue count") 

#Distribution 0f residue count for top 15 classifications
pdb %>% filter(classification %in% name_class) %>%
  ggplot(aes(residueCount,fill=classification)) + 
  geom_density() +
  ggtitle("Distribution 0f residue count for top 15 classifications")

#Distribution of resolution
pdb %>% ggplot(aes(resolution)) + 
  geom_density()

#Distribution of resolution for top 15 classifications
pdb %>% filter(classification %in% name_class) %>% 
  ggplot(aes(resolution,fill=classification)) + 
  geom_density()

#Distribution of molecular weight
pdb  %>% filter(classification %in% name_class) %>%
  ggplot(aes(structureMolecularWeight,fill=classification)) + 
  geom_density() +
  ggtitle("Distribution of molecular weight in top 15 classification")

#we have to filter out some outlier that make graph unreadable
#plot the distribution of molecular weight under 1 million
pdb %>% filter(structureMolecularWeight<1*10^6) %>%
  filter(classification %in% name_class) %>%
  ggplot(aes(structureMolecularWeight,fill=classification)) + 
  geom_density() +
  ggtitle("Distribution of molecular weight under 1 million")

#plot the distribution of molecular weight under 500,000
pdb %>% filter(structureMolecularWeight<5*10^5) %>%
  filter(classification %in% name_class) %>%
  ggplot(aes(structureMolecularWeight,fill=classification)) + 
  geom_density() +
  ggtitle("Distribution of molecular weight under 500,000")

#Density Matthew
pdb %>% ggplot(aes(densityMatthews)) + 
  geom_density() +
  ggtitle("Distribution of density matthew")

#Distribution of Density matthew of top 15 classification
pdb %>% filter(classification %in% name_class) %>%
  ggplot(aes(densityMatthews,fill=classification)) + 
  geom_density() +
  ggtitle("Distribution of Density matthew of top 15 classification")

#density percent sol
pdb %>% ggplot(aes(densityPercentSol)) + 
  geom_density() +
  ggtitle("Distribution of Density percent sol")

pdb %>% filter(classification %in% name_class) %>%
  ggplot(aes(densityPercentSol,fill=classification)) + 
  geom_density() +
  ggtitle("distribution of density percent sol in top 15 classifications")

#Pearson correlation matrix
pdb_cor <- cor(pdb %>%
                 subset(select = -c(classification,experimentalTechnique,macromoleculeType,sequence)))
pdb_cor
corrplot(pdb_cor)


#we will keep only residue count, resolution, molecular weight, density matthew, density percent sol, phvalue for training purpose
temp_1 <- pdb %>%
  filter(classification %in% name_class) %>%
  subset(select=-c(sequence,experimentalTechnique,macromoleculeType))

# test set will be 20% of pdb data
set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)` instead
test_index <- createDataPartition(y = temp_1$classification, times = 1, p = 0.2, list = FALSE)
train_set <- temp_1[-test_index,]
temp_2 <- temp_1[test_index,]

# Make sure classification in test set are also in train set
test_set <- temp_2 %>% 
  semi_join(train_set, by = c("classification"))

# Add rows removed from test set back into train set
removed <- anti_join(temp_2, test_set)
train_set <- rbind(train_set, removed)

#knn
control_knn <- trainControl(method="cv",number = 10,p=0.9)
train_knn <- train(classification ~ . ,data = train_set,method="knn",trControl=control_knn,tuneGrid=data.frame(k=seq(3,9,2)))
ggplot(train_knn)
confusionMatrix(predict(train_knn,test_set),as.factor(test_set$classification))
results <- data.frame(method="knn",Accuracy=confusionMatrix(predict(train_knn,test_set),as.factor(test_set$classification))$overall[["Accuracy"]])
results

#Decision tree
train_rpart <- train(classification ~ . ,data = train_set,method="rpart",tuneGrid=data.frame(cp=seq(0.0,0.1,len=25)))
confusionMatrix(predict(train_rpart,test_set),as.factor(test_set$classification))
results <- rbind(results,data.frame(method="Decision tree",Accuracy=confusionMatrix(predict(train_rpart,test_set),as.factor(test_set$classification))$overall[["Accuracy"]]))
results
varImp(train_rpart)

#random forest
control_rf <- trainControl(method = "cv",number = 5)
train_rf <- train(classification ~ . ,data = train_set,method="rf",ntree=150,trControl=control_rf)
confusionMatrix(predict(train_rf,test_set),as.factor(test_set$classification))
results <- rbind(results,data.frame(method="Random Forest",Accuracy=confusionMatrix(predict(train_rf,test_set),as.factor(test_set$classification))$overall[["Accuracy"]]))
results
varImp(train_rf)






