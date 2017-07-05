source("halt.R")

# this is short script to check the data and suggest minimal support

options <- commandArgs(trailingOnly = T)
WD <- getwd();
if (!is.null(WD)) setwd(WD);

input  <- options[1]

if (is.na(input)) {
  print("No data given")
  halt()
}

# read the matrix and extract only meaningful sites and species

  data <- read.table(input,sep=",")
  my_matrix <- data.matrix(data)

# take only the sites with >1 bacteria and bacteria that occur in >2 different sites

  while (length(which(apply(my_matrix,2,sum)<3))+length(which(apply(my_matrix,1,sum)<2))>0) {
    my_matrix <- my_matrix[,which(apply(my_matrix,2,sum)>2)]
    my_matrix <- my_matrix[which(apply(my_matrix,1,sum)>1),]
  }
  places <- dimnames(my_matrix)[[1]]
  names <- dimnames(my_matrix)[[2]]

# provide rows and columns names is needed 

  if (length(places)==0) {
    places <- vector("character",length(my_matrix[,1]))
    for (i in 1:length(my_matrix[,1])) {
      places[i] <- paste("SITE",as.character(i),sep="_")
    }
    dimnames(my_matrix)[[1]] <- places
  }

  if (length(names)==0) {
    names <- vector("character",length(my_matrix[1,]))
    for (i in 1:length(my_matrix[1,])) {
      names[i] <- paste("SPECIE",as.character(i),sep="_")
    }
	dimnames(my_matrix)[[2]] <- names
  }  

# calculate the minimal meaningfull support

  sparsity <- length(which(my_matrix==1))/length(my_matrix)
  freq <- apply(my_matrix,2,sum)/length(places)
  abundancy <- apply(my_matrix,1,sum)/length(names)
  Fr <- mean(freq)
  min.support <- round(round(Fr*Fr,digits = 3) * length(places),0)
  if (min.support<3)
    min.support = 3
  print(paste("Sparsity of the matrix: ", round(sparsity*100,digits=2), "%",sep=""))
  print(paste("Mean frequency of species: ",round(Fr*100,digits=2),"%",sep=""))
  print(paste("Mean presence of species in sampling sites: ",round(mean(abundancy)*100,digits=2),"%",sep=""))
  print(paste("Recommended min.support: ",min.support," common sites",sep=""))

