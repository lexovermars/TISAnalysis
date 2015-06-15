args <- commandArgs(TRUE)
file_name <- args[1]
file_name_input <- paste("output/",file_name,"_3_max_length.txt",sep="")

#Do the initial PCA
data <- read.table(file_name_input,header=TRUE)

#Check here which columns you want to include (dependent on vector length in PCA_analysis.py)
#complete with startnt
fit <- prcomp(data[,7:154],center=TRUE)

all_data <- cbind(as.character(data$Start_label),fit$x[,1:2])
plot(all_data[which(all_data[,1]=="annotated"),2:3],pch='.',cex=3,col="red")
points(all_data[which(all_data[,1]=="downstream"),2:3],pch='.',cex=3,col="blue")
points(all_data[which(all_data[,1]=="upstream"),2:3],pch='.',cex=3,col="black")

#Make a matrix with the loadings and write them to a file
library(MASS)
complete <- cbind(as.character(data$Gene),as.character(data$Pos_rel),all_data)
file_name_pca_scores <- paste("R_output/",file_name,"_pca_scores.txt",sep="")
file_name_pca_scores
write.matrix(complete,file=file_name_pca_scores,sep="\t")
#Write PCA loadings
file_name_loadings <- paste("R_output/",file_name,"_pca_loadings.txt",sep="")
loading_matrix <- cbind(as.character(rownames(fit$rotation)),fit$rotation[,1:3])
write.matrix(loading_matrix,file=file_name_loadings,sep="\t")

#project all protential starts on PC1 of the initial PCA
file_name_input <- paste("output/",file_name,"_all_starts.txt",sep="")
data_all_starts <- read.table(file_name_input,header=TRUE)
projected_scores <- predict(fit,data_all_starts)
file_name_projected_output <- paste("R_output/",file_name,"_pca_scores_projected.txt",sep="")
projected_matrix <- cbind(as.character(data_all_starts$Gene),data_all_starts$Pos_start,data_all_starts$Pos_rel,as.character(data_all_starts$Start_label),projected_scores[,1:3])
write.matrix(projected_matrix,file=file_name_projected_output,sep="\t")
