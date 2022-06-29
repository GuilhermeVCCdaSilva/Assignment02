#!/usr/bin/env Rscript


install_BiocManager = function() {
    install.packages("BiocManager")
    BiocManager::install("pcaMethods")
    
}


read_data = function() {
    snp_data = read.csv(file = 'camels_project/camels_outfiles/camels_pca.str',
                    header=FALSE, 
                    sep=" ") 

    return(snp_data)
}


create_pca = function(snp_data) {
    camels_pca = pca(snp_data[,-c(1,2)], 
               scale="none", 
               center=T, 
               nPcs=4, 
               method="nipals") 
    str(camels_pca) 
    print(camels_pca@R2) 
    return(camels_pca)
}


plot_pca = function(snp_data,camels_pca) {
    slplot(camels_pca, 
       scol=snp_data[, c(2)], 
       scoresLoadings=c(TRUE,FALSE),
       sl=NULL, 
       spch=20) 

    legend("bottomright",
        legend=sort(unique(snp_data[, c(2)])),
        col=sort(unique(snp_data[, c(2)])),
        pch=20) 

}


install_BiocManager()
library(pcaMethods)
snp_data <- read_data()
camels_pca <- create_pca(snp_data)
plot_pca(snp_data,camels_pca)
