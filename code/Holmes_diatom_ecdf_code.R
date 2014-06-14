library(vegan)

#load data: holmesdiatomsabundbylake.csv
Diatoms<-read.csv(file.choose(), header=T)
str(Diatoms)
DiatomMatrix<-as.matrix(Diatoms[,-1])
str(DiatomMatrix)

#eCDF - drop zeros, plot remaining species
v<-NULL
for(i in 1:dim(DiatomMatrix)[1]){
	tvec<-DiatomMatrix[i,]
	ind<-tvec!=0
	nozeroesvec<-tvec[ind]
	if(i==1) {plot(ecdf(nozeroesvec))}
	else {lines(ecdf(nozeroesvec))}
	}

for(i in 1:1){
	tvec<-DiatomMatrix[i,]
	ind<-tvec!=0
	nozeroesvec<-tvec[ind]
	if(i==1) {plot(ecdf(nozeroesvec))}
	else {lines(ecdf(nozeroesvec))}
	}