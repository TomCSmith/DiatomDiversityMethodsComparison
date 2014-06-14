library(vegan)

### Exploratory analysis of Holmes data...looking at the ecdfs ###

# Read in Holmes data...I think this is Windows compatible
Diatoms <- read.csv(file.path("..", "data", "archival", "HolmesDiatomsAbundByLake.csv"),
                    header=T)

DiatomMatrix <- as.matrix(Diatoms[,-1])
S_vec = array(NA, )

# eCDF - drop zeros, plot remaining species
for(i in 1:dim(DiatomMatrix)[1]) {

    # Drop zeros from SAD
	tvec <- DiatomMatrix[i,]
	ind <- tvec != 0
	nozeroesvec <- tvec[ind]
    nozeroesvec = sort(nozeroesvec)

    ecdf_fxn = ecdf(nozeroesvec)
	if(i==1) {
        plot(log(nozeroesvec), ecdf_fxn(nozeroesvec), type="b", pch=21, col=i,
                    xlab="log(Abundance)", ylab="CDF")
    }
	else {
        lines(log(nozeroesvec), ecdf_fxn(nozeroesvec), type="b", pch=21, col=i)
    }
}


