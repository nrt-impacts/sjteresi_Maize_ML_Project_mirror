# load libraries
library(BGLR)

geno_fpath = "US_NAM_markers_biallelic_SNPs.final.tsv"
pheno_fpath = "US_NAM_phenotypes.final.tsv"
train_prop = 0.6
nIter = 1500
burnIn = 1000

# load marker data
geno <- as.matrix(
            read.csv(
            file=geno_fpath,
            sep='\t',
            row.names=1,
            header=T
        )
    )

# load phenotype data
pheno <- as.matrix(
    read.csv(
        file = pheno_fpath,
        sep='\t',
        row.names=1,
        header=T,
    )
)

# get dimensions of matrix
dim(geno)
dim(pheno)

# set random seed for data row selection reproducibility
set.seed(112519)

# generate a vector of sample indices
selections <- sample(
    1:nrow(geno),
    as.integer(ceiling(train_prop * nrow(geno))),
    replace = FALSE
)

# extract train, test genotypes
train_geno <- geno[selections,]
test_geno <- geno[-selections,]

# extract train, test phenotypes
train_pheno <- pheno[selections,]
test_pheno <- pheno[-selections,]

# conduct BayesA regression
fmBA <- BGLR(
    y = as.vector(train_pheno[,1]),
    ETA = list(
        list(
            X=train_geno,
            model='BayesA'
        )
    ),
    nIter=nIter,
    burnIn=burnIn,
    saveAt='ba_'
)

# calculate estimators for testing data
test_pheno_hat = as.vector(test_geno %*% fmBA$ETA[[1]]$b + fmBA$mu)

# plot the correlation between testing predictions and true values
plot(test_pheno_hat, test_pheno[,1])

# calculate the correlation between testing predictions and true values
cor(test_pheno_hat, test_pheno[,1])

# code snippets from demo code
#plot(abs(fmBA$ETA[[1]]$b),col=4,cex=.5, type='o',main='BayesA')
#abline(v=QTL,col=2,lty=2)
