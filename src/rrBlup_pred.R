# load libraries
library(rrBLUP)

geno_fpath = "US_NAM_markers_biallelic_SNPs.final.tsv"
pheno_fpath = "US_NAM_phenotypes.final.tsv"
train_prop = 0.6

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

# set random seed
set.seed(112519)

# generate a vector of sample indices
selections <- sample(
    1:(nrow(geno)-1),
    as.integer(ceiling(train_prop * nrow(geno))),
    replace = FALSE
)

# extract train, test genotypes
train_geno <- geno[selections,]
test_geno <- geno[-selections,]

# extract train, test phenotypes
train_pheno <- pheno[selections,]
test_pheno <- pheno[-selections,]

# mixed.solve is the key function of rrBLUP
anthesis_fit <- mixed.solve(
    as.vector(train_pheno[,1]),
    Z=train_geno,
    K=NULL,
    SE=FALSE,
    return.Hinv=FALSE
)

# extract random effects coefficients for our marker data
coeff <- as.matrix(anthesis_fit$u)

# make predictions on our test data
blue_test_pheno <- test_geno %*% coeff
blup_test_pheno <- as.vector(blue_test_pheno[,1]) + as.vector(anthesis_fit$beta)

# make accuracy calculations
anthesis_accuracy <- cor(
    blup_test_pheno,
    test_pheno[,1],
    use="complete"
)

print(anthesis_accuracy)
