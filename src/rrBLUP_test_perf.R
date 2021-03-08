# load libraries
library(rrBLUP)

# set working directory
setwd("C:/Users/rs14/Documents/classes/hrt891a_NRT-IMPACTS_Foundations/final/US_NAM/")

# input genotype file path
geno_fpath = "US_NAM_markers_biallelic_SNPs.final.tsv"
# input phenotype file path
pheno_fpath = "US_NAM_phenotypes.final.tsv"
# training proportion
train_prop = 0.8
# number of tests for the rrBLUP method
nTest = 13
# output file path (stores correlation coefficients)
out_fpath = "US_NAM_rrBLUP_perf.tsv"

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
        header=T
    )
)

# seed rng using clock time
set.seed(as.integer(Sys.time()))

# make an empty data frame to store data
out_df <- data.frame(
    trial = integer(),
    train_cor_1 = double(),
    test_cor_1 = double(),
    train_cor_2 = double(),
    test_cor_2 = double(),
    train_cor_3 = double(),
    test_cor_3 = double()
)

# for each trial
for(i in 1:nTest) {
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

    # concatenate empty vector to new row in df
    # this will be modified by later steps
    out_df[i,] <- c(i, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # for each trait
    for(j in 1:3) {
        # conduct rrBLUP regression
        rrBLUP_fit <- mixed.solve(
            as.vector(train_pheno[,j]),
            Z=train_geno,
            K=NULL,
            SE=FALSE,
            return.Hinv=FALSE
        )

        # extract random effects coefficients for our marker data
        coeff <- as.matrix(rrBLUP_fit$u)

        # get intercept
        intercept <- as.vector(rrBLUP_fit$beta)

        # save the correlations
        # save training correlations
        out_df[i,2*j] <- cor(
            as.vector(train_pheno[,j]),
            as.vector((train_geno %*% coeff)[,1]) + intercept
        )
        # save testing correlations
        out_df[i,(2*j)+1] <- cor(
            as.vector(test_pheno[,j]),
            as.vector((test_geno %*% coeff)[,1]) + intercept
        )

        # write progress messages
        cat(j)
        flush.console()
    }
    # write data frame to file
    write.table(
        out_df,
        file=paste(
            "rr/",
            as.character(as.integer(Sys.time())),
            "_",
            out_fpath,
            sep=''
        ),
        sep='\t',
        row.names=FALSE
    )

    # print progress of processes
    cat("...", i, '\n', sep='')
    flush.console()
}
