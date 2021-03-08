import numpy
import pandas
import matplotlib.pyplot as plt
import seaborn

from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, r2_score

# Scott
import os
import argparse
import logging
import coloredlogs
#--------------------------------------------

def def_seed(n=111519):
    """Define random seed for training data partitioning"""
    seed = n
    return seed

def load_geno_or_pheno(input_dir, my_file):
    """Load genotype and phenotype data as pandaframe"""
    loaded_panda_data = pandas.read_csv(
        os.path.join(input_dir, my_file),
        sep='\t',
        header=0,
        index_col=0)
    return loaded_panda_data

def convert_to_numpy(input_pandaframe):
    if isinstance(input_pandaframe, pandas.core.frame.DataFrame):
        return input_pandaframe.values

def def_frac(n=0.8):
    """Define fraction of data partitioning"""
    train_frac = n
    return train_frac

def build_rbf():
    """Build radial basis function kernel (rbf) model"""
    svr_rbf = SVR(
        kernel='rbf',
        C=100,
        gamma=0.1,
        epsilon=.1
    )
    return svr_rbf

def build_lin():
    """build linear kernel model"""
    svr_lin = SVR(
        kernel='linear',
        C=100,
        gamma='auto'
    )
    return svr_lin

def build_poly():
    """build polynomial kernel model"""
    svr_poly = SVR(
        kernel='poly',
        C=100,
        gamma='auto',
        degree=3,
        epsilon=.1,
        coef0=1
)

def train_rbf(svr_rbf, genotype, phenotype):
    """train radial basis function kernel (rbf) model"""
    svr_rbf = svr_rbf.fit(genotype, phenotype)
    return svr_rbf

def train_lin(svr_lin, genotype, phenotype):
    """train linear kernel (lin) model"""
    svr_lin = svr_lin.fit(genotype, phenotype)
    return svr_lin

def train_poly(svr_poly, genotype, phenotype):
    """train polynomial kernel (poly) model"""
    svr_poly = svr_poly.fit(genotype, phenotype)
    return svr_poly

def test_func():
    pass
    """






########################################
# step 5) test support vector machines
########################################

rbf_pred = svr_rbf.predict(test_geno_USNAM)
lin_pred = svr_lin.predict(test_geno_USNAM)
poly_pred = svr_poly.predict(test_geno_USNAM)

predictions = [
    ["Method", "MSE", "r_sq"],
    
    ["rbf",
     mean_squared_error(test_pheno_USNAM, rbf_pred
     r2_score(test_pheno_USNAM, rbf_pred)],
     
    ["linear",
     mean_squared_error(test_pheno_USNAM, lin_pred),
     r2_score(test_pheno_USNAM, lin_pred)],
     
    ["polynomial",
     mean_squared_error(test_pheno_USNAM, poly_pred),
     r2_score(test_pheno_USNAM, poly_pred)]
]

########################################
# step 6) write results to file
########################################

with open("svr_results.tsv", "w") as fout:
    for item in predictions:
        fout.write("%s\t%s\t%s\n" % (item[0], item[1], item[2]))

    """
def validate_args(args, logger):
    """Raise if an input argument is invalid."""
    if not os.path.isdir(args.input_dir):
        logger.critical("argument 'input_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))

if __name__ == '__main__':
    # FILE NAMES
    geno_fname = "US_NAM_markers_biallelic_SNPs.final.tsv"
    pheno_fname = "US_NAM_phenotypes.final.tsv"
    history_fname = "US_NAM_NN_history.tsv"

    parser = argparse.ArgumentParser(description="Start building models")
    path_main = os.path.abspath(__file__)
    parser.add_argument('input_dir', type=str,
                        help='Parent directory input data')
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.join(path_main, '../Results/'),
                        help='Parent directory to output results')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='set debugging level to DEBUG')
    args = parser.parse_args()
    args.output_dir = os.path.abspath(args.output_dir)
    args.input_dir = os.path.abspath(args.input_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    logger.info("Starting code... '%s'"%(args.input_dir))
    for argname, argval in vars(args).items():
        logger.debug("%-12s: %s"%(argname, argval))
    validate_args(args, logger)

    # LOAD DATA:
    # load the files; use first row and column for indexing.
    logger.info('Importing genotype, this may take a moment...')
    geno_df_USNAM = load_geno_or_pheno(args.input_dir, geno_fname)
    logger.info('Importing phenotype, this may take a moment...')
    pheno_df_USNAM = load_geno_or_pheno(args.input_dir, pheno_fname)

    # DEFINE VARIABLES
    seed = def_seed()
    train_frac = def_frac()

    # Split data into train/testing sets
    train_geno_df_USNAM = geno_df_USNAM.sample(frac=train_frac, random_state=seed)
    train_pheno_df_USNAM = pheno_df_USNAM.loc[train_geno_df_USNAM.index,:]
    test_geno_df_USNAM = geno_df_USNAM.drop(train_geno_df_USNAM.index)
    test_pheno_df_USNAM = pheno_df_USNAM.drop(train_geno_df_USNAM.index)

    # Convert the data to a numpy array
    train_geno_USNAM = convert_to_numpy(train_geno_df_USNAM)
    train_pheno_USNAM = convert_to_numpy(train_pheno_df_USNAM)
    test_geno_USNAM = convert_to_numpy(test_geno_df_USNAM)
    test_pheno_USNAM = convert_to_numpy(test_pheno_df_USNAM)

    # Build the models
    build_rbf = build_rbf()
    build_lin = build_lin()
    build_poly = build_poly()

    # Train the models
    train_rbf = train_rbf(build_rbf, train_geno_USNAM, train_pheno_USNAM)
    train_lin = train_lin(build_lin, train_geno_USNAM, train_pheno_USNAM)
    train_poly = train_poly(build_poly, train_geno_USNAM, train_pheno_USNAM)

    # Test
