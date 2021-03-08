from __future__ import absolute_import, division, print_function, unicode_literals
import pathlib
import matplotlib.pyplot as plt
import pandas
import seaborn
import numpy

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# Scott
import os
import argparse
import logging
import coloredlogs
import time
#--------------------------------------------------------------
# progress bar class for tensorflow output
class progressbar(keras.callbacks.Callback):
      def on_epoch_end(self, epoch, logs):
          # if epic
          if epoch % 10 == 0:
              print('.', end='')
          if epoch % 200 == 0:
              print(' '+str(epoch))

class PrintDot(keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs):
        if epoch % 100 == 0: print('')
        print('.', end='')

def build_NN(data_shape):
    """
    ARGS:
            data_shape (tuple): a tuple with the first number filled
                Ex: (1124,)

    Define a function to build our neural network.
    Creates the object for training.
    """
    model = keras.Sequential([
        layers.Dense(
            64,
            activation='relu',
            input_shape=data_shape
            ),
            layers.Dropout(0.1), # 10% dropout rate to prevent overfitting
            layers.Dense(
                64,
                activation='relu'
            ),
            layers.Dense(1) #final layer is 1 for the value we're predicting
        ])

      # use RMSprop optimizer with learning rate of 0.001
      # TODO modify?
    optimizer = tf.keras.optimizers.RMSprop(0.001)

      # Compile our model to prepare for use
    model.compile(loss='mse',
        optimizer=optimizer,
        metrics=['mae', 'mse'])
    return model

def train_NN(model_USNAM, train_geno_USNAM, train_pheno_USNAM,
            col_number, epochs, validation_split = 0.2, verbose = 0,
            my_callback = [progressbar()]):
    """
    Train the neural network
    Train on genotype data to predict column 'col_number' of the phenotype data
    """
    history_model_USNAM = model_USNAM.fit(
                            train_geno_USNAM,           # genotype data
                            train_pheno_USNAM[:,col_number], # phenotype data (1 column)
                            epochs = train_epoch,       # number of training epochs
                            validation_split = validation_split, # split the training data into 20% validation
                            verbose = verbose,                # no verbose messages
                            callbacks = my_callback # print progress bar
    )
    return history_model_USNAM

def save_models():
    pass
def def_seed(n=111519):
    """Define random seed for training data partitioning"""
    seed = n
    return seed

def def_frac(n=0.8):
    """Define fraction of data partitioning"""
    train_frac = n
    return train_frac

def def_epochs(n=300):
    """Define number of training epochs"""
    train_epoch = n
    return train_epoch


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

def validate_args(args, logger):
    """Raise if an input argument is invalid."""
    if not os.path.isdir(args.input_dir):
        logger.critical("argument 'input_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))

def return_corr(model_to_predict, train_geno, train_pheno, pheno_col):
    """
    ARGS:
        model_to_predict:
        train_geno:
        train_pheno:
        pheno_col: (int) of column number to subset

    Description:
    """
    corr = numpy.corrcoef(model_to_predict.predict(train_geno).flatten(),
                   train_pheno[:,pheno_col])[0][1]
    return corr

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
    logger.info('Importing genotype, this may take a moment...')
    geno_df_USNAM = load_geno_or_pheno(args.input_dir, geno_fname)
    logger.info('Importing phenotype, this may take a moment...')
    pheno_df_USNAM = load_geno_or_pheno(args.input_dir, pheno_fname)

    # DEFINE VARIABLES
    seed = def_seed()


    corr_columns = ['Train_Anthesis', 'Train_Kernel', 'Train_Height',
                    'Test_Anthesis', 'Test_Kernel', 'Test_Height']

    correlation_dataframe = pandas.DataFrame(columns=corr_columns)
    ultimate_data = []
    for i in range(1,51): #for each trial
        logger.info('Beginning trial loops')
        train_frac = def_frac()
        # TODO change as desired
        train_epoch = def_epochs(n=100)

        # Split data into train/testing sets
        logger.info('Fraction the data')
        train_geno_df_USNAM = geno_df_USNAM.sample(frac=train_frac, random_state=seed)
        train_pheno_df_USNAM = pheno_df_USNAM.loc[train_geno_df_USNAM.index,:]
        test_geno_df_USNAM = geno_df_USNAM.drop(train_geno_df_USNAM.index)
        test_pheno_df_USNAM = pheno_df_USNAM.drop(train_geno_df_USNAM.index)

        # Convert the data to a numpy array
        logger.info('Converting to Numpy array')
        train_geno_USNAM = convert_to_numpy(train_geno_df_USNAM)
        train_pheno_USNAM = convert_to_numpy(train_pheno_df_USNAM)
        test_geno_USNAM = convert_to_numpy(test_geno_df_USNAM)
        test_pheno_USNAM = convert_to_numpy(test_pheno_df_USNAM)

        #----------------------------------------------------------
        # Build the neural network
            # The argument passed to build_NN is just for the number of rows
            # This command creates the neural network object
            # Get the first axis of the shape and feed it to build_NN
            # MAGIC NUMBER
        logger.info('Building Neural Network')
        model_USNAM = build_NN((train_geno_USNAM.shape[1],))
        early_stop = keras.callbacks.EarlyStopping(monitor='val_loss', patience=20)

        #----------------------------------------------------------
        # Train the neural network
            # Train on genotype data to predict column 0 of the phenotype data

        # TODO get this more paramterized and use a for loop to generate the tables
        logger.info('Looping on the Phenotypes')
        data = []
        for i in range(0,3):
            history_model_USNAM = train_NN(model_USNAM,
                                        train_geno_USNAM,
                                        train_pheno_USNAM,
                                        i,
                                        train_epoch,
                                        my_callback = [early_stop, PrintDot()])

            # convert history to dataframe
            #history_df_USNAM = pandas.DataFrame(history_model_USNAM.history)

            # write history to file
            # TODO update te name to cover more filetypes and columns
            #history_df_USNAM.to_csv(args.output_dir+'/'+'Col_'+str(i)+'_'+history_fname, sep='\t')

            #print(model_USNAM.predict(train_geno_USNAM))
            #print(model_USNAM.predict(test_geno_USNAM))
            val = return_corr(model_USNAM, train_geno_USNAM, train_pheno_USNAM, i)
            data.append(val)
            val = return_corr(model_USNAM, test_geno_USNAM, test_pheno_USNAM, i)
            data.append(val)
        # MAGIC NUMBER
        assert len(data) == 6
        ultimate_data.append(data)
        #print(data)
        #print(correlation_dataframe.head()) 
    correlation_dataframe = pandas.DataFrame(ultimate_data,columns = corr_columns)
    correlation_dataframe.to_csv('Results/Correlation.tsv',
                            sep='\t',
                            index=False
                           )

            
