# By Scott Teresi
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


import os
import logging
import coloredlogs
import argparse
import glob

def validate_args(args, logger):
    """Raise if an input argument is invalid."""
    if not os.path.isdir(args.input_dir):
        logger.critical("argument 'input_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))

def load_pd(my_file):
    """
    Takes a pandaframe history object as argument.
    Previously generated in nn_pred.py
    """
    history = pd.read_csv(my_file, sep='\t')
    history.rename(columns={'Unnamed: 0': 'epoch'}, inplace = True)
    return history

def plot_history(hist, my_file, ylims_mae, ylims_mse):
    # TODO parameterize the ylabel for each phenotype
    #list_glob = glob.glob
    pheno_dict = {0: 'Anthesis GDD', 1: 'Kernel Row Number', 2: 'Plant Height'}


    plt.figure()
    plt.xlabel('Epoch')
    if my_file[4] == '0':
        plt.ylabel('Mean Abs Error: ' + str(pheno_dict[0]))
    elif my_file[4] == '1':
        plt.ylabel('Mean Abs Error: ' + str(pheno_dict[1]))
    elif my_file[4] == '2':
        plt.ylabel('Mean Abs Error: ' + str(pheno_dict[2]))
    else:
        plt.ylabel('Mean Abs Error: [Phenotype]')
    plt.plot(hist['epoch'], hist['mae'],
            label='Train Error')
    plt.plot(hist['epoch'], hist['val_mae'],
            label = 'Val Error')
    plt.yscale('log')

    plt.ylim(ylims_mae)
    plt.legend()
    plt.title('Neural Network')
    plt.savefig('Results/' + 'MAE_' + my_file)

    plt.figure()
    plt.xlabel('Epoch')
    if my_file[4] == '0':
        plt.ylabel('Mean Squared Error: ' + str(pheno_dict[0]))
    elif my_file[4] == '1':
        plt.ylabel('Mean Squared Error: ' + str(pheno_dict[1]))
    elif my_file[4] == '2':
        plt.ylabel('Mean Squared Error: ' + str(pheno_dict[2]))
    else:
        plt.ylabel('Mean Abs Error: [$Phenotype^2$]')
    plt.plot(hist['epoch'], hist['mse'],
            label='Train Error')
    plt.plot(hist['epoch'], hist['val_mse'],
            label = 'Val Error')
    plt.yscale('log')

    plt.ylim(ylims_mse)
    plt.legend()
    plt.title('Neural Network')
    plt.savefig('Results/' + 'MSE_' + my_file)
    plt.show()


if __name__ == '__main__':
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

    file_list = glob.glob(args.input_dir + '/*.tsv')
    for my_file in file_list:
        history = load_pd(my_file)
        ylims_mae = [history.mae.min(), history.mae.max()]
        ylims_mse = [history.mse.min(), history.mse.max()]
        print(history.tail())
        my_file = my_file.split('/')[-1]
        my_file = my_file.strip('.tsv')
        plot_history(history, my_file, ylims_mae, ylims_mse)

