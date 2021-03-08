# By Scott Teresi
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_pd(my_file):
    """
    ARGS:
        my_file: (string) Filename with folder name extension in front. The
        dataframe has rows that are individual trials with each column being
        a phenotype and training or testing status. Each value is a correlation
        value.
    """
    data = pd.read_csv(my_file, sep='\t', header = 'infer')
    return data


def seaborn_multi(data, show=False):
    """
    ARGS:
        my_data: (pandas.core.frame.DataFrame) Pandas DataFrame object.
        Description of dataframe in load_pd.
    OBJ:
        Compare the performance between all of the methods
    """
    sns.boxplot(x = 'Traits', y = 'Cor', hue='Method', data =
                data.loc[data['Set'] == 'Train'])
    plt.xlabel('Phenotype')
    plt.ylabel('Correlation')
    plt.title('Training Model Performance')
    plt.savefig('Correlations/Train_MultiBox.png')
    if show:
        plt.show()
    plt.clf()

    sns.boxplot(x = 'Traits', y = 'Cor', hue='Method', data =
                data.loc[data['Set'] == 'Test'])
    plt.xlabel('Phenotype')
    plt.ylabel('Correlation')
    plt.title('Testing Model Performance')
    plt.savefig('Correlations/Test_MultiBox.png')
    if show:
        plt.show()

def seaborn_singular(data, show=False):
    """
    ARGS:
        my_data: (pandas.core.frame.DataFrame) Pandas DataFrame object.
        Description of dataframe in load_pd.
    OBJ:
        Generate a box plot of the performance of a single method
    """
    plt.figure(figsize=[5,7])
    sns.boxplot(x = 'Set', y = 'Cor', hue='Traits', data =
                data.loc[data['Method'] == 'rrBLUP'])
    plt.xlabel('Phenotype')
    plt.ylabel('Correlation')
    plt.title('rrBLUP Performance')
    plt.savefig('Correlations/rrBLUP_Sin.png')
    if show:
        plt.show()
    plt.clf()

    plt.figure(figsize=[5,7])
    sns.boxplot(x = 'Set', y = 'Cor', hue='Traits', data =
                data.loc[data['Method'] == 'BayesA'])
    plt.xlabel('Phenotype')
    plt.ylabel('Correlation')
    plt.title('BayesA Performance')
    plt.savefig('Correlations/BayesA_Sin.png')
    if show:
        plt.show()
    plt.clf()

    plt.figure(figsize=[5,7])
    sns.boxplot(x = 'Set', y = 'Cor', hue='Traits', data =
                data.loc[data['Method'] == 'NN'])
    plt.xlabel('Phenotype')
    plt.ylabel('Correlation')
    plt.title('NN Performance')
    plt.savefig('Correlations/NN_Sin.png')
    if show:
        plt.show()
    plt.clf()

if __name__ == '__main__':
    my_filename = 'Correlations/All_Corr.tsv'
    data = load_pd(my_filename)
    seaborn_multi(data, show=False)
    seaborn_singular(data, show=False)
