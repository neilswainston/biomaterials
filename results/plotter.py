'''
synbiochem (c) University of Manchester 2019

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot(filename):
    '''Plot.'''
    sns.set(style='ticks', palette='pastel')
    sns.set_style('whitegrid', {'grid.linestyle': '--'})

    xls = pd.ExcelFile(filename)

    for sheet_name in xls.sheet_names:
        df = _get_df(xls, sheet_name)
        _boxplot(df)


def _get_df(xls, sheet_name):
    '''Get df.'''
    df = pd.read_excel(xls, sheet_name)
    rep_cols = [col for col in df.columns if col.startswith('Rep')]
    reps_df = pd.DataFrame([[idx, val]
                            for idx, vals in df[rep_cols].iterrows()
                            for val in vals.values],
                           columns=['idx', 'target conc']).set_index('idx')

    reformat_df = df.drop(rep_cols, axis=1).join(reps_df)
    reformat_df.name = sheet_name

    return reformat_df


def _boxplot(df):
    '''Box plot.'''
    g = sns.catplot(data=df, x='id', y='target conc',
                    hue='substrate concentration', col='target',
                    sharex=False, kind='box')

    g.set_titles('{col_name}')

    plt.savefig('%s.png' % df.name)


def main(args):
    '''main method.'''
    for filename in args:
        plot(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
