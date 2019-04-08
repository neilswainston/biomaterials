'''
synbiochem (c) University of Manchester 2019

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_style('whitegrid', {'grid.linestyle': '--'})


def plot(filename, out_dir='out'):
    '''Plot.'''
    xls = pd.ExcelFile(filename)
    name, _ = os.path.splitext(os.path.basename(filename))

    for sheet_name in xls.sheet_names:
        df = _get_df(xls, sheet_name)
        _boxplot(df, os.path.join(out_dir, name))


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


def _boxplot(df, out_dir):
    '''Box plot.'''
    g = sns.catplot(data=df, x='id', y='target conc',
                    hue='substrate concentration', col='target',
                    kind='box', palette='pastel')

    g.set_titles('{col_name}')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    plt.savefig(os.path.join(out_dir, '%s.png' % df.name))


def main(args):
    '''main method.'''
    for filename in args:
        plot(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
