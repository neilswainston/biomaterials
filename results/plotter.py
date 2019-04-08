'''
synbiochem (c) University of Manchester 2018

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
    fig, axs = plt.subplots(ncols=len(df['target'].unique()),
                            sharey=True)

    for idx, (substrate, group_df) in enumerate(df.groupby('target')):
        ax = axs[idx]
        ax.grid(True)
        ax.set_title(substrate)

        sns.boxplot(x='id', y='target conc',
                    data=group_df,
                    hue='substrate concentration',
                    ax=ax)

    lgd = fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0,
                     title='substrate conc.')

    plt.savefig('%s.png' % df.name,
                bbox_extra_artists=(lgd,),
                bbox_inches='tight')


def main(args):
    '''main method.'''
    for filename in args:
        plot(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
