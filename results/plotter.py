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
    xls = pd.ExcelFile(filename)

    for sheet_name in xls.sheet_names:
        _, ax = plt.subplots()

        df = pd.read_excel(xls, sheet_name)
        rep_cols = [col for col in df.columns if col.startswith('Rep')]
        reps_df = pd.DataFrame([[idx, val]
                                for idx, vals in df[rep_cols].iterrows()
                                for val in vals.values],
                               columns=['idx', 'target conc']).set_index('idx')

        reformat_df = df.drop(rep_cols, axis=1).join(reps_df)

        sns.boxplot(x='id', y='target conc',
                    data=reformat_df,
                    hue='substrate concentration',
                    ax=ax)

        # sns.swarmplot(x='id', y='target conc',
        #             data=reformat_df,
        #              hue='substrate concentration',
        #              color='.25',
        #              ax=ax)

        plt.savefig(filename + '.png')


def main(args):
    '''main method.'''
    for filename in args:
        plot(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
