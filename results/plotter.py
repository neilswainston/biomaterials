'''
synbiochem (c) University of Manchester 2019

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import glob
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


sns.set_style('whitegrid', {'grid.linestyle': '--'})


def plot(filename, out_dir='out'):
    '''Plot.'''
    xls = pd.ExcelFile(filename)
    project_name, _ = os.path.splitext(os.path.basename(filename))

    for sheet_name in xls.sheet_names:
        df = _get_df(xls, sheet_name)
        project_dir = os.path.join(out_dir, project_name)
        sheet_dir = os.path.join(project_dir, sheet_name)

        if df is not None:
            if sheet_name.startswith('enzyme'):
                df['Sample description'] = \
                    df.apply(_get_enzyme_screen_desc, axis=1)

                for group_id, group_df in df.groupby(['substrate', 'target']):
                    group_df.name = ' to '.join(group_id) + ' enzyme screen'

                    _boxplot(group_df,
                             os.path.join(sheet_dir, '_'.join(group_id)))
            else:
                df['Sample description'] = \
                    df.apply(_get_pathway_screen_desc, axis=1)

                for group_id, group_df in df.groupby(['substrate']):
                    group_df = group_df.sort_values('Analyte order')
                    group_df.name = group_id + ' pathway screen'

                    _barplot(group_df,
                             os.path.join(sheet_dir, group_id))


def _get_df(xls, sheet_name):
    '''Get df.'''
    xls_df = pd.read_excel(xls, sheet_name, dtype={'plasmid id': object,
                                                   'host id': object})

    if xls_df.empty:
        return None

    xls_df.dropna(how='all', inplace=True)
    xls_df['substrate'].fillna('None', inplace=True)
    xls_df['plasmid id'].fillna('None', inplace=True)

    rep_cols = [col for col in xls_df.columns if col.startswith('Rep')]
    reps_df = pd.DataFrame([[idx, val]
                            for idx, vals in xls_df[rep_cols].iterrows()
                            for val in vals.values
                            if pd.notnull(val)],
                           columns=['idx', 'target conc']).set_index('idx')

    df = xls_df.drop(rep_cols, axis=1).join(reps_df)
    df.name = sheet_name

    return df


def _get_enzyme_screen_desc(row):
    '''Get enzyme screen description.'''
    return ('I' if row['induced'] else 'NI') + \
        (', %0.1f mM' % row['substrate concentration']) + \
        (', %0.1f hr' % row['incubation time'])


def _get_pathway_screen_desc(row):
    '''Get pathway screen description.'''
    return '%s, host: %s' % (row['target'], row['host id'])


def _boxplot(df, out_dir):
    '''Box plot.'''
    num_cols = len(df['plasmid id'].unique()) * \
        len(df['Sample description'].unique())

    fig, ax = _get_fig_ax(num_cols)

    g = sns.boxplot(data=df, x='plasmid id', y='target conc',
                    hue='Sample description',
                    palette='pastel', ax=ax)

    g.set_title(df.name)
    ax.set(xlabel='Plasmid(s)', ylabel='Titre (mg/l)')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fig.savefig(os.path.join(out_dir, '%s.png' % df.name))
    plt.close('all')


def _barplot(df, out_dir):
    '''Bar plot.'''
    num_cols = len(df['plasmid id'].unique()) * \
        len(df['Sample description'].unique())

    fig, ax = _get_fig_ax(num_cols)

    g = sns.barplot(data=df, x='plasmid id', y='target conc',
                    hue='Sample description',
                    palette='pastel', ax=ax)

    g.set_title(df.name)
    ax.set(xlabel='Plasmid(s)', ylabel='Titre (mg/l)')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fig.savefig(os.path.join(out_dir, '%s.png' % df.name))
    plt.close('all')


def _get_fig_ax(num_cols):
    '''Get fig and ax based on number of columns.'''
    fig, ax = plt.subplots()
    fig.set_size_inches(0.6 * num_cols + 2.0, 5.0)

    return fig, ax


def main(args):
    '''main method.'''
    for filename in glob.glob(os.path.join(args[0], '*.xlsx')):
        if '~' not in filename:
            plot(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
