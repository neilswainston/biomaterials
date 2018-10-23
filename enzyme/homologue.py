'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

from Bio import Phylo
from cogent import LoadSeqs
from cogent.core.moltype import PROTEIN
from cogent.evolve.models import JTT92
from cogent.phylo import distance, nj
from synbiochem.utils import seq_utils


def get_homologues(id_seqs):
    '''Get homologues.'''
    clustal_file = 'clustal.fasta'
    seq_utils.do_clustal(id_seqs, result_file=clustal_file)
    aln = LoadSeqs(clustal_file, moltype=PROTEIN)
    dists = distance.EstimateDistances(aln, submodel=JTT92())
    dists.run()
    njtree = nj.nj(dists.getPairwiseDistances())
    njtree = njtree.balanced()
    print njtree.asciiArt()

    newick = Phylo.convert(clustal_file, in_format='fasta',
                           out_file='newick.txt', out_format='newick')
    print newick.name


def main(args):
    '''main method.'''
    id_seqs = seq_utils.read_fasta(args[0])
    get_homologues(id_seqs)


if __name__ == '__main__':
    main(sys.argv[1:])
