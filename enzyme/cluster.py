'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from synbiochem.utils import seq_utils


def get_positions(id_seqs, ref_seq_id, target_seq_id, positions):
    '''Get aligned positions relative to reference sequence.'''
    clustal_file = 'clustal.fasta'
    seq_utils.do_clustal(id_seqs, result_file=clustal_file)

    cluster_seqs = seq_utils.read_fasta(clustal_file)

    ref_seq = cluster_seqs[ref_seq_id]
    ref_seq_aln, _ = _get_pos_map(ref_seq)

    trgt_seq = cluster_seqs[target_seq_id]
    _, trgt_aln_seq = _get_pos_map(trgt_seq)

    return [[pos,
             ref_seq[ref_seq_aln[pos]],
             trgt_aln_seq[ref_seq_aln[pos]],
             trgt_seq[ref_seq_aln[pos]]]
            for pos in positions]


def _get_pos_map(seq):
    '''Get 0-indexed position sequence index to alignment index map.'''
    seq_aln_pos_map = {}
    aln_seq_pos_map = {}
    seq_idx = 0

    for align_idx, align_chr in enumerate(seq):
        if align_chr != '-':
            seq_aln_pos_map[seq_idx] = align_idx
            aln_seq_pos_map[align_idx] = seq_idx
            seq_idx += 1

    return seq_aln_pos_map, aln_seq_pos_map


def main(args):
    '''main method.'''
    id_seqs = seq_utils.read_fasta(args[0])
    print get_positions(id_seqs, args[1], args[2],
                        [int(pos) for pos in args[3:]])


if __name__ == '__main__':
    main(sys.argv[1:])
