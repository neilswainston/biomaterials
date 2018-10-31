'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from xml.dom.minidom import parse

from synbiochem.utils.seq_utils import write_fasta


def parse_blast(ids_filename, blast_filename, out_filename='ac_seqs.fasta'):
    '''Parse blast XML.'''
    with open(ids_filename) as fle:
        acs = [line.split('|')[1] for line in fle]

    dom = parse(blast_filename)

    id_seqs = {hit.attributes['ac'].value:
               hit.getElementsByTagName('matchSeq')[
        0].firstChild.nodeValue.replace('-', '')
        for hit in dom.getElementsByTagName('hit')
        if hit.attributes['ac'].value in acs}

    write_fasta(id_seqs, out_filename)


def main(args):
    '''main method.'''
    parse(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
