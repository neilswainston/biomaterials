'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from xml.dom.minidom import parse

from synbiochem.utils.seq_utils import write_fasta


def main(args):
    '''main method.'''
    with open(args[0]) as fle:
        acs = [line.split('|')[1] for line in fle]

    dom = parse(args[1])

    id_seqs = {hit.attributes['ac'].value:
               hit.getElementsByTagName('matchSeq')[
        0].firstChild.nodeValue.replace('-', '')
        for hit in dom.getElementsByTagName('hit')
        if hit.attributes['ac'].value in acs}

    write_fasta(id_seqs, 'ac_seqs.fasta')


if __name__ == '__main__':
    main(sys.argv[1:])
