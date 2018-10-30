'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from synbiochem.utils import ice_utils


def duplicate(url, username, password, ice_id, mutations):
    '''duplicate.'''
    ice_client = ice_utils.ICEClient(url, username, password)
    ice_entry = ice_client.get_ice_entry(ice_id).copy()
    dna = ice_entry.get_dna()
    cds = [feature for feature in dna['features']
           if feature['typ'] == 'http://purl.obolibrary.org/obo/SO_0000316']

    print dna['seq']

    for aa_pos, mutation in mutations.iteritems():
        dna_pos = cds[0]['start'] - 1 + (aa_pos * 3)
        dna['seq'] = dna['seq'][:dna_pos] + mutation + dna['seq'][dna_pos + 3:]

    print dna['seq']


def main(args):
    '''main.'''
    duplicate(args[0], args[1], args[2], args[3], {0: 'DBK', 1: 'NTN'})


if __name__ == '__main__':
    main(sys.argv[1:])
