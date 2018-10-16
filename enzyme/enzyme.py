'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from xml.dom.minidom import parse


def main(args):
    '''main method.'''
    with open(args[0]) as fle:
        acs = [line.split('|')[1] for line in fle]

    dom = parse(args[1])

    hits = [[hit.attributes['ac'].value,
             hit.getElementsByTagName('querySeq')[0].firstChild.nodeValue,
             hit.getElementsByTagName('pattern')[0].firstChild.nodeValue,
             hit.getElementsByTagName('matchSeq')[0].firstChild.nodeValue]
            for hit in dom.getElementsByTagName('hit')
            if hit.attributes['ac'].value in acs]

    for hit in hits:
        print '\n'.join(hit[1:]) + '\n'


if __name__ == '__main__':
    main(sys.argv[1:])
