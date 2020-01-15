#!/usr/bin/env python

import sys
from Bio import SeqIO
import os

if __name__ == '__main__':

    with open(sys.argv[1], 'rU') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            aro = record.id.split('|')[4].replace('ARO:', '')
            with open(os.path.join('card_fasta', aro), 'w') as out_fh:
                SeqIO.write(record, out_fh, 'fasta')

