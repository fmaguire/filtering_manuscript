#!/usr/bin/env python
import json
import collections
import math
import itertools
import numpy as np
import glob
import re
import os
import sys
from Bio import SeqIO

class CARD():
    """
    Parser for the CARD database
    """
    def __init__(self, card_json_fp):
        with open(card_json_fp) as fh:
            with open(card_json_fp) as fh:
                self.card = json.load(fh)
            self.version = self.card['_version']

            # to avoid having to except them later when parsing the other
            # entries
            del self.card['_version']
            del self.card['_timestamp']
            del self.card['_comment']

            self.supported_models = ['protein homolog model',
                                     'protein variant model',
                                     'protein overexpression model',
                                     'protein knockout model']

            self.proteins, self.nucleotides = self.get_sequences()

            self.aro_to_gene_family = self.build_aro_to_gene_family()
            self.gene_family_to_aro = self.build_gene_family_to_aro()
            self.family_sizes = self.calculate_family_sizes()


    def get_sequences(self):
        """
        Gather list of accession, prot and nucleotide sequence tuples from
        the card.json in a pair of dictionaries

        The output format is {ARO: (accession, sequence)}
        """

        data = {'protein_sequence': {}, 'dna_sequence': {}}

        for card_item in self.card.values():
            if card_item['model_type'] in self.supported_models:

                aro = card_item['ARO_accession']
                aro_name = card_item['ARO_name']
                sequences = card_item['model_sequences']['sequence']

                for seq_ix in sequences:
                    for sequence_type in sequences[seq_ix]:

                        # to ignore the taxonomy data
                        if sequence_type in ['protein_sequence', 'dna_sequence']:
                            sequence = sequences[seq_ix][sequence_type]
                            acc = ">gb|{}|{}|{}|".format(sequence['accession'],
                                                         aro,
                                                         aro_name.replace(' ', '_'))
                            data[sequence_type].update({aro: (acc,
                                                              sequence['sequence'])})
        return data['protein_sequence'], data['dna_sequence']


    def build_aro_to_gene_family(self):
        """
        Build a dictionary mapping each ARO to a single AMR gene family
        """
        aro_to_gene_family = {}
        for card_item in self.card.values():
            if card_item['model_type'] not in self.supported_models:
                continue

            aro = card_item['ARO_accession']
            # as multiple gene families are possible per ARO
            # although they are relatively rare so maybe I should talk with
            # Andrew about making them unique
            gene_families = []
            for category in card_item['ARO_category'].values():
                if category['category_aro_class_name'] == 'AMR Gene Family':
                    gene_families.append(category['category_aro_name'])

            # fix the situations where there are multiple gene families manually
            # all glycopeptide resistant gene clusters have 2 gene families one
            # indicating that it is a grgc and the other with its class
            # for now we are just using the cluster level and will deal with
            # the specifics at ARO level.
            if "glycopeptide resistance gene cluster" in gene_families and len(gene_families) > 1:
                gene_families.remove('glycopeptide resistance gene cluster')

            # this is a fusion protein so can be assigned to a new class
            if aro == '3002598':
                gene_families = ["AAC(6')_ANT(3'')"]
            if aro == '3002597':
                gene_families = ["APH(2'')_AAC(6')"]
            if aro in ['3002546', '3002600']:
                gene_families = ["AAC(3)_AAC(6')"]

            # also a fusion so assigned a new fusion class
            if 'class C LRA beta-lactamase' in gene_families and \
                    'class D LRA beta-lactamase' in gene_families:
                gene_families = ['class D/class C beta-lactamase fusion']

            # 23S with multiple resistance classes
            if aro == '3004181':
                gene_families = ['23S rRNA with mutation conferring resistance '
                                 'to macrolide and streptogramins antibiotics']

            # additional self-resistance class to indicate resistance genes made by
            # antibiotic producer removing the self resistant term
            if "fluoroquinolone resistant parC" in gene_families and \
                    "fluoroquinolone self resistant parC" in gene_families:
                gene_families = ['fluoroquinolone resistant parC']

            if "kirromycin self resistant EF-Tu" in gene_families and \
                    'elfamycin resistant EF-Tu' in gene_families:
                gene_families = ['elfamycin resistant EF-Tu']

            if 'aminocoumarin self resistant parY' in gene_families and \
                    'aminocoumarin resistant parY' in gene_families:
                gene_families = ['aminocoumarin resistant parY']

            # efflux components
            if aro in ['3000263', '3000833', '3003382', '3000832',
                           '3000815', '3003896', '3000823', '3003511',
                           '3003381', '3000817', '3003895', '3000676',
                           '3003383', '3003585', '3004107', '3003820']:
                gene_families = ['efflux regulator']
            if aro == '3000237':
                gene_families = ['efflux component']


            # They are homologous parts of topo IV and II but it looks like this is actually parC
            if 'fluoroquinolone resistant parC' in gene_families and \
                    'fluoroquinolone resistant gyrA' in gene_families:
                gene_families = ['fluoroquinolone resistant parC']

            # things that need fixed
            # this looks like a mistake and is only UhpA
            if 'UhpT' in gene_families and 'UhpA' in gene_families:
                gene_families = ['UhpA']

            # missing families
            if aro == "3004450":
                gene_families = ['TRU beta-lactamase']
            if aro == "3004294":
                gene_families = ['BUT beta-lactamase']

            aro_to_gene_family.update({aro: gene_families})

        # check dict is valid and makes unique aro:amr family relationship
        for aro, gene_family in aro_to_gene_family.items():
            if len(gene_family) != 1:
                print(aro, gene_family)
                raise ValueError("ARO and gene families don't map 1:1 {}, {}".format(aro, gene_family))
            else:
                aro_to_gene_family[aro] = gene_family[0]

            if aro not in self.proteins:
                raise ValueError("ARO not in proteins {}".format(aro))
            if aro not in self.nucleotides:
                raise ValueError("ARO not in nucleotides {}".format(aro))

        # clean out AROs not in metagenome
        for aro in ["3003023", "3004253", "3003022", "3003021", "3003820", "3004056", "3002855"]:
            del aro_to_gene_family[aro]


        return aro_to_gene_family

    def build_gene_family_to_aro(self):
        """
        Reverse the dictionary mapping AROs to gene families
        """
        gene_family_to_aro = {}
        for aro, gene_family in self.aro_to_gene_family.items():
            if gene_family not in gene_family_to_aro:
                gene_family_to_aro.update({gene_family: [aro]})
            else:
                gene_family_to_aro[gene_family].append(aro)

        return gene_family_to_aro

    def calculate_family_sizes(self):
        family_sizes = collections.Counter()
        for aro, seq in self.nucleotides.items():
            if aro in self.aro_to_gene_family:
                amr_family = self.aro_to_gene_family[aro]
                family_sizes[amr_family] += 1
        return family_sizes


    def get_seqs_per_family(self, seq_type):
        """
        Write sequences to per family fasta files
        """
        if seq_type == 'protein':
            seq_folder = 'family_fasta_aa'
            seq_dict = self.proteins
        elif seq_type == 'nucleotide':
            seq_folder = 'family_fasta_nt'
            seq_dict = self.nucleotides
        else:
            raise ValueError("Seq Type must be 'protein' or 'nucleotide'")

        if not os.path.exists(seq_folder):
            os.mkdir(seq_folder)

        for aro, seq_data in seq_dict.items():
            acc, seq = seq_data
            amr_family = self.aro_to_gene_family[aro]
            family_fasta_fp = self.convert_amr_family_to_filename(amr_family) + '.fasta'

            with open(seq_folder + "/" + family_fasta_fp, 'a') as out_fh:
                out_fh.write(acc + '\n' + seq + '\n')

    def convert_amr_family_to_filename(self, family):
        fp = os.path.join(family.replace(' ', '_').replace('/', '_'))
        return fp

    def get_amr_family_from_filenames(self):
        family_to_filenames = {}
        for family in self.aro_to_gene_family.values():
            family_to_filenames[family] = self.convert_amr_family_to_filename(family)

        filename_to_family = {}
        for family, filename in family_to_filenames.items():
            filename_to_family[filename] = family
        return filename_to_family


if __name__ == '__main__':

    card = CARD('../../data/CARD_canonical/card.json')
    card.get_seqs_per_family('protein')
    card.get_seqs_per_family('nucleotide')
