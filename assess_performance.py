#!/usr/bin/env python

import pickle
import json
import pysam
import os
import gzip
import pandas as pd

import utils

from collections import Counter

def check_aro(aro):
    if aro[0] != 3:
        int(aro)
    else:
        raise ValueError(aro)

class FilterPerformance():
    """
    Performance for filter only tools that don't predict specific AROs
    """
    def __init__(self, tool_name, param):
        self.true_positive = 0
        self.false_positive = 0

        self.tool_name = tool_name
        self.params = param

    def to_pandas(self, aro_counts):
        data = {}
        data['false_positive'] = self.false_positive
        data['true_positive'] = self.true_positive
        data['tool'] = self.tool_name
        data['params'] = self.params
        data['total'] = pd.Series(aro_counts).sum()

        df = pd.Series(data)
        df.to_csv('tool_runs/' + self.tool_name + "_" + self.params + ".csv")


class familyPerformance():
    def __init__(self, tool_name, param):

        self.correct_family = {}
        self.wrong_family = {}
        self.false_positive = {}

        self.tool_name = tool_name
        self.params = param

    def add_correct_family(self, family):
        if family not in self.correct_family:
            self.correct_family[family] = 1
        else:
            self.correct_family[family] += 1

    def add_wrong_family(self, family):
        if family not in self.wrong_family:
            self.wrong_family[family] = 1
        else:
            self.wrong_family[family] += 1

    def add_false_positive(self, family):
        if family not in self.false_positive:
            self.false_positive[family] = 1
        else:
            self.false_positive[family] += 1

    def to_pandas(self, family_counts):

        fp = pd.Series(self.false_positive, name='false_positive')
        correct = pd.Series(self.correct_family, name='correct_family')
        wrong = pd.Series(self.wrong_family, name='wrong_family')

        df = pd.concat([correct, wrong, fp],
                axis=1, join='outer')

        df = df.fillna(0)
        df['tool'] = self.tool_name
        df['params'] = self.params

        df['totals'] = pd.Series(family_counts)

        df['missed'] = df['totals'] - (df['correct_family'] \
                + df['wrong_family'])

        df.to_csv('tool_runs/' + self.tool_name + "_" + self.params + ".csv")



class Performance():
    def __init__(self, tool_name, param):
        """
        Got a hit:
            - hit has label
                - label is right
                - label family is right but aro is wrong
                - label family and aro are wrong
            - hit doesn't have label
                - false positive
        No hit:
            - hit has label
                - false negative/miss tally from total correct hits - total hits
            - hit has no label
                - true negative

        """
        self.correct_aro = {}
        self.false_positive = {}
        self.wrong_aro_correct_family = {}
        self.wrong_aro_wrong_family = {}

        self.tool_name = tool_name
        self.params = param

    def add_false_positive(self, aro):
        """
        i.e. a hit that isn't any card sequence
        # allows family for hmmsearch
        """
        if aro not in self.false_positive:
            self.false_positive[aro] = 1
        else:
            self.false_positive[aro] += 1

    def add_wrong_aro_correct_family(self, aro):
        """
        wrong aro but correct family
        """
        if aro not in self.wrong_aro_correct_family:
            self.wrong_aro_correct_family[aro] = 1
        else:
            self.wrong_aro_correct_family[aro] += 1

    def add_correct_aro(self, aro):
        """
        ARO predicted matches label
        """
        if aro not in self.correct_aro:
            self.correct_aro[aro] = 1
        else:
            self.correct_aro[aro] += 1

    def add_wrong_aro_wrong_family(self, aro):
        """
        Hits an ARO but wrong aro from the wrong family
        """
        if aro not in self.wrong_aro_wrong_family:
            self.wrong_aro_wrong_family[aro] = 1
        else:
            self.wrong_aro_wrong_family[aro] += 1

    def to_pandas(self, aro_counts):

        fp = pd.Series(self.false_positive, name='false_positive')
        correct = pd.Series(self.correct_aro, name='correct_aro')
        correct_family = pd.Series(self.wrong_aro_correct_family, name='wrong_aro_correct_family')
        wrong = pd.Series(self.wrong_aro_wrong_family, name='wrong_aro_wrong_family')

        df = pd.concat([correct, correct_family, wrong, fp],
                axis=1, join='outer')
        df = df.fillna(0)
        df['tool'] = self.tool_name
        df['params'] = self.params

        df['totals'] = pd.Series(aro_counts)

        df['missed'] = df['totals'] - (df['correct_aro'] \
                + df['wrong_aro_correct_family'] \
                + df['wrong_aro_wrong_family'])

        df.to_csv('tool_runs/' + self.tool_name + "_" + self.params + ".csv")

class Labels():
    """
    Labelling information and contextualisation
    """

    def __init__(self, label_fp, card_json_fp):

        self.card = utils.CARD(card_json_fp)

        self.labels = self.get_labels(label_fp)

        self.valid_aros = set(self.labels.values())

        self.aro_counts = self.get_aro_counts_in_labels()
        self.aro_to_family = self.card.aro_to_gene_family

        self.valid_families = set(self.aro_to_family.values())

        all_labels = set(self.labels.values())
        for aro in all_labels:
            if aro not in self.aro_to_family:
                print(aro)
                assert False

        for aro in self.aro_to_family:
            if aro not in all_labels:
                print(aro)
                assert False

        self.family_counts = self.get_family_counts_in_labels()

    def get_labels(self, label_fp):
        """
        Get label TSV as dictionary
        """
        if os.path.exists('labels.pkl'):
            with open('labels.pkl', 'rb') as fh:
                labels = pickle.load(fh)
        else:
            labels = {}
            with open(label_fp) as fh:
                for line in fh:
                    line = line.split('\t')
                    # to deal with aac multilabels
                    # not a big deal as we have canonical/prevalence
                    # so nothing is getting under-represented
                    label_aro = line[2].split('; ')[0]

                    # remove private models not in dataset
                    if label_aro in ['3000489', '3000309', '3002818']:
                        pass
                    else:
                        labels.update({line[0]: label_aro})

            with open('labels.pkl', 'wb') as fh:
                pickle.dump(labels, fh)
        return labels


    def get_aro_counts_in_labels(self):
        """
        Get prevalence of each ARO label
        """
        if os.path.exists('aro_counts.pkl'):
            with open('aro_counts.pkl', 'rb') as fh:
                aro_counts = pickle.load(fh)
        else:
            aro_counts = dict(Counter(self.labels.values()))

            # 3003893 and 3003900 have duplicate sequences in CARD so counts
            #
            # grep "^>" card_nucleotides.fna | sort | uniq -d
            #>gb|HG738867|-|2786398-2788945|ARO:3003900|Escherichia_coli_CyaA_with_mutation_conferring_resistance_to_fosfomycin_[Escherichia_coli_str._K-12_substr._MC4100]
            # >gb|HG738867|+|2930707-2931298|ARO:3003893|Escherichia_coli_UhpA_with_mutation_conferring_resistance_to_fosfomycin_[Escherichia_coli_str._K-12_substr._MC4100]
            # need adjusted to get the correct number and not end up with
            # negative missed
            # they aren't next to one another so shouldn't be a major problem
            #aro_counts['3003900'] = 12140
            #aro_counts['3003893'] = 40

            with open('aro_counts.pkl', 'wb') as fh:
                pickle.dump(aro_counts, fh)
        return aro_counts

    def get_family_counts_in_labels(self):
        family_counts = Counter()
        for label in self.labels.values():
            family_counts[self.aro_to_family[label]] += 1
        return family_counts


# aro level performance calculated based on first hit for multimaps as all 30 mapq
# family level performance only
def assess_groot(fp, labels, tool_name, params):

    default_perf = Performance(tool_name, params + "_default")

    samfile = pysam.AlignmentFile(fp, 'rb')
    seen = set()
    for record in samfile:

        read = record.qname
        if read in seen:
            continue
        else:
            if record.reference_id != -1:

                hit_aro = record.reference_name.split('|')[4].replace('ARO:', '')
                check_aro(hit_aro)
                if hit_aro not in labels.valid_aros:
                    continue

                family_hit = labels.aro_to_family[hit_aro]

                if read in labels.labels:
                    true_aro = labels.labels[read]
                    true_family = labels.aro_to_family[true_aro]

                    if hit_aro == true_aro:
                        default_perf.add_correct_aro(true_aro)
                    elif hit_aro != true_aro and family_hit == true_family:
                        default_perf.add_wrong_aro_correct_family(true_aro)
                    else:
                        default_perf.add_wrong_aro_wrong_family(true_aro)

                else:
                    # if read gets hit and not in labels
                    default_perf.add_false_positive(hit_aro)

            seen.add(read)

    default_perf.to_pandas(labels.aro_counts)


def assess_blast_table(fp, labels, tool_name, params, protein=False, nt_db=False):

    default_perf = Performance(tool_name, params + "_default")

    min_1e5_perf = Performance(tool_name, params + "_1e-5")
    min_1e10_perf = Performance(tool_name, params + "_1e-10")

    min_50_perf = Performance(tool_name, params + "_min50")
    min_100_perf = Performance(tool_name, params + "_min100")

    fh = gzip.open(fp, 'rb')

    seen = set()
    for line in fh:
        line = line.decode('utf-8').split('\t')

        if protein:
            read = "_".join(line[0].split('_')[:-3])
        else:
            read = line[0]

        # diamond will output multiple hits to the same sequence even if
        # constraining to 1 top hit therefore if the previous read is the same
        # then skip to the next read
        if read in seen:
            continue
        else:
            if nt_db:
                aro = line[1].split('|')[4].replace('ARO:', '')
            else:
                aro = line[1].split('|')[2].replace('ARO:', '')

            check_aro(aro)
            if aro not in labels.valid_aros:
                continue


            evalue = float(line[10])
            bitscore = float(line[11])

            if read in labels.labels:

                true_aro = labels.labels[read]

                # aro matches label
                if aro == true_aro:
                    default_perf.add_correct_aro(true_aro)

                    if evalue <= 1e-5:
                        min_1e5_perf.add_correct_aro(true_aro)

                    if evalue <= 1e-10:
                        min_1e10_perf.add_correct_aro(true_aro)

                    if bitscore >= 50:
                        min_50_perf.add_correct_aro(true_aro)

                    if bitscore >= 100:
                        min_100_perf.add_correct_aro(true_aro)

                # wrong aro but correct family
                elif labels.aro_to_family[aro] == labels.aro_to_family[true_aro]:

                    default_perf.add_wrong_aro_correct_family(true_aro)

                    if evalue <= 1e-5:
                        min_1e5_perf.add_wrong_aro_correct_family(true_aro)

                    if evalue <= 1e-10:
                        min_1e10_perf.add_wrong_aro_correct_family(true_aro)

                    if bitscore >= 50:
                        min_50_perf.add_wrong_aro_correct_family(true_aro)

                    if bitscore >= 100:
                        min_100_perf.add_wrong_aro_correct_family(true_aro)

                # wrong aro and wrong family
                else:
                    default_perf.add_wrong_aro_wrong_family(true_aro)
                    if evalue <= 1e-5:
                        min_1e5_perf.add_wrong_aro_wrong_family(true_aro)

                    if evalue <= 1e-10:
                        min_1e10_perf.add_wrong_aro_correct_family(true_aro)

                    if bitscore >= 50:
                        min_50_perf.add_wrong_aro_correct_family(true_aro)

                    if bitscore >= 100:
                        min_100_perf.add_wrong_aro_correct_family(true_aro)


            else:
                # if read gets hit and not in labels
                default_perf.add_false_positive(aro)

                if evalue <= 1e-5:
                    min_1e5_perf.add_false_positive(aro)

                if evalue <= 1e-10:
                    min_1e10_perf.add_false_positive(aro)

                if bitscore >= 50:
                    min_50_perf.add_false_positive(aro)

                if bitscore >= 100:
                    min_100_perf.add_false_positive(aro)

            # add read to 'seen'
            seen.add(read)


    fh.close()
    default_perf.to_pandas(labels.aro_counts)
    min_1e5_perf.to_pandas(labels.aro_counts)
    min_1e10_perf.to_pandas(labels.aro_counts)
    min_50_perf.to_pandas(labels.aro_counts)
    min_100_perf.to_pandas(labels.aro_counts)


def assess_sam(fp, labels, tool_name, params, protein=False):
    """
    Might be unfair to sam mappers as only taking the first hit
    """

    default_perf = Performance(tool_name, params + "_default")

    samfile = pysam.AlignmentFile(fp, 'rb', header=None)
    seen = set()
    for record in samfile:
        if record.reference_id != -1:
            if protein:
                read = ':'.join(record.qname.split(':')[3:])
            else:
                read = record.qname

            if read in seen:
                continue
            else:

                if protein:
                    hit_aro = record.reference_name.split('|')[2].replace('ARO:', '')
                else:
                    hit_aro = record.reference_name.split('|')[4].replace('ARO:', '')

                check_aro(hit_aro)
                if hit_aro not in labels.valid_aros:
                    continue

                if read in labels.labels:
                    true_aro = labels.labels[read]

                    hit_family = labels.aro_to_family[hit_aro]
                    true_family = labels.aro_to_family[true_aro]

                    if hit_aro == true_aro:
                        # aro matches label
                        default_perf.add_correct_aro(true_aro)
                    elif hit_aro != true_aro and hit_family == true_family:
                        # wrong aro but correct family
                        default_perf.add_wrong_aro_correct_family(true_aro)
                    else:
                        # wrong aro and wrong family
                        default_perf.add_wrong_aro_wrong_family(true_aro)

                else:
                    # if read gets hit and not in labels we want to add one to
                    # number of false hits against that ARO therefore use hit_aro
                    # and not true_aro
                    #if hit_aro == '3003917':
                    #    for key in dir(record):
                    #        print("\n####")
                    #        print(key, ":" , str(getattr(record, key)))
                    default_perf.add_false_positive(hit_aro)
                seen.add(read)

    default_perf.to_pandas(labels.aro_counts)


def assess_hmmsearch(fp, labels, tool_name, params, protein=False):
    """
    """

    default_perf = familyPerformance(tool_name, params + "_default")

    min_1e5_perf =  familyPerformance(tool_name, params + "_1e-5")
    min_1e10_perf = familyPerformance(tool_name, params + "_1e-10")
    min_1e20_perf = familyPerformance(tool_name, params + "_1e-20")
    min_1e40_perf = familyPerformance(tool_name, params + "_1e-40")


    min_50_perf =  familyPerformance(tool_name, params + "_min50")
    min_100_perf = familyPerformance(tool_name, params + "_min100")
    min_150_perf = familyPerformance(tool_name, params + "_min150")
    min_200_perf = familyPerformance(tool_name, params + "_min200")

    fh = open(fp, 'r')

    # sort the tbl file so its read ordered not family ordered
    seen = set()

    filename_to_family = labels.card.get_amr_family_from_filenames()



    for line in fh:
        if line.startswith('#'):
            continue
        else:

            line = line.strip().split()

            if protein:
                read = "_".join(line[0].split('_')[:-3])
            else:
                read = line[0]

            # hmmsearch will output multiple hits to the same sequence even if
            # constraining to 1 top hit therefore if the previous read is the same
            # then skip to the next read
            if read in seen:
                continue
            else:
                # skip the vanU extra
                if line[2] in filename_to_family:
                    family_hit = filename_to_family[line[2]]
                else:
                    continue

                if family_hit not in labels.valid_families:
                    print(family_hit)
                    assert False

                evalue = float(line[4])
                bitscore = float(line[5])


                if read in labels.labels:
                    true_aro = labels.labels[read]

                    # wrong aro but correct family
                    true_family = labels.aro_to_family[true_aro]
                    if true_family == family_hit:

                        default_perf.add_correct_family(true_family)

                        if evalue <= 1e-5:
                            min_1e5_perf.add_correct_family(true_family)

                        if evalue <= 1e-10:
                            min_1e10_perf.add_correct_family(true_family)

                        if bitscore >= 50:
                            min_50_perf.add_correct_family(true_family)

                        if bitscore >= 100:
                            min_100_perf.add_correct_family(true_family)

                    # wrong aro and wrong family
                    else:
                        default_perf.add_wrong_family(true_family)
                        if evalue <= 1e-5:
                            min_1e5_perf.add_wrong_family(true_family)

                        if evalue <= 1e-10:
                            min_1e10_perf.add_wrong_family(true_family)

                        if bitscore >= 50:
                            min_50_perf.add_wrong_family(true_family)

                        if bitscore >= 100:
                            min_100_perf.add_wrong_family(true_family)


                else:
                    # if read gets hit and not in labels
                    # add family name
                    default_perf.add_false_positive(family_hit)

                    if evalue <= 1e-5:
                        min_1e5_perf.add_false_positive(family_hit)

                    if evalue <= 1e-10:
                        min_1e10_perf.add_false_positive(family_hit)

                    if bitscore >= 50:
                        min_50_perf.add_false_positive(family_hit)

                    if bitscore >= 100:
                        min_100_perf.add_false_positive(family_hit)

                # add read to 'seen'
                seen.add(read)


    fh.close()
    default_perf.to_pandas(labels.family_counts)
    min_1e5_perf.to_pandas(labels.family_counts)
    min_1e10_perf.to_pandas(labels.family_counts)
    min_50_perf.to_pandas(labels.family_counts)
    min_100_perf.to_pandas(labels.family_counts)



def assess_biobloom(fp, labels, tool_name, params):

    default_perf = FilterPerformance(tool_name, params + "_default")
    any_perf = FilterPerformance(tool_name, params + '_any')
    high_perf = FilterPerformance(tool_name, params + '_min0.5')
    ver_high_perf = FilterPerformance(tool_name, params + '_min0.8')

    fh = gzip.open(fp, 'rb')
    seen = set()
    for line in fh:
        line = line.decode('utf-8')
        if line.startswith('>'):

            read = ' '.join(line.replace('>', '').split(' ')[:-1])
            if read in seen:
                continue
            else:


                score = float(line.split(' ')[-1])

                if read in labels.labels:
                    #default min is 0.15
                    any_perf.true_positive += 1

                    if score >= 0.15:
                        default_perf.true_positive += 1

                    if score >= 0.5:
                        high_perf.true_positive += 1

                    if score >= 0.8:
                        ver_high_perf.true_positive += 1

                else:
                    # if read gets hit and not in labels
                    default_perf.false_positive += 1
                    high_perf.false_positive += 1
                    any_perf.false_positive += 1
                    ver_high_perf.false_positive += 1

                seen.add(read)

    fh.close()
    any_perf.to_pandas(labels.aro_counts)
    default_perf.to_pandas(labels.aro_counts)
    high_perf.to_pandas(labels.aro_counts)
    ver_high_perf.to_pandas(labels.aro_counts)


if __name__ == '__main__':

    # read labels as dict
    print('Loading labels')
    labels = Labels('../data/test_metagenome/test_clean_labels.tsv',
                    '../data/CARD_canonical/card.json')

    #### nucleotide query aa database methods
    print('nt domain, aa co-domain')
    blastx_default = assess_blast_table('nt_to_aa/blastx/default.top_hits.gz', labels, 'blastx', 'default')

    diamond_blastx = assess_blast_table('nt_to_aa/diamond_blastx/default.top_hits.gz', labels, 'diamond_blastx', 'default')
    diamond_sensitive_blastx = assess_blast_table('nt_to_aa/diamond_blastx/sensitive.top_hits.gz', labels, 'diamond_blastx', 'sensitive')
    diamond_more_sensitive_blastx = assess_blast_table('nt_to_aa/diamond_blastx/more_sensitive.top_hits.gz', labels, 'diamond_blastx', 'more_sensitive')

    paladin_min5 = assess_sam('nt_to_aa/paladin/min5_clean.bam', labels, 'paladin', 'min5', protein=True)
    paladin_min15 = assess_sam('nt_to_aa/paladin/min15_clean.bam', labels, 'paladin', 'default', protein=True)
    paladin_min25 = assess_sam('nt_to_aa/paladin/min25_clean.bam', labels, 'paladin', 'min25', protein=True)

    ##### amino acid space methods
    print('aa domain, aa co-domain')
    blastp_default = assess_blast_table('aa/orfM/blastp/default.top_hits.gz', labels, 'blastp', 'default', protein=True)

    diamond_blastp = assess_blast_table('aa/orfM/diamond_blastp/default.top_hits.gz', labels, 'diamond_blastp', 'default', protein=True)
    diamond_sensitive_blastp = assess_blast_table('aa/orfM/diamond_blastp/sensitive.top_hits.gz', labels, 'diamond_blastp', 'sensitive', protein=True)
    diamond_more_sensitive_blastp = assess_blast_table('aa/orfM/diamond_blastp/more_sensitive.top_hits.gz', labels, 'diamond_blastp', 'more_sensitive', protein=True)

    hmmsearch_aa = assess_hmmsearch('aa/orfM/hmmsearch_aa/aa_hmmsearch_default_sorted_evalue.tbl',
                                    labels,
                                    'hmmsearch_aa',
                                    'default', protein=True)

    ####### nucleotide space methods
    print('nt domain, nt co-domain')
    megablast_default = assess_blast_table('nt/blastn/megablast_default.top_hits.gz', labels, 'blastn', 'megablast', nt_db=True)
    #blastn_default = assess_blast_table('nt/blastn/blastn_default.out6.gz', labels, 'blastn', 'default', nt_db=True)

    bwa_default = assess_sam('nt/bwa/default_clean.bam', labels, 'bwa', 'default')
    bwa_high = assess_sam('nt/bwa/high_thresh_clean.bam', labels, 'bwa', 'high')
    bwa_low = assess_sam('nt/bwa/low_thresh_clean.bam', labels, 'bwa', 'low')

    bowtie2 = assess_sam('nt/bowtie2/default_clean.bam',
                         labels, 'bowtie2', 'default')
    bowtie2 = assess_sam('nt/bowtie2/e2e_very_fast_clean.bam',
                         labels, 'bowtie2', 'e2e_very_fast')
    bowtie2 = assess_sam('nt/bowtie2/e2e_fast_clean.bam',
                         labels, 'bowtie2', 'e2e_fast')
    bowtie2 = assess_sam('nt/bowtie2/e2e_very_sensitive_clean.bam',
                         labels, 'bowtie2', 'e2e_very_sensitive')
    bowtie2 = assess_sam('nt/bowtie2/l_very_fast_clean.bam',
                         labels, 'bowtie2', 'l_very_fast')
    bowtie2 = assess_sam('nt/bowtie2/l_fast_clean.bam',
                         labels, 'bowtie2', 'l_fast')
    bowtie2 = assess_sam('nt/bowtie2/l_sensitive_clean.bam',
                         labels, 'bowtie2', 'l_sensitive')
    bowtie2 = assess_sam('nt/bowtie2/l_very_sensitive_clean.bam',
                         labels, 'bowtie2', 'l_very_sensitive')

    ###biobloom
    biobloom = assess_biobloom('nt/biobloom/output/5_card_5.fa.gz', labels, 'biobloom', 'k5')
    biobloom = assess_biobloom('nt/biobloom/output/9_card_9.fa.gz', labels, 'biobloom', 'k9')
    biobloom = assess_biobloom('nt/biobloom/output/15_card_15.fa.gz', labels, 'biobloom', 'k15')
    biobloom = assess_biobloom('nt/biobloom/output/25_card_25.fa.gz', labels, 'biobloom', 'k25')



    groot = assess_groot('nt/groot/align/card_k5_s128_j090_default_clean.bam', labels, 'groot',   'card_k5_s128_j090')
    groot = assess_groot('nt/groot/align/card_k5_s128_j099_default_clean.bam', labels, 'groot',   'card_k5_s128_j090')
    groot = assess_groot('nt/groot/align/card_k5_s256_j090_default_clean.bam', labels, 'groot',   'card_k5_s256_j090')
    groot = assess_groot('nt/groot/align/card_k5_s256_j099_default_clean.bam', labels, 'groot',   'card_k5_s256_j090')
    groot = assess_groot('nt/groot/align/card_k5_s64_j090_default_clean.bam', labels, 'groot',   'card_k5_s64_j090')
    groot = assess_groot('nt/groot/align/card_k5_s64_j099_default_clean.bam', labels, 'groot',   'card_k5_s64_j099')
    groot = assess_groot('nt/groot/align/card_k7_s128_j099_default_clean.bam', labels, 'groot',   'card_k7_s128_j090')
    groot = assess_groot('nt/groot/align/card_k7_s256_j099_default_clean.bam', labels, 'groot',   'card_k7_s256_j090')
    groot = assess_groot('nt/groot/align/card_k7_s64_j099_default_clean.bam', labels, 'groot',   'card_k7_s64_j0990')
    groot = assess_groot('nt/groot/align/card_k9_s128_j090_default_clean.bam', labels, 'groot',   'card_k9_s128_j090')
    groot = assess_groot('nt/groot/align/card_k9_s128_j099_default_clean.bam', labels, 'groot',   'card_k9_s128_j090')
    groot = assess_groot('nt/groot/align/card_k9_s256_j090_default_clean.bam', labels, 'groot',   'card_k9_s256_j090')
    groot = assess_groot('nt/groot/align/card_k9_s256_j099_default_clean.bam', labels, 'groot',   'card_k9_s256_j090')
    groot = assess_groot('nt/groot/align/card_k9_s64_j090_default_clean.bam', labels, 'groot',   'card_k9_s64_j0900')
    groot = assess_groot('nt/groot/align/card_k9_s64_j099_default_clean.bam', labels, 'groot',   'card_k9_s64_j0990')


    # hmmsearch
    hmmsearch_nt = assess_hmmsearch('nt/hmmsearch_nt/nt_hmmsearch_default_sorted_evalue.tbl',
                                    labels,
                                    'hmmsearch_nt',
                                    'default')

