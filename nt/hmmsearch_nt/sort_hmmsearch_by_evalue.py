import tqdm

if __name__ == '__main__':

    #with open('nt_hmmsearch_default_sorted.tbl') as in_fh:
    with open('nt_hmmsearch_default_sorted.tbl') as in_fh:
        with open('nt_hmmsearch_default_sorted_evalue.tbl', 'w') as out_fh:

            #skip first four lines
            next(in_fh)
            next(in_fh)
            next(in_fh)
            next(in_fh)

            first_line = next(in_fh).strip().split()

            current_read = first_line[0]
            hits = [first_line]

            for line in tqdm.tqdm(in_fh):
                while line.startswith('#'):
                    try:
                        line = next(in_fh)
                    except StopIteration:
                        break



                line = line.strip().split()
                read = line[0]
                if read != current_read:
                    # sort and output reads
                    try:
                        hits = sorted(hits, key=lambda x: float(x[4]))
                    except IndexError:
                        print(hits)
                        assert False
                    for hit in hits:
                        out_fh.write("\t".join(hit) + '\n')
                    # make a new store
                    hits = [line]
                    current_read = read
                else:
                    hits.append(line)

            # for final set of reads
            try:
                hits = sorted(hits, key=lambda x: float(x[4]))
                for hit in hits:
                    out_fh.write("\t".join(hit) + '\n')
            except IndexError:
                print(hits)
