# goal: read in prefetch file(s); recalc ANI if necessary
import sys
import argparse
import csv
#import pandas as pd
import numpy as np
from sourmash.logging import notify
from sourmash.distance_utils import containment_to_distance


def main(args):

    # handle file input
    prefetch_csvs = args.prefetch_csvs
    if args.from_file:
        for inF in args.from_file:
            ff_csvs = [x.strip() for x in open(inF, 'r')]
            prefetch_csvs += ff_csvs

    # read in each file and load into table. there should only be a single line in each
    writer=None
    with open(args.output_csv, 'w') as outF:
        for inF in prefetch_csvs:
            with open(inF, 'r') as pf:
                pf_r = csv.DictReader(pf)
                if writer is None:
                    writer = csv.DictWriter(outF, fieldnames = pf_r.fieldnames)
                    writer.writeheader()

                for row in pf_r:
                    if not row["query_containment_ani"]:
                        # grab required columns
                        q_containment = float(row['f_match_query'])
                        m_containment = float(row['f_query_match'])
                        ksize = int(row['ksize'])
                        scaled = int(row['scaled'])
                        n_unique_kmers = int(row['query_bp'])
                        # recalculate containment ani
                        query_ani_res = containment_to_distance(q_containment, ksize, scaled, n_unique_kmers=n_unique_kmers)
                        match_ani_res = containment_to_distance(m_containment, ksize, scaled, n_unique_kmers=n_unique_kmers)
                        # don't let any ANI values get zeroed out --> estimate independtly
                        query_ani = 1-query_ani_res.dist
                        match_ani = 1-match_ani_res.dist
                        avg_ani = np.mean([query_ani, match_ani])
                        max_ani = max(query_ani, match_ani)

                        row["query_containment_ani"] = query_ani
                        row["match_containment_ani"] = match_ani
                        row["average_containment_ani"] = avg_ani
                        row["max_containment_ani"] = max_ani

                        writer.writerow(row)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('prefetch_csvs', nargs='*')
    p.add_argument('--from-file', '--prefetch-from-file', nargs="*", help="file(s) containing paths to prefetch csvs")
    p.add_argument('-o', '--output-csv', required=True, help='output csv')
    args = p.parse_args()
    sys.exit(main(args))
