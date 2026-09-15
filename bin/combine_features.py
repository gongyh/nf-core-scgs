import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cov', required=True, help='Coverage matrix TSV')
    parser.add_argument('--kmer4', required=True, help='4-mer frequency TSV')
    parser.add_argument('--kmer5', required=True, help='5-mer frequency TSV')
    parser.add_argument('--seqstats', required=True, help='Sequence stats TSV')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    cov = pd.read_csv(args.cov, sep='\t', index_col=0)
    k4 = pd.read_csv(args.kmer4, sep='\t', index_col=0)
    k5 = pd.read_csv(args.kmer5, sep='\t', index_col=0)
    stats = pd.read_csv(args.seqstats, sep='\t', index_col=0)


    combined = cov.join([k4, k5, stats], how='inner')
    combined.fillna(0, inplace=True)
    combined.to_csv(args.output, sep='\t')

if __name__ == '__main__':
    main()
