import sys
import pandas as pd


def main():
    df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col = 0)

    if sys.argv[2] == '0':
        df.drop('-', axis=1, inplace=True)

    sum = df.sum(axis=1)

    df_proportion = df.div(sum, axis=0)

    df_proportion.to_csv(sys.argv[3], sep='\t')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python convert_ps_site_count_to_proportion.py <ps_site_count> <remove_zero_column> <outfile>')
        exit()

    main()
