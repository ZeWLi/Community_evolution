import pandas as pd
import sys


def main():
    if len(sys.argv) == 5:
        p_value = pd.read_table(sys.argv[1],sep='\t',engine='python')
        omega = pd.read_table(sys.argv[2],sep='\t',engine='python')
        anno = pd.read_table(sys.argv[3],sep='\t',engine='python')

        out1 = pd.merge(p_value,omega,on='cluster')
        out = pd.merge(out1,anno,on='cluster')

        out.to_csv(sys.argv[4],sep = '\t',index=False)
    else:
        p_value = pd.read_table(sys.argv[1], sep='\t',engine='python')
        omega = pd.read_table(sys.argv[2], sep='\t',engine='python')

        out = pd.merge(p_value, omega, on='cluster')
        out.to_csv(sys.argv[3], sep='\t', index=False)


if __name__ == '__main__':
    if len(sys.argv) != 5 and len(sys.argv) != 4:
        print('Usage: python 15.py p.txt omega.txt (anno.txt) out.txt')
        exit()

    main()
