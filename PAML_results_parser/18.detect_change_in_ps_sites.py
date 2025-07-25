import sys
import os
from Bio import SeqIO

def detect_pair_ps_sites(file):
    sequence = []
    with open(file,'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence.append(str(record.seq))

    out = {}
    for i,j in enumerate(sequence[0]):
        lis = []
        for k,l in enumerate(sequence):
            lis.append(l[i])

        ps_change = list(set(lis))
        if len(ps_change) == 2:
            if f'{ps_change[0]}{ps_change[1]}' not in out:
                out[f'{ps_change[0]}{ps_change[1]}'] = 1
            else:
                out[f'{ps_change[0]}{ps_change[1]}'] += 1

    return out


def main():
    pwd = os.getcwd()
    file = os.listdir(sys.argv[1])
    output1 = {}

    os.chdir(sys.argv[1])
    for i in file:
        out = detect_pair_ps_sites(i)
        for j in out:
            if j not in output1:
                output1[j] = out[j]
            else:
                output1[j] += out[j]

    os.chdir(pwd)

    # output processing
    output = {}
    for i in output1:
        if not f'{i[1]}{i[0]}' in output and not i in output:
            output[i] = output1[i]
        else:
            output[f'{i[1]}{i[0]}'] += output1[i]

    with open(sys.argv[2],'w') as f:
        f.write('ps_site_change\tcount\n')
        for i in output:
            f.write(f'{i}\t{output[i]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage')
        exit()

    main()