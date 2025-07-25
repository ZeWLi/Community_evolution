import os
import sys
from Bio import SeqIO
from collections import Counter


def count_ps_site_species(file):
    sequence = []
    with open(file,'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence.append(str(record.seq))

    out = []

    for i,j in enumerate(sequence[0]):
        lis = []
        for k,l in enumerate(sequence):
            lis.append(l[i])

        out.append(len(set(lis)))

    return file.split('_')[0],out


def main():
    pwd = os.getcwd()
    lis = os.listdir(sys.argv[1])
    os.chdir(sys.argv[1])
    final_output = {}
    avg = {}

    pwd1 = os.getcwd()

    for i in lis:
        output = {}

        os.chdir(i)

        #if not os.path.exists('ps_site_species_count.txt'):
        pwd2 = os.getcwd()
        lis1 = os.listdir('ps_site')
        os.chdir('ps_site')

        for j in lis1:
            output[count_ps_site_species(j)[0]] = count_ps_site_species(j)[1]

        os.chdir(pwd2)

        temp = []

        with open('ps_site_species_count.txt','w') as f:
            for s in output:
                f.write(f'{s}\t{output[s]}\n')
                for q in output[s]:
                    temp.append(q)

        final_output[i] = dict(Counter(temp).most_common())

        avg[i] = sum(temp) / len(temp) if len(temp) != 0 else 'NA'

        # with open('ps_site_species_count.txt','r') as f:
        #     for line in f:
        #         if i not in final_output:
        #             final_output[i] = []
        #
        #         for q in eval(line.split('\t')[1]):
        #             final_output[i].append(q)

        os.chdir(pwd1)

    os.chdir(pwd)

    with open(sys.argv[2],'w') as f:
        f.write('cluster\tps_site_count\taverage\n')
        for i in final_output:
            f.write(f'{i}\t{final_output[i]}\t{avg[i]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage')
        exit()

    main()
