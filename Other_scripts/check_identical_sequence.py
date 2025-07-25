import sys
from Bio import SeqIO
import os


def check_identical_sequence(dir):
    path = []
    for root,dirs,files in os.walk(dir):
        for file in files:
            if file.endswith('.fa'):
                path.append(os.path.join(root,file))

    out = []
    count = 0
    for i in path:
        sequence = []
        with open(i, 'r') as f:
            count += 1
            for record in SeqIO.parse(f, 'fasta'):
                sequence.append(str(record.seq))

        if len(set(sequence)) == 1:
            out.append(i.split('/')[-1].split('.')[0])

    return out,count


def main():
    lis = os.listdir(sys.argv[1])
    out = []
    pwd = os.getcwd()
    count = 0

    os.chdir(sys.argv[1])
    pwd1 = os.getcwd()

    for i in lis:
        os.chdir(i)
        output,count1 = check_identical_sequence('msa_nuc')[0],check_identical_sequence('msa_nuc')[1]
        count += count1
        if len(output) != 0:
            for j in output:
                out.append(f'{i}\t{j}\n')
        os.chdir(pwd1)

    os.chdir(pwd)


    with open(sys.argv[2],'w') as f:
        f.write('cluster\tgene_cluster\n')
        for i in out:
            f.write(i)

    print(count)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage')
        exit()

    main()
