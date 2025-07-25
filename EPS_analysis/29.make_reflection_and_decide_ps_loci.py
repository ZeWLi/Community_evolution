import os
import sys
from Bio import SeqIO


def parse_phobius_line(file, dir_msa):
    with open(file, 'r') as f:
        lines = f.readlines()[1:]

    sequence = {}
    location = {}

    for line in lines:
        split = line.rstrip('\n').split('\t')
        cluster, seq = split[0].split(',')
        os.system(f"seqkit grep --pattern '{seq}' {dir_msa}/{cluster}.mafft.fa > temp.fa")
        with open('temp.fa', 'r') as f:
            sequence[cluster] = SeqIO.read(f, 'fasta')

        os.system('rm temp.fa')

        if len(split[1]) != 0:
            for i in split[1].split(',')[:-1]:
                interval = i.split('-')
                if cluster not in location:
                    location[cluster] = {}

                location[cluster][(eval(interval[0]), eval(interval[1]))] = 'inside'

        if len(split[2]) != 0:
            for i in split[2].split(',')[:-1]:
                interval = i.split('-')
                if cluster not in location:
                    location[cluster] = {}
                location[cluster][(eval(interval[0]), eval(interval[1]))] = 'outside'

        if len(split[3]) != 0:
            for i in split[3].split(',')[:-1]:
                interval = i.split('-')
                if cluster not in location:
                    location[cluster] = {}

                location[cluster][(eval(interval[0]), eval(interval[1]))] = 'transmembrane'

        if len(split[4]) != 0:
            for i in split[4].split(',')[:-1]:
                start, end = map(int, i.split('-'))
                location.setdefault(cluster, {})[(start, end)] = 'signal'


    return sequence, location


def main():
    lis = os.listdir(sys.argv[1])
    pwd = os.getcwd()
    os.chdir(sys.argv[1])
    count = {}

    pwd1 = os.getcwd()
    for i in lis:
        os.chdir(i)
        out = []
        # sequence, location = parse_phobius_line('phobius.parse.out', 'msa')
        sequence, location = parse_phobius_line('phobius/phobius.reparse.out', 'msa')

        with open('ps_site_summary.txt','r') as f:
            lines = f.readlines()[1:]

        inside, outside, trans, signal = 0, 0, 0, 0

        for line in lines:
            split = line.split('\t')
            cluster = split[0]
            loci = split[1]
            seq = sequence[cluster].seq

            gap_count = seq[0:eval(loci)].count("-")
            modified_loci = eval(loci) - gap_count

            for interval, attribute in location[cluster].items():
                if interval[0] <= modified_loci <= interval[1]:
                    out.append(f'{cluster}\t{loci}\t{attribute}\n')
                    if attribute == 'inside':
                        inside += 1
                    elif attribute == 'outside':
                        outside += 1
                    elif attribute == 'transmembrane':
                        trans += 1
                    elif attribute == 'signal':
                        signal += 1
                    break

        count[i] = [inside, outside, trans, signal]

        #with open('ps_loci_location.txt', 'w') as f:
            #f.write('cluster\tloci\tlocation\n')
            #for s in out:
                #f.write(s)

        os.chdir(pwd1)

    os.chdir(pwd)

    with open(sys.argv[2], 'w') as f:
        #f.write('cluster\tinside\toutside\ttransmembrane\tsignal\n')
        f.write('cluster\tintracellular\textracellular\ttransmembrane\tsignal\n')
        for i in count:
            f.write(f'{i}\t{count[i][0]}\t{count[i][1]}\t{count[i][2]}\t{count[i][3]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python 29.make_reflection_and_decide_ps_loci.py <path/to/ps_site> <outfile>')
        exit()

    main()
