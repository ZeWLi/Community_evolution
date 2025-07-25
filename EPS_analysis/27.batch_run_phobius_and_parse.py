#!~/anaconda3/bin/python

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_longest_seq_of_cluster_and_check_length(file):
    with open(file,'r') as f:
        sequences = [(seq_record.id, len(seq_record.seq),seq_record) for seq_record in SeqIO.parse(f, "fasta")]

    sorted_sequences = sorted(sequences, key=lambda x: x[1], reverse=True)

    if sorted_sequences[0][1] - sorted_sequences[1][1] > 33:
        selected = sorted_sequences[1][2]
        flag = 1
    else:
        selected = sorted_sequences[0][2]
        flag = 0

    selected.id = file.split('/')[-1].split('.')[0] + ',' + selected.id
    selected.description = ''
    #selected_seqrecord = selected, id = file.split('/')[-1].split('.')[0] + ',' + selected.id, description='')
    '''
    selected.id = file.split('/')[-1].split('.')[0] + ',' + selected.id
    selected.description = ''


    longest = ''
    for seq_id,seq_record in records.items():
        if len(seq_record.seq) > len(longest):
            longest = seq_record

    longest.id = file.split('/')[-1].split('.')[0] + ',' + longest.id
    longest.description = ''

    return longest
    '''
    return selected,flag


def parse_phobius_out(file):
    out = {}
    with open(file,'r') as f:
        for line in f:
            if line.startswith('ID'):
                id = line.split()[1]
                outside = ''
                inside = ''
                trans = ''
                signal = ''
            elif line.startswith('FT'):
                if line.split()[1] == 'TRANSMEM':
                    trans += line.split()[2] + '-' + line.split()[3] + ','
                elif line.split()[1] == 'DOMAIN':
                    if line.split()[4] ==  'CYTOPLASMIC.':
                        inside += line.split()[2] + '-' + line.split()[3] + ','
                    elif line.split()[4] == 'NON':
                        outside += line.split()[2] + '-' + line.split()[3] + ','
                # elif line.startswith('SIGNAL'):
                elif line.split()[1] == 'SIGNAL':
                    signal += line.split()[2] + '-' + line.split()[3] + ','
            elif line.startswith('//'):
                out[id] = [inside,outside,trans,signal]

    return out


def main():
    pwd = os.getcwd()
    lis = os.listdir(sys.argv[1])
    os.chdir(sys.argv[1])
    pwd1 = os.getcwd()

    for i in lis:
        os.chdir(i)
        sequence = []
        longer_sequence = []
        ps_clstr = os.listdir('ps_clstr')
        for j in ps_clstr:
            seq,flag = get_longest_seq_of_cluster_and_check_length(f'ps_clstr/{j}')
            sequence.append(seq)
            if flag == 1:
                longer_sequence.append(j)


        with open('phobius.fa','w') as f:
            SeqIO.write(sequence,f,'fasta')

        if len(longer_sequence) != 0:
            with open('unusual_long_cluster.txt','w') as f:
                for k in longer_sequence:
                    f.write(k + '\n')

        # phobius
        os.system('sed -i "s/*//g" phobius.fa')
        os.system('/data/software/metaomics/phobius/1.01/phobius.pl -long phobius.fa > phobius_out.txt')

        # parse
        out = parse_phobius_out('phobius_out.txt')
        with open('phobius.parse.out','w') as f:
            f.write('ID\tinside\toutside\ttransmembrane\tsignal\n')
            for k in out:
                f.write(f'{k}\t{out[k][0]}\t{out[k][1]}\t{out[k][2]}\t{out[k][3]}\n')

        try:
            os.mkdir('phobius')
        except FileExistsError:
            os.system('rm -r phobius')
            os.mkdir('phobius')

        os.system('mv phobius.fa phobius_out.txt unusual_long_cluster.txt phobius')

        os.chdir(pwd1)

    os.chdir(pwd)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('wrong args!')
        exit()

    main()
