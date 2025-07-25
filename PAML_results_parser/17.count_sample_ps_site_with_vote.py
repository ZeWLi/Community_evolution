import os
import sys
from Bio.Seq import Seq
from Bio import SeqIO


def count_sample_ps_site(file):
    dic = {}
    with open(file,'r') as f:
        dict_seq = SeqIO.to_dict(SeqIO.parse(f,'fasta'))

    # with voting
    for i in dict_seq:
        sample = f'{i.split("_")[1]}_{i.split("_")[0]}'
        if sample not in dic:
            dic[sample] = [str(dict_seq[i].seq)]
        else:
            dic[sample].append(str(dict_seq[i].seq))

    for i in dic:
        if len(dic[i]) == 1:
            dic[i] = dic[i][0]
        else:
            # vote for every position
            vote = {}
            for j in range(len(dic[i][0])):
                vote[j] = {}
                for k in dic[i]:
                    if k[j] not in vote[j]:
                        vote[j][k[j]] = 1
                    else:
                        vote[j][k[j]] += 1
                vote[j] = sorted(vote[j].items(),key=lambda x:x[1],reverse=True)[0][0]

            # output
            seq = ''
            for j in vote:
                seq += vote[j]
            dic[i] = seq

    return dic

def main():
    ps_site_all = os.listdir(sys.argv[1])

    sample_all = {}
    for i in ps_site_all:
        if i.endswith('.fa') or i.endswith('.fasta'):
            ps_count = count_sample_ps_site(os.path.join(sys.argv[1],i))
            for j in ps_count:
                if j not in sample_all:
                    sample_all[j] = ''
                    sample_all[j] += ps_count[j]
                else:
                    sample_all[j] += ps_count[j]

    # count frequency
    frequency = {}
    amino_acids = "ACDEFGHIKLMNPQRSTVWY-N" # 20 amino acids + lost(-)

    for i in sample_all:
        freq = {}
        for aa in amino_acids:
            count = Seq(sample_all[i]).count(aa)
            freq[aa] = count
        frequency[i] = freq

    with open(sys.argv[2],'w') as f:
        if sys.argv[2].endswith('.csv'):
            f.write('sample,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,-\n')
            for i in frequency:
                f.write(f'{i},{frequency[i]["A"]},{frequency[i]["C"]},{frequency[i]["D"]},{frequency[i]["E"]},{frequency[i]["F"]},{frequency[i]["G"]},{frequency[i]["H"]},{frequency[i]["I"]},{frequency[i]["K"]},{frequency[i]["L"]},{frequency[i]["M"]},{frequency[i]["N"]},{frequency[i]["P"]},{frequency[i]["Q"]},{frequency[i]["R"]},{frequency[i]["S"]},{frequency[i]["T"]},{frequency[i]["V"]},{frequency[i]["W"]},{frequency[i]["Y"]},{frequency[i]["-"]}\n')
        else:
            f.write('sample\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\t-\n')
            for i in frequency:
                f.write(f'{i}\t{frequency[i]["A"]}\t{frequency[i]["C"]}\t{frequency[i]["D"]}\t{frequency[i]["E"]}\t{frequency[i]["F"]}\t{frequency[i]["G"]}\t{frequency[i]["H"]}\t{frequency[i]["I"]}\t{frequency[i]["K"]}\t{frequency[i]["L"]}\t{frequency[i]["M"]}\t{frequency[i]["N"]}\t{frequency[i]["P"]}\t{frequency[i]["Q"]}\t{frequency[i]["R"]}\t{frequency[i]["S"]}\t{frequency[i]["T"]}\t{frequency[i]["V"]}\t{frequency[i]["W"]}\t{frequency[i]["Y"]}\t{frequency[i]["-"]}\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage:')
        exit()
    main()
