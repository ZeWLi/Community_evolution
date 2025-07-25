import os
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_ps_site(file):
    if file.endswith('.rst'):
        command = f"cat -n {file} | grep 'Positively selected sites' | awk '{{print $1}}'"
        res = subprocess.run(command, capture_output=True, text=True,shell= True)
        indice = res.stdout.split('\n')

        indice_ps = []
        indice_ps.append(eval(indice[1]))
        indice_ps.append(eval(indice[4]))

        with open(file, 'r') as f:
            m2 = f.readlines()[indice_ps[0] + 3:]
            f.seek(0)
            m8 = f.readlines()[indice_ps[1] + 3:]

        for j,k in enumerate(m2):
            if k == '\n':
                end_m2 = j
                break

        if m8 == []:
            end_m8 = 0
        else:
            for s,q in enumerate(m8):
                if q == '\n':
                    end_m8 = s
                    break

        m2_ps = m2[0:end_m2]
        m2_sig_ps = {}
        for i in m2_ps:
            split = i.split()
            if '*' in split[2]:
                m2_sig_ps[split[0]] = split[3]+ split[4] + split[5]

        m8_ps = m8[0:end_m8]
        m8_sig_ps = {}
        for i in m8_ps:
            split = i.split()
            if '*' in split[2]:
                m8_sig_ps[split[0]] = split[3] + split[4] + split[5]

        sig_ps = m2_sig_ps

        if m8_sig_ps != {}:
            for i in m8_sig_ps:
                if i in m2_sig_ps:
                    if eval(m2_sig_ps[i].split('+-')[0]) >= eval(m8_sig_ps[i].split('+-')[0]):
                        sig_ps[i] = m2_sig_ps[i]
                else:
                    sig_ps[i] = m8_sig_ps[i]

    return sig_ps


def output_ps_site(dic,file,outfile):
    with open(file,'r') as f:
        records = SeqIO.to_dict(SeqIO.parse(f,'fasta'))

    out = {}

    for seq_id,seq_record in records.items():
        seq = ''
        for i in dic:
            index = int(i)
            seq += str(seq_record.seq[index- 1])

        seq_record = Seq(seq)
        out[seq_id] = SeqRecord(seq_record,id = seq_id,description='')

    if seq != '':
        SeqIO.write(list(out.values()), outfile, 'fasta')


def main():
    # get ps cluster
    name = []

    ps_clstr = os.listdir(sys.argv[1])
    for i in ps_clstr:
        name.append(i.rstrip('.fa'))

    pwd = os.getcwd()
    path = os.path.abspath(sys.argv[3])
    out_path = os.path.abspath(sys.argv[4])

    os.chdir(sys.argv[2])

    try:
        os.mkdir(out_path)
    except FileExistsError:
        pass

    all = {}
    for i in name:
        if not 'mafft' in i:
            ps_site = parse_ps_site(f'{i}.rst')
            for j in ps_site:
                all[f'{i}_{j}'] = ps_site[j]
            output_ps_site(ps_site,f'{path}/{i}.mafft.fa',f'{out_path}/{i}_ps.fa')

    os.chdir(pwd)

    with open(sys.argv[5],'w') as f:
        f.write('cluster\tsite\tomega\n')
        for i in all:
            f.write(f'{i.split("_")[0]}\t{i.split("_")[1]}\t{all[i]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print('Usage: python 16.py ')
        exit()
    main()
