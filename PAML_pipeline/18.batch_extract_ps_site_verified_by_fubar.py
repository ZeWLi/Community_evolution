import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def output_ps_site(lis,file,outfile):
    with open(file,'r') as f:
        records = SeqIO.to_dict(SeqIO.parse(f,'fasta'))

    out = {}

    for seq_id,seq_record in records.items():
        seq = ''
        for i in lis:
            index = int(i)
            seq += str(seq_record.seq[index- 1])

        seq_record = Seq(seq)
        out[seq_id] = SeqRecord(seq_record,id = seq_id,description='')

    if seq != '':
        SeqIO.write(list(out.values()), outfile, 'fasta')


## parse ps site from output file
codeml_ps = {}
count_codeml = 0
with open(sys.argv[1]) as f:
    for line in f:
        split = line.strip().split('\t')
        if len(split) == 2:
            codeml_ps[split[0]] = split[1].split(',')
            count_codeml += len(codeml_ps[split[0]])


## parse ps site from hyphy fubar file
fubar_ps = {}
count_fubar = 0
with open(sys.argv[2]) as f:
    for line in f:
        split = line.strip().split('\t')
        if not split[0] in fubar_ps:
            fubar_ps[split[0]] = []

        fubar_ps[split[0]].append(split[1])
        count_fubar += 1


## verify whether it is also ps site in hyphy fubar result
verified_ps_loci = {}
count_verified = 0
for clstr in codeml_ps:
    if clstr in fubar_ps:
        for site in codeml_ps[clstr]:
            if site in fubar_ps[clstr]:
                count_verified += 1
                if not clstr in verified_ps_loci:
                    verified_ps_loci[clstr] = []

                verified_ps_loci[clstr].append(site)

print(f'Total PS sites in codeml: {count_codeml}')
print(f'Total PS sites in verified by hyphy FUBAR: {count_fubar}')
print(f'Total PS sites in verified by both codeml and hyphy FUBAR: {count_verified}')



pwd = os.getcwd()
out_path = os.path.abspath(sys.argv[4])

try:
    os.mkdir(out_path)
except FileExistsError:
    pass

for clstr in verified_ps_loci:
    temp_lis = verified_ps_loci[clstr]

    if temp_lis:
        output_ps_site(temp_lis,os.path.abspath(sys.argv[3]) + f'/{clstr}.mafft.fa',out_path + f'/{clstr}.ps.site.fa')


os.chdir(pwd)

with open(sys.argv[5],'w') as f:
    for clstr in verified_ps_loci:
        if verified_ps_loci[clstr]:
            f.write(f'{clstr}\t{",".join(verified_ps_loci[clstr])}\n')
