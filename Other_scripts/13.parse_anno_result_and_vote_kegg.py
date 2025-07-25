import os
import sys
import re
from collections import Counter
import subprocess


# former version, you need to do gene prediction
def parse_anno_out(file):
    ko = []
    with open(file, 'r') as f:
        out = f.readlines()
        for line in out:
            if not line.startswith('#'):
                ko.append(re.findall(r'\bK\d+\b', line)[0])

    count = Counter(ko)
    try:
        maxi = max(count.values())
        vote_out = [k for k, v in count.items() if v == maxi][0]
        indice = [i for i, x in enumerate(ko) if x == vote_out][0]
        split = re.split(r'\ +', out[indice + 2], maxsplit=5)
        output = f'{split[1]}; {split[-1]}'.rstrip('\n')
    except ValueError:
        output = 0

    return file.rstrip('.kofam.rm.repeats.out'), output


def kegg_anno(file):
    kegg = []
    with open(file, 'r') as f:
        out = f.readlines()
        for line in out:
            if not line.startswith('#'):
                kegg.append(re.findall(r'\bK\d+\b', line)[0])

    count = Counter(kegg)
    try:
        maxi = max(count.values())
        vote_out = [k for k, v in count.items() if v == maxi][0]
        indice = [i for i, x in enumerate(kegg) if x == vote_out][0]
        split = re.split(r'\ +', out[indice+2], maxsplit=5)
        output = f'{vote_out}; {split[-1]}'.rstrip('\n')
    except ValueError:
        output = 0

    return file.split('/')[-1].rstrip('.kofam.rm.repeats.out'), output


def parse_anno_from_kegg_all(file, dir):
    gene_name = []

    res = subprocess.run(['grep', '>', file], capture_output=True, text=True)
    result = res.stdout.split('\n')
    for line in result:
        if line != '':
            gene_name.append(line.split(' # ')[0].lstrip('>'))

    ko = []
    out = []
    for gene in gene_name:
        gene_sep = gene.split('_')
        res1 = subprocess.run(['grep', '-w', f'{gene}', f'{dir}/{gene_sep[0]}_{gene_sep[1]}_{gene_sep[2]}_{gene_sep[3]}.pro.kofam.rm.repeats.out'],capture_output=True, text=True)
        if res1.stdout != '':
            ko.append(re.findall(r'\bK\d+\b', res1.stdout)[0])
            out.append(res1.stdout)

    count = Counter(ko)
    try:
        maxi = max(count.values())
        vote_out = [k for k, v in count.items() if v == maxi][0]
        indice = [i for i, x in enumerate(ko) if x == vote_out][0]
        split = re.split(r'\ +', out[indice], maxsplit=5)
        output = f'{split[1]}; {split[-1]}'.rstrip('\n')
    except ValueError:
        output = 0

    return file.split('/')[-1].rstrip('.fa'), output


def main():
    vote_out = {}

    to_be_anno_kegg = []
    for root, dirs, files in os.walk(sys.argv[1]):
        for file in files:
            key, value = parse_anno_from_kegg_all(os.path.join(root, file), sys.argv[2])
            vote_out[key] = value
            if value == 0:
                to_be_anno_kegg.append(os.path.join(root, file))

    # kegg
    try:
        os.mkdir('temp3')
    except FileExistsError:
        os.system('rm -r temp3')
        os.mkdir('temp3')

    for path in to_be_anno_kegg:
        os.system(f'cp {path} temp3')

    if os.path.exists("kegg"):
        os.system('rm -r kegg')

    os.system('perl /data/script/batch_run_everything/batch_run_func_annotation_v4.pl -i temp3 -d kegg -o kegg -e 1e-5 -p 30')

    os.chdir('kegg/kegg')
    os.system('mv */*.out ..')
    os.chdir('..')
    os.system('rm -r kegg')
    os.chdir('..')

    to_be_anno_pfam = []
    for root, dirs, files in os.walk('kegg'):
        for file in files:
            key, value = kegg_anno(os.path.join(root, file))
            vote_out[key] = value
            if value == 0:
                to_be_anno_pfam.append(file.rstrip('.kofam.rm.repeats.out'))

    os.system('rm -r temp3 kegg')

    with open(sys.argv[3], 'w') as f1:
        # f1.write('cluster\tfunction\n')
        for i in vote_out:
            if vote_out[i] == 0:
                f1.write(f'{i}\tHyp protein\n')
            else:
                f1.write(f'{i}\t{vote_out[i]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python 13.parse_anno_result_and_vote.py anno out.txt')
        exit()

    main()
