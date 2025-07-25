'''
Date:20230830
nohup python3 /data/lizw/script/paml_pipeline/5.parse_mcl_out_for_cluster_revised.py -i out.blast_abc_mcl.out.I14 -o clstr_test -g all_gene.fa &
'''
import os
import argparse
from Bio import SeqIO

def parse_mcl_out(file):
    clst = {}
    with open(file,'r') as f:
        count = 0
        for line in f:
            count += 1
            clst[f'cluster{count}'] = line.strip('\n').split('\t')
    return clst

def filt_small_clstr(clst):
    clstr_done = {}
    for i in clst:
        if len(clst[i]) > 4:
            clstr_done[i] = clst[i]
    return clstr_done

def get_clstr(dict,genefile,dir,profile):
    for i in dict:
        with open('wxy.txt','w') as f:
            for j in dict[i]:
                f.write(f'{j}\n')
        os.system(f'seqkit grep -f wxy.txt {genefile} > {dir}/{i}.nuc.fa')
        os.system(f'seqkit grep -f wxy.txt {profile} > {dir}/{i}.fa')
        os.system('rm wxy.txt')
    return 0

def clstr_selection(file,cutoff):
    records = list(SeqIO.parse(file,'fasta'))
    out = []
    rm = []
    length = []
    for i in records:
        length.append(len(i.seq))
    max_len = float(max(length))
    for j,k in enumerate(length):
        if (k / max_len) >= cutoff:
            out.append(records[j])
        else:
            rm.append(records[j].id)
    SeqIO.write(out,file,'fasta')
    return rm

def rm_small_clstr(dir):
    for root,dirs,files in os.walk(dir):
        os.chdir(dir)
        for file in files:
            if file.endswith('.fa'):
                num = eval(os.popen(f'grep ">" {file} | wc -l').readlines()[0])
                if num < 5:
                    os.system(f'rm {file}')
                    os.system(f'rm {file.rstrip(".fa")}.rm.txt')

def main():
    try:
        os.mkdir(args.o)
    except:
        pass

    dic = parse_mcl_out(args.i)
    dic_done = filt_small_clstr(dic)
    get_clstr(dic_done,args.g,args.o,args.p)

    for root,dirs,files in os.walk(args.o):
        os.chdir(args.o)
        for file in files:
            rm = clstr_selection(file,args.c)
            if rm != []:
                with open(f"{file.strip('.fa')}.rm.txt", 'w') as f:
                    for i in rm:
                        f.write(f'{i}\n')
        os.chdir('..')

    rm_small_clstr(args.o)
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='mcl out')
    parser.add_argument('-g',help='genefile')
    parser.add_argument('-p', help='protein file')
    parser.add_argument('-o', help='output clsr directory')
    parser.add_argument('-c', help='cutoff',type=int,default=0.5)
    args = parser.parse_args()

    main()