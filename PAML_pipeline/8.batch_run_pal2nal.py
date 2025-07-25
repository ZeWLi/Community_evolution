import os
import sys

def run_pal2nal(nuc,alignment,out):
    os.system(f'perl /data/software/metaomics/pal2nal.v14/pal2nal.pl {alignment} {nuc} -output fasta > {out}')
    return 0

def batch_get_path(alignment,out):
    path_align = []
    name = []
    for root,dirs,files in os.walk(alignment):
        for file in files:
            path_align.append(os.path.join(root,file))
            name.append(file.rstrip('.mafft.fa'))
    try:
        os.mkdir(out)
    except:
        os.system(f'rm -r {out}')
        os.mkdir(out)
    return path_align,name

def main():
    (path,name) = batch_get_path(sys.argv[2],sys.argv[3])
    for i,j in enumerate(path):
        run_pal2nal(f'{sys.argv[1]}/{name[i]}.nuc.fa',j,f'{sys.argv[3]}/{name[i]}.aligned.nuc.fa')
    return 0

if __name__ =='__main__':
    if len(sys.argv) != 4:
        print('wrong argv')
        exit()

    main()
