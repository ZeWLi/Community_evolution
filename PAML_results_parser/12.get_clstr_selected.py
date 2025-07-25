import os
import sys


def parse_p_adj_file(file,flag):
    selected_clstr = []
    with open(file,'r') as f:
        all_line = f.readlines()
        if flag == 0:
            for line in all_line[1:]:
                split = line.split('\t')
                if eval(split[-1]) < 0.05 or eval(split[-2]) < 0.05:
                    selected_clstr.append(split[0])

        elif flag != 0:
            for line in all_line[1:]:
                split = line.split('\t')
                if eval(split[-1]) < 0.05 and eval(split[-2]) < 0.05:
                    selected_clstr.append(split[0])

    return selected_clstr


def get_clstr_file(list,dir,outdir):
    try:
        os.mkdir(outdir)
    except:
        os.system(f'rm -r {outdir}')
        os.mkdir(outdir)

    for i in list:
        os.system(f'cp {dir}/{i}.fa {outdir}')


def main():
    if len(sys.argv) == 5:
        clstr = parse_p_adj_file(sys.argv[1],eval(sys.argv[4]))
        get_clstr_file(clstr,sys.argv[2],sys.argv[3])

    elif len(sys.argv) == 4:
        clstr = parse_p_adj_file(sys.argv[1], 0)
        get_clstr_file(clstr, sys.argv[2], sys.argv[3])

    else:
        print('Usage: python 12.get_clstr_selected.py p.txt clstr selected_clstr 0')


if __name__ == '__main__':
    main()

