import subprocess
import os
import sys


def get_w_from_mlc(file):
    res = subprocess.run(['grep','omega',file], capture_output=True, text=True) # M0 output
    result = res.stdout.split('\n')
    w = result[-2].split()[-1]

    return file.rstrip('.mlc'),w


def main():
    ps_clstr = []
    for root,dirs,files in os.walk(sys.argv[1]):
        for file in files:
            ps_clstr.append(file.rstrip('.fa'))

    pwd = os.path.abspath('.')
    os.chdir(sys.argv[2])
    omega = {}
    for clstr in ps_clstr:
        if not 'mafft' in clstr:
            out = get_w_from_mlc(f'{clstr}.mlc')
            omega[out[0]] = out[1]
        else:
            out = get_w_from_mlc(f'{clstr.split(".")[0]}.mlc')
            omega[out[0]] = out[1]
    os.chdir(pwd)

    with open(sys.argv[3],'w') as f:
        f.write('cluster\tomega\n')
        for i in omega:
            f.write(f'{i}\t{omega[i]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python 14.py dir_ps_clstr dir_mlc out.txt")
        exit()

    main()
