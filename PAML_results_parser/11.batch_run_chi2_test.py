import subprocess
import sys
import os
import re


def get_mlc_info(dir):
    path = []
    name = []

    for root,dirs,files in os.walk(dir):
        for file in files:
            if file.endswith('.mlc'):
                path.append(os.path.join(root,file))
                name.append(file.rstrip('.mlc'))

    return path,name


def get_lnL_and_run_chi2_test(mlcfile):
    res = subprocess.run(['grep','lnL',mlcfile],capture_output=True, text=True)
    result = res.stdout.split('\n')

    lnL = []
    for line in result:
        if line != '':
            try:
                lnL.append(eval(line.split()[-2]))
            except NameError:
                return [2,2]

    # run chi2
    res1 = subprocess.run(['chi2','2',f'{abs(2*(lnL[1]-lnL[2]))}'], capture_output=True, text=True)
    res2 = subprocess.run(['chi2','2',f'{abs(2*(lnL[4]-lnL[5]))}'], capture_output=True, text=True)

    patt = r'prob\s*=\s*([\d.e+-]+)'
    p = []
    p.append(re.findall(patt,res1.stdout)[0])
    p.append(re.findall(patt,res2.stdout)[0])

    return p


def main():
    p_dic = {}
    path,name = get_mlc_info(sys.argv[1])
    for i,j in enumerate(path):
        p_dic[name[i]] = get_lnL_and_run_chi2_test(j)

    with open('temp.txt','w') as f:
        f.write('cluster\tp_12\tp_78\n')
        for i in p_dic:
            f.write(f'{i}\t{p_dic[i][0]}\t{p_dic[i][1]}\n')

    os.system(f"Rscript /data/lizw/script/paml_pipeline/run_adjust_p.R temp.txt {sys.argv[2]}")

    os.system('rm temp.txt')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python 11.batch_run_chi2_test.py dir outfile with_R_file')
        exit()

    main()

