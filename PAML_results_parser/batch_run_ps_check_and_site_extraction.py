import os
import sys


def ps_pipeline(cluster):
    pwd = os.getcwd()
    os.chdir(cluster)

    a = os.system(f'python /data/lizw/script/paml_pipeline/11.batch_run_chi2_test.py {cluster}_paml p_adj.txt')
    b = os.system('python /data/lizw/script/paml_pipeline/12.get_clstr_selected.py p_adj.txt clstr ps_clstr 1')
    c = os.system(f'python /data/lizw/script/paml_pipeline/14.get_w_for_ps_gene.py ps_clstr {cluster}_paml ps_omega.txt')
    d = os.system('python /data/lizw/script/paml_pipeline/15.merge_p_anno_and_omega.py p_adj.txt ps_omega.txt ps_summary.txt')
    e = os.system(f'python /data/lizw/script/paml_pipeline/16.parse_ps_site_from_rst.py ps_clstr {cluster}_paml msa ps_site ps_site_summary.txt')
    f = os.system('python /data/lizw/script/paml_pipeline/17.count_sample_ps_site_with_vote.py ps_loci_extracellular ps_loci_extracellular.txt')
    g = os.system(f'python /data/lizw/script/paml_pipeline/18.detect_change_in_ps_sites.py ps_site ps_change_count.txt')

    os.chdir(pwd)

    lis = [a, b, c, d, e, f, g]
    err = []
    if sum(lis) != 0:
        for i,j in enumerate(lis):
            if j >= 0 :
                err.append(i+1)

    if err != []:
        return err
    else:
        return 0


def main():
    lis = os.listdir(sys.argv[1])

    err_raised = {}

    pwd = os.getcwd()
    os.chdir(sys.argv[1])
    for clstr in lis:
        err = ps_pipeline(clstr)
        if err != 0:
            err_raised[clstr] = err

    os.chdir(pwd)

    os.system(f'python /data/lizw/script/paml_pipeline/batch_count_ps_clstr.py {sys.argv[1]} summary_{sys.argv[1]}.txt')

    if err_raised != {}:
        print(err_raised)
        with open(f'err_{sys.argv[1]}.txt','w') as f:
            for i in err_raised:
                f.write(f'{i}\t{err_raised[i]}\n')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage:')
        exit()

    main()
