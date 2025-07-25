#!/usr/bin/python
'''
from a dRep cluster of genomes to paml output.
python /data/lizw/script/paml_pipeline/batch_run_paml_pipeline.py 0.genomes
'''
import sys
import os

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('wrong argv!')
        exit()

    # step 1: gene prediction and combination
    os.system(f'perl /data/script/batch_run_everything/batch_run_gene_prediction.pl {sys.argv[1]}')
    os.chdir(f'{sys.argv[1]}.pro')
    os.system('cat *.pro.fa > ../all_pro.fa')
    os.chdir('..')
    os.system('perl /data/script/assemble.statistic/seqs.trim.pl all_pro.fa all_pro.trim33.fa 33')
    os.chdir(f'{sys.argv[1]}.gene')
    os.system('cat *.gene.fa > ../all_gene.fa')
    os.chdir('..')

    # step 2: all-to-all blastp
    os.system(f'python /data/lizw/script/paml_pipeline/3.batch_run_all_to_all_blast_pro_revised.py -t 50')

    # step 3: mcl
    os.system('perl /data/lizw/script/paml_pipeline/4.get_rBBH_from_table_m8_identity.pl blast.out blast_90.out')
    os.system('cut -f 1,2,11 blast_90.out > blast_abc.out')
    os.system("mcxload -abc blast_abc.out --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o blast_abc_mcl.out -write-tab  blast_abc_tbl.out")
    os.system('mcl blast_abc_mcl.out -I 1.4 -use-tab blast_abc_tbl.out')

    # step 4: get cluster and select cluster
    os.system('python3 /data/lizw/script/paml_pipeline/5.parse_mcl_out_for_cluster_revised.py -i out.blast_abc_mcl.out.I14 -o clstr -g all_gene.fa -p all_pro.fa')
    os.mkdir('nuc')
    os.mkdir('rm')
    os.system('mv clstr/*.nuc.fa nuc')
    os.system('mv clstr/*.rm.txt rm')

    # step 5: tree building
    os.system('python /data/lizw/script/paml_pipeline/6.batch_run_mafft.py clstr 50')
    os.mkdir('msa_file')
    os.system('mv clstr/*.mafft.fa msa_file')
    os.system('cp -r msa_file msa')

    os.system('python /data/lizw/script/paml_pipeline/7.batch_run_iqtree.py msa_file 50 tree')

    # step 6: pal2nal
    os.system('python /data/lizw/script/paml_pipeline/8.batch_run_pal2nal.py nuc msa msa_nuc')
    os.system('rm -r msa_file')

    # step 7: paml
    #os.system('python /data/lizw/script/paml_pipeline/batch_run_paml_with_Bio_pipeline_revised.py msa_nuc tree paml 50')

    # step 8: check
    #os.system('python /data/lizw/script/paml_pipeline/check_paml_failed.py msa paml')