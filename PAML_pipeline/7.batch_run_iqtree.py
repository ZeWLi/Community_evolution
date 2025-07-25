'''
nohup python /data/lizw/script/paml_pipeline/7.batch_run_iqtree.py file_tree 50 tree_out &
'''

import os
import sys

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('specify tree directory')
        exit()

    for root,dirs,files in os.walk(sys.argv[1]):
        for file in files:
            if file.endswith('.fa') or file.endswith('.phy'):
                exit = os.system(f'iqtree -s {os.path.join(root,file)} -alrt 1000 -bb 1000 -nt {sys.argv[2]}')
                if exit == 512:
                    os.system(f'iqtree -s {os.path.join(root,file)} -alrt 1000 -nt {sys.argv[2]}')


    try:
        os.mkdir(sys.argv[3])
    except:
        os.system(f'rm -r {sys.argv[3]}')
        os.mkdir(sys.argv[3])

    os.system(f'cp {sys.argv[1]}/*.treefile {sys.argv[3]}')


