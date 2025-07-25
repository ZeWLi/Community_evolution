import sys
import os
import numpy as np
from Bio.Phylo.PAML import codeml
from mpi4py import MPI


def parse_paml_file(dir1):
    alignment = []
    cluster = []
    for root, dirs, files in os.walk(dir1):
        for file in files:
            if file.endswith('.fa'):
                alignment.append(os.path.abspath(os.path.join(root, file)))
                cluster.append(file.rstrip('.aligned.nuc.fa'))

    return alignment, cluster


def codeml_with_para(alignment, tree, outfile):
    cml = codeml.Codeml()
    cml.alignment = alignment
    cml.tree = tree
    cml.out_file = outfile
    cml.set_options(noisy=9, verbose=1, runmode=0, seqtype=1, CodonFreq=2, clock=0, aaDist=0, model=0,
                    NSsites=[0, 1, 2, 3, 7, 8], icode=0, Mgene=0)
    cml.set_options(fix_kappa=0, kappa=2, fix_omega=0, omega=1, getSE=0, RateAncestor=0, Small_Diff=.5e-6, cleandata=0,
                    fix_blength=0)
    try:
        cml.run(verbose=True)
    except:
        return -1
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('wrong args!')
        sys.exit()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    try:
        os.mkdir(sys.argv[3])
    except:
        if rank == 0:
            pass

    cwd_tree = os.path.abspath(sys.argv[2])
    out_cwd = os.path.abspath(sys.argv[3])


    def process_paml_task(index):
        path = os.path.abspath(paml_parse[0][index])
        name = paml_parse[1][index]

        os.mkdir(f'paml_{name}')

        if not os.path.exists(f'paml_{name}'):
            os.mkdir(f'paml_{name}')

        try:
            os.chdir(f'paml_{name}')
            a = codeml_with_para(path, f'{cwd_tree}/{name}.mafft.fa.treefile', f'{name}.mlc')
            os.chdir('..')
            if a == -1:
                os.system(f'/bin/rm -r paml_{name}')
            os.chdir(f'paml_{name}')
            os.system(f'/bin/mv {name}.mlc {out_cwd}')
            os.system(f'/bin/mv rst {name}.rst')
            os.system(f'/bin/mv {name}.rst {out_cwd}')
            os.chdir('..')
            os.system(f'/bin/rm -r paml_{name}')
        except FileNotFoundError:
            pass

        return 0


    paml_parse = parse_paml_file(sys.argv[1])

    # Asynchronous non-blocking communication
    request_send = []
    request_recv = []

    if rank == 0:
        for dest_rank in range(1, size):
            # Send a signal to other ranks to start receiving
            data_to_send = np.array([1], dtype='i')
            request_send.append(comm.Isend(data_to_send, dest=dest_rank, tag=1))

    data_to_recv = np.empty(1, dtype='i')
    if rank != 0:
        # Receive the signal from rank 0 to start processing
        request_recv.append(comm.Irecv(data_to_recv, source=0, tag=1))

    # Wait for all requests to complete
    MPI.Request.waitall(request_send)
    MPI.Request.waitall(request_recv)

    # Process PAML tasks asynchronously
    for i in range(rank, len(paml_parse[0]), size):
        process_paml_task(i)

    MPI.Finalize()