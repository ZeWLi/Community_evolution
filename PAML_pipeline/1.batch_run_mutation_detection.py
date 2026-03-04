import os, sys
from multiprocessing import Pool, cpu_count
import time


def process_single_file(args):
    file, clstr_path, nuc_dir = args

    try:
        # 构建绝对路径
        input_file = os.path.join(clstr_path, file)
        base_name = file.strip(".fa")

        muscle_output = os.path.join(clstr_path, f'{base_name}.muscle5.fa')
        tree_output = os.path.join(clstr_path, f'{base_name}.muscle5.pal2nal.fasttree.newick')
        pal2nal_output = os.path.join(clstr_path, f'{base_name}.muscle5.pal2nal.fa')
        nuc_file = os.path.join(nuc_dir, f'{base_name}.nuc.fa')

        # muscle5 - 
        os.system(f'muscle5 -super5 {input_file} -output {muscle_output} -amino -threads 1')

        # pal2nal - 
        if os.path.exists(nuc_file):
            os.system(
                f'perl /data/software/metaomics/pal2nal.v14/pal2nal.pl {muscle_output} {nuc_file} -output fasta > {pal2nal_output}')

        # fasttree -
        os.system(f'FastTree -gtr {pal2nal_output} > {tree_output}')

        return True

    except Exception as e:
        print(f"Error processing {file}: {e}")
        return False


def collect_all_tasks(base_dir):

    tasks = []
    cluster_dirs = []

    all_clusters = [d for d in os.listdir(base_dir)
                    if os.path.isdir(os.path.join(base_dir, d))]

    for cluster_name in all_clusters:
        cluster_path = os.path.abspath(os.path.join(base_dir, cluster_name))
        clstr_path = os.path.join(cluster_path, 'clstr')
        nuc_path = os.path.join(cluster_path, 'nuc')

        if not os.path.exists(clstr_path) or not os.path.exists(nuc_path):
            continue

        align_files = [f for f in os.listdir(clstr_path) if f.endswith('.fa')]

        if align_files:
            cluster_dirs.append(cluster_path)
            nuc_dir_abs = os.path.abspath(nuc_path)
            clstr_path_abs = os.path.abspath(clstr_path)

            for file in align_files:
                tasks.append((file, clstr_path_abs, nuc_dir_abs))

    return tasks, cluster_dirs


def organize_all_files(cluster_dirs):
    print("Creating directories and organizing files...")

    for cluster_dir in cluster_dirs:
        cluster_name = os.path.basename(cluster_dir)
        print(f"Organizing {cluster_name}...")

        msa_muscle_dir = os.path.join(cluster_dir, 'msa_muscle')
        msa_nuc_muscle_dir = os.path.join(cluster_dir, 'msa_nuc_muscle')
        fasttree_dir = os.path.join(cluster_dir, 'fasTtree')

        os.makedirs(msa_muscle_dir, exist_ok=True)
        os.makedirs(msa_nuc_muscle_dir, exist_ok=True)
        os.makedirs(fasttree_dir, exist_ok=True)


        clstr_path = os.path.join(cluster_dir, 'clstr')

        muscle_files = [f for f in os.listdir(clstr_path) if f.endswith('.muscle5.fa')]
        for f in muscle_files:
            src = os.path.join(clstr_path, f)
            dst = os.path.join(msa_muscle_dir, f)
            if os.path.exists(src):
                os.rename(src, dst)

        tree_files = [f for f in os.listdir(clstr_path) if f.endswith('.muscle5.pal2nal.fasttree.newick')]
        for f in tree_files:
            src = os.path.join(clstr_path, f)
            dst = os.path.join(fasttree_dir, f)
            if os.path.exists(src):
                os.rename(src, dst)

        pal2nal_files = [f for f in os.listdir(clstr_path) if f.endswith('.muscle5.pal2nal.fa')]
        for f in pal2nal_files:
            src = os.path.join(clstr_path, f)
            dst = os.path.join(msa_nuc_muscle_dir, f)
            if os.path.exists(src):
                os.rename(src, dst)


def process_all_parallel(base_dir, num_processes=None):
    if num_processes is None:
        num_processes = cpu_count()

    base_dir = os.path.abspath(base_dir)

    print(f"Collecting tasks from {base_dir}...")

    tasks, cluster_dirs = collect_all_tasks(base_dir)

    if not tasks:
        print("No tasks found!")
        return

    print(f"Found {len(tasks)} files in {len(cluster_dirs)} clusters")
    print(f"Processing with {num_processes} processes...")

    start_time = time.time()

    with Pool(processes=num_processes) as pool:
        results = pool.map(process_single_file, tasks,chunksize = 1)

    end_time = time.time()

    successful = sum(results)
    failed = len(results) - successful

    print(f"Processing completed in {end_time - start_time:.2f} seconds")
    print(f"Successful: {successful}, Failed: {failed}")

    # 整理文件
    organize_all_files(cluster_dirs)
    print("All tasks completed!")


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        base_directory = sys.argv[1]
        num_processes = int(sys.argv[2]) if len(sys.argv) > 2 else 1

        if not os.path.exists(base_directory):
            print(f"Error: Directory {base_directory} does not exist")
            sys.exit(1)

        process_all_parallel(base_directory, num_processes)
    else:
        print("Usage: python script.py <base_directory> [num_processes]")