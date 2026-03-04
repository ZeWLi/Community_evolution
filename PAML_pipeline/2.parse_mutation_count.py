import os, sys
from ete3 import Tree
from Bio import AlignIO
import pandas as pd


def parse_tree_total_length(tree):
    treefile = Tree(tree)
    return sum(node.dist for node in treefile.traverse())


def parse_alignment_total_length(alignment):
    align = AlignIO.read(alignment, "fasta")
    return align.get_alignment_length()


def main():
    pwd = os.getcwd()
    base_dir = sys.argv[1]

    cluster_dirs = [d for d in os.listdir(base_dir)
                    if os.path.isdir(os.path.join(base_dir, d))]

    print(f"Found {len(cluster_dirs)} cluster directories")
    results = []

    for cluster_name in cluster_dirs:
        cluster_path = os.path.abspath(os.path.join(base_dir, cluster_name))
        align_path = os.path.join(cluster_path, 'msa_nuc_muscle')  # pal2nal file
        muscle_path = os.path.join(cluster_path, 'msa_muscle')  # muscle file
        tree_path = os.path.join(cluster_path, 'fasTtree')

        print(f"Processing cluster: {cluster_name}")

        if not os.path.exists(align_path):
            print(f"  Warning: {align_path} not found")
            continue

        if not os.path.exists(tree_path):
            print(f"  Warning: {tree_path} not found")
            continue

        align_files = [f for f in os.listdir(align_path) if f.endswith('.muscle5.pal2nal.fa')]
        print(f"  Found {len(align_files)} alignment files")

        for file in align_files:
            base_name = file.replace(".muscle5.pal2nal.fa", "")

            pal2nal_file = os.path.join(align_path, file) 
            muscle_file = os.path.join(muscle_path, f'{base_name}.muscle5.fa')
            tree_file = os.path.join(tree_path, f'{base_name}.muscle5.pal2nal.fasttree.newick')

            print(f"    Processing: {base_name}")
            print(f"      Alignment: {pal2nal_file}")
            print(f"      Tree: {tree_file}")

            if not os.path.exists(pal2nal_file):
                print(f"      Error: Alignment file not found")
                continue

            if not os.path.exists(tree_file):
                print(f"      Error: Tree file not found")
                continue

            try:
                tree_length = parse_tree_total_length(tree_file)
                align_length = parse_alignment_total_length(pal2nal_file)

                results.append({
                    'cluster': cluster_name, 
                    'gene_cluster': base_name,
                    'tree_length': tree_length,
                    'alignment_length': align_length,
                    'mutation': tree_length * align_length
                })

                print(f"      Success: tree_len={tree_length:.4f}, align_len={align_length}")

            except Exception as e:
                print(f"      Error processing {file}: {e}")

    print(f"\nTotal results: {len(results)}")

    if results:
        # output to csv files
        results_df = pd.DataFrame(results)
        output_file = os.path.join(pwd, sys.argv[2])
        results_df.to_csv(output_file, index=False)
        print(f"Results saved to: {output_file}")
        print(results_df.head())
    else:
        print("No results to save!")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python parse_mutation_count.py <base_directory> <output_csv>')
        sys.exit(1)

    main()