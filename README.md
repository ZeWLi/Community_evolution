# Community_evolution
Scripts that used for the analysis of community evolutionary rates and natural selection

## PAML pipeline

`0.batch_run_paml_pipeline_for_tree_and_msa.py`

The purpose of this script is to use scripts/software to prepare all the gene clusters in a genome cluster, nucleotide alignment, and all the phylogenetic trees as inputs for the PAML.

```bash
python /data/lizw/script/paml_pipeline/batch_run_paml_pipeline.py 0.genomes
```

`batch_run_gene_prediction.pl`

The purpose of this script is to run Prodigal for gene prediction.

```bash
perl /data/script/batch_run_everything/batch_run_gene_prediction.pl 0.genomes
```

`seqs.trim.pl`

The purpose of this script is to filt proteins < 33 aa (approximately 100 bp).

```bash
perl /data/script/assemble.statistic/seqs.trim.pl all_pro.fa all_pro.trim33.fa 33
```

`3.batch_run_all_to_all_blast_pro_revised.py`

The purpose of this script is to run all-to-all blast.

```bash
python /data/lizw/script/paml_pipeline/3.batch_run_all_to_all_blast_pro_revised.py -t 50
```

`4.get_rBBH_from_table_m8_identity.pl`

The purpose of this script is to select blast hits with identity ≥ 90%, e-value < 1e-5, and shorter sequence with ≥ 50% length of the longer sequence.

```bash
perl /data/lizw/script/paml_pipeline/4.get_rBBH_from_table_m8_identity.pl blast.out blast_90.out
```

`5.parse_mcl_out_for_cluster_revised.py`

The purpose of this script is to parse mcl output with genes ≥ 5.

```bash
python3 /data/lizw/script/paml_pipeline/5.parse_mcl_out_for_cluster_revised.py -i out.blast_abc_mcl.out.I14 -o clstr -g all_gene.fa -p all_pro.fa
```

`6.batch_run_mafft.py`

The purpose of this script is to run mafft for all the gene clusters.

```bash
python /data/lizw/script/paml_pipeline/6.batch_run_mafft.py clstr 50
```

`7.batch_run_iqtree.py`

The purpose of this script is to run iqtree for all msa files.

```bash
python /data/lizw/script/paml_pipeline/7.batch_run_iqtree.py msa_file 50 tree
```

`8.batch_run_pal2nal.py`

The purpose of this script is to run pal2nal to prepare the nucleotide alignments.

```bash
python /data/lizw/script/paml_pipeline/8.batch_run_pal2nal.py nuc msa msa_nuc
```

`final_version_mpi.py`

The purpose of this script is to run PAML on super computing clusters.

```bash
python3 final_version_mpi.py msa_nuc tree paml
```

## PAML_results_parser

`batch_run_ps_check_and_site_extraction.py`

The purpose of this script is to use scripts to parse the PAML results, including positively selected (PS) signal, PS loci, and omega value, and summarize data.

```bash
nohup python /data/lizw/script/paml_pipeline/batch_run_ps_check_and_site_extraction.py parse_1_2 &
```

`11.batch_run_chi2_test.py`

The purpose of this script is to run chi-squared test and run p-value adjustment.

```bash
python /data/lizw/script/paml_pipeline/11.batch_run_chi2_test.py cluster419_paml p_adj.txt
```

`run_adjust_p.R`

The purpose of this script is to run p-value correction using Bonferroni method.

```bash
Rscript /data/lizw/script/paml_pipeline/run_adjust_p.R temp.txt merge.txt
```

`12.get_clstr_selected.py`

The purpose of this script is to find gene clusters detected with significantly positive selection signals with both M1a vs M2a and M7 vs M8 tests with corrected p-value < 0.05.

```bash
python /data/lizw/script/paml_pipeline/12.get_clstr_selected.py p_adj.txt clstr ps_clstr 1
```

`14.get_w_for_ps_gene.py`

The purpose of this script is to parse the omega output of M0.

```bash
python /data/lizw/script/paml_pipeline/14.get_w_for_ps_gene.py cluster419_and_selected/cluster419 ps_omega.txt
```

`15.merge_p_anno_and_omega.py`

The purpose of this script is to summarize data of omega and PS signal.

```bash
python /data/lizw/script/paml_pipeline/15.merge_p_anno_and_omega.py p_adj.txt ps_omega.txt ps_summary.txt
```

`16.parse_ps_site_from_rst.py`

The purpose of this script is to parse the PS loci from the output of “codeml”.

```bash
python /data/lizw/script/paml_pipeline/16.parse_ps_site_from_rst.py ps_clstr cluster419_paml msa ps_site ps_site_summary.txt
```

`17.count_sample_ps_site_with_vote.py`

The purpose of this script is to count the PS loci, and if given a gene cluster with sequences from identical samples, the script vote for the majority of amino acid usage.

```bash
python /data/lizw/script/paml_pipeline/17.count_sample_ps_site_with_vote.py ps_loci_extracellular ps_loci_extracellular.txt
```

`18.detect_change_in_ps_sites.py`

The purpose of this script is to detect the amino acid changes of the PS loci.

```bash
python /data/lizw/script/paml_pipeline/18.detect_change_in_ps_sites.py ps_site ps_change_count.txt
```

`batch_count_ps_clstr.py`

The purpose of this script is to summary the data mentioned above.

```bash
python /data/lizw/script/paml_pipeline/batch_count_ps_clstr.py parse_1_2 summary_parse_1_2.txt
```

## EPS_analysis

`27.batch_run_phobius_and_parse.py`

The purpose of this script is to find the longest sequences in the gene cluster, and use it to run Phobius.

```bash
nohup python /data/lizw/script/paml_pipeline/27.batch_run_phobius_and_parse.py all &
```

`convert_ps_site_count_to_proportion.py`

The purpose of this script is to convert the PS loci values to proportional data, also remove the gaps in PS loci summary.

```bash
python /data/lizw/script/paml_pipeline/convert_ps_site_count_to_proportion.py 20240228_ps_count_site.txt 0 20240228_ps_count_rm_gap.txt
```

`29.make_reflection_and_decide_ps_loci.py`

The purpose of this script is to summarize the total count of PS loci for all genome clusters.

```bash
nohup python /data/lizw/script/paml_pipeline/29.make_reflection_and_decide_ps_loci.py all 20240826_ps_location_rerun.txt &
```

`31.batch_count_aa_species_ps_site.py` 

The purpose of this script is to count and summarize the number of amino acid used of all the PS loci.

```bash
nohup python /data/lizw/script/paml_pipeline/31.batch_count_aa_species_ps_site.py all 20240304_ps_site_aa_species.txt &
```

`36.parse_phobius_for_gene_sublocation.py`

The purpose of this script is to parse the protein regions of subcellular location.

```bash
nohup python /data/lizw/script/paml_pipeline/36.parse_phobius_for_gene_sublocation.py all 20240611_gene_sublocation.txt &
```

## Other scripts

`13.parse_anno_result_and_vote_kegg.py`

The purpose of the script is to parse the KEGG annotations of PS genes (we have annotated all MAGs), and vote for the majority of the function for gene clusters.

```bash
nohup python /data/lizw/script/paml_pipeline/13.parse_anno_result_and_vote_kegg.py ps_clstr /data/lizw/0.all_sample_ecology/new_paml/all_kofam ps_clstr out.txt &
```

`check_identical_sequence.py`

The purpose of this script is to find the gene clusters with identical sequences that considered as “no variations”.

```bash
nohup python /data/lizw/script/paml_pipeline/check_identical_sequence.py middle_done_all 20240124_identical_cluster.txt &
```

## Computing_scripts
`unweighted_unifrac.R`
The purpose of this script is to calculate unweighted unifrac distance and perform NMDS and ANOSIM analysis.

`evolutionary_rates.R`
The purpose of this script is to calculate the evolutionary rates of samples based on their branch lengths in a phylogenetic tree.

`all_aa_gravy.R`
The purpose of the script is to calculate GRAVY value of community-level amino acid composition.
