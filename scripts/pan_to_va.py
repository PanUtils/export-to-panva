#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:
# Script is used to pre-process PanTools pangenomes for variant analysis in PanVA.

# Imports:
from datetime import datetime
import glob
import logging
import os
import sys
import ast
import pandas as pd
import shutil
from configparser import ConfigParser
from multiprocessing import Pool
from functools import reduce
from itertools import repeat
import tqdm
# in git imports
from genome_annot import *
from create_homologies_json import *
from alignment_stats import *
from homgroup_seq_info import *
from run_per_group import *
from pheno_specific_var import *
from include_pheno_info import *
from make_alignments import *
from collect_pantools_data import *
from process_newick_trees import *
from config_logger import *

# Use logger from config_logger.py
logger = logging.getLogger("export-to-panva")

# Functions:
def main():
    """
    Main function of the pre-processing scripts. It is used to pre-process a PanTools pangenome for
    visualisation in PanVA.
    usage: python3 pan_to_va.py config.ini
    :return:
    """
    # Load configurations & prepare log file
    start_run = datetime.now()
    # Load configfile
    config = ConfigParser()
    config.read(str(sys.argv[1]))
    # Create log file
    log_file = str(config.get('GENERAL', 'logfile'))
    # Set run log file format
    logging.basicConfig(filename=log_file, format='%(asctime)s [%(levelname)s] - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', filemode='a', level=logging.DEBUG)
    logging.info('Pre-processing started.')
    logging.info("Configuration file used: '{}'.".format(str(sys.argv[1])))

    # GENERAL settings:
    # Set PanTools database location from config
    pangenome_path = config.get('GENERAL', 'pangenome_path')
    pangenome_name = os.path.basename(os.path.normpath(pangenome_path))

    group_ext = config.get('GENERAL', 'group_extension')
    group_name = os.path.basename(os.path.normpath(group_ext))

    # Setting up file locations for output
    panva_path = config.get('GENERAL', 'panva_dir')
    panva_dir_name = pangenome_name + "_" + group_name + "_panva_input"
    panva_dir_path = os.path.join(panva_path, panva_dir_name)
    try:
        os.makedirs(panva_dir_path, exist_ok=True)
    except FileExistsError:
        pass
    panva_path = os.path.join(panva_dir_path, "homology")
    # panva_path = os.path.join(panva_path, panva_dir_name)
    logger.info(f"Output path: {panva_path}")
    try:
        os.mkdir(panva_path)
    except FileExistsError:
        pass
    pan_grp_path = os.path.join(pangenome_path, group_ext)

    # Max Number of cores allowed to be used at a time
    core_count = config.getint('GENERAL', 'core_count')
    # check if the given number of cores is possible on system else overwrite to the max possible
    if core_count > os.cpu_count():
        core_count = os.cpu_count()
        logging.warning("Number of cores exceeded available. Run continued using '{}' cores.".format(os.cpu_count()))
    else:
        logging.info("Number of cores used '{}'".format(core_count))

    # MSA SETTINGS:
    min_align = config.getint('MSA', 'min_align_len')
    logger.info("Minimum alignment length of {}".format(min_align))

    min_num_memb = config.getint('MSA', 'min_num_members')
    logger.info("Minimum number of members required in a homology group {}".format(min_num_memb))
    min_num_uniq = config.getint('MSA', 'min_uniq_members')

    logger.info("Minimum number of unique members required in a homology group {}".format(min_num_uniq))

    msa_type = config.get('MSA', 'msa_type')
    if config.getboolean('MSA', 'trimmed'):
        trimmed = "_trimmed"
    else:
        trimmed = ""

    if msa_type == 'msa_per_group_var':
        seqtyping = 'var'
    else:
        seqtyping = 'nuc'

    # HOMOLOGY settings:
    selection = config.get('HOMOLOGY', 'selection')

    # PHENOMETA settings:
    inc_pheno = config.getboolean('PHENOMETA', 'pheno_info')
    pheno_spvar = config.get('PHENOMETA', 'pheno_var')
    if pheno_spvar == '' or pheno_spvar == ' ':
        pheno_spvar = None
        logger.info("pheno_var is 'None'")

    acc_meta = config.get('PHENOMETA', 'reseq_meta')

    # GENE CLASSIFICATION settings:
    gene_class = str(config.get('GENE_CLASS', 'classification'))
    gene_class_path = os.path.join(pangenome_path, gene_class)
    gene_class_path = os.path.join(gene_class_path, "classified_groups.csv")

    # Functions
    grp_id_list = collect_data(pan_grp_path, selection)

    # start a multiprocess run get align_info of every groups
    with Pool(core_count) as pool_1:
        # get_alignment_info found in: alignment_stats.py
        alignment_info = pool_1.starmap(get_alignment_info, zip(grp_id_list, repeat(seqtyping), repeat(trimmed)))
    pool_1.close()

    # report settings used during the run
    logging.info("Filtering with minimal alignment length of '{}'.".format(min_align))
    logging.info("Filtering with minimal number of members per homology group of '{}'.".format(min_num_memb))
    logging.info("Filtering with minimal number of unique members per homology group '{}'.".format(min_num_uniq))

    # Apply the filters
    # alignment_info found in alignment_stats.py
    df_aligned, id_list_filt = filter_on_alignment(alignment_info, pan_grp_path, panva_path, min_align, min_num_memb,
                                                  min_num_uniq)

    df_aligned.to_csv(os.path.join(panva_path, 'tmp_filtered_groups.csv'), index=False)

    filter_time = datetime.now()
    logging.info("Time spend gathering and filtering the pangenome data {}.".format(filter_time - start_run))

    if inc_pheno == True:
        logger.info('Include phenotype and metadata')
        # makes df for all genomes linking meta and pheno data.
        df_phenos = pheno_meta_maker(pangenome_path)
    else:
        df_phenos = None
    # check if accession meta data is given/exists
    if os.path.isfile(acc_meta):
        logger.info('Include phenotype and metadata of accessions')
        acc_meta_df = pd.read_csv(acc_meta)
        # check if the column accession_id exists as is expected
    elif acc_meta == '' or acc_meta == ' ':
        acc_meta_df = None
    else:
        raise ValueError("the given path to accession meta/pheno data does not exist or is not empty in config.")

    print("*Note progress bar: Updates on next group start not on group done")
    with Pool(core_count) as pool_2:
        # prep_group_pheno found in run_per_group.py
        meta_info = pool_2.starmap(prep_group_pheno, tqdm.tqdm(zip(id_list_filt, repeat(panva_path), repeat(seqtyping),
                                                         repeat(pheno_spvar), repeat(df_phenos), repeat(acc_meta_df),
                                                                   repeat(trimmed)), total=len(id_list_filt)))
    pool_2.close()

    indiv_time = datetime.now()
    logging.info('Time passed making all individual homology group files: {}'.format(indiv_time - filter_time))
    logger.info("Creating homologies.json")
    id_class = gene_classification(gene_class_path, id_list_filt)
    if pheno_spvar is not None:
        create_homologies(pangenome_path, panva_path, df_aligned, id_class, pheno_var=meta_info)
    else:
        create_homologies(pangenome_path, panva_path, df_aligned, id_class, pheno_var=None)
    json_time = datetime.now()
    logging.info('Time passed making homologies.json: {}'.format(json_time - indiv_time))

    # Annotations
    logger.info("Making annotations.csv, last pre-processing step.")
    print("*Note progress bar: Updates on next group start not on group done")
    with Pool(core_count) as pool_3:
        pool_3.starmap(pool3wrap, tqdm.tqdm(zip(id_list_filt, repeat(panva_path), repeat(msa_type), repeat(trimmed)),
                                            total=len(id_list_filt)))
    pool_3.close()
    annot_time = datetime.now()
    logging.info('Time passed making annotation(s) files for the homology groups: {}'.format(annot_time-json_time))
    logger.info("Finished making the annotation(s) files")

    # Trees
    logger.info("Include trees if present")
    # Dictionary mapping tree types to their respective subdirectories and filenames
    tree_mappings = {
        'gene_distance': ('gene_classification', 'genome_gene_distance.tree', 'gene_distance.txt'),
        'kmer_distance': ('kmer_classification', 'genome_kmer_distance.tree', 'kmer_distance.txt'),
        'ani': ('ANI/MASH', 'ANI.newick', 'ANI.txt'),
        'core_snp': ('core_snp_tree', 'informative.fasta.treefile', 'core_snp.txt'),
    }
    for tree_name in config['TREES']:
        if config.getboolean('TREES', tree_name):
            if tree_name in tree_mappings:
                mapping = tree_mappings[tree_name]
                worker(pangenome_path, panva_path, mapping)

    final_time = datetime.now()
    logger.info("Process completed")
    logging.info("Total time passed: {} \n".format(final_time - start_run))
    return


if __name__ == "__main__":
    main()
