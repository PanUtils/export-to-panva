#!/usr/bin/python

# Author: Sander Vlugter
# Start date v1: 21-12-2022
# WUR - Bioinformatics
# Functionality:

# Imports:
from homgroup_seq_info import *
import pheno_specific_var
from make_alignments import *
from include_pheno_info import *


# ran in pool
def prep_group_pheno(hom_id, panva_p, seqtype='nuc', pheno_var=None, meta=None, var_meta=None, trimmed="_trimmed"):
    """
    Flashes out the alignments for a homology group, by adding additional information if present.

    :param var_meta: reseq phenotypes
    :param hom_id: str - path to single homology grp id in alignments/../.. .
    :param panva_p: str - full path to output dir (serves as input for PanVA).
    :param seqtype: str - define the type of sequence default 'nuc_trimmed'.
    :param pheno_var: bool - indicates if there is data on meta/phenotype specific variant positions.
    :param meta: optional - can be None or dataframe containing meta/phenotype data on the aligned sequences.
    :return: meta_info - Dataframe - basis of the homologies.json
    """

    df_seq_info = create_seq_info(hom_id, panva_p)

    # makes sequences.csv
    df_all_seq = merge_sequences(hom_id, panva_p, method=seqtype, trimmed=trimmed)

    all_info_seq = pd.merge(df_all_seq, df_seq_info, on=['mRNA_id'])

    # if both not none ADD phenometa data AND pheno specific var
    if meta is not None:
        # makes metadata.csv
        group_meta = hom_group_pheno(hom_id, panva_p, meta, df_seq_info, var_meta)

        meta_info = alignment_posinfo(hom_id, panva_p, all_info_seq, phe_var=pheno_var, meta=group_meta,
                                      seqtype=seqtype, trimmed=trimmed)

    # if meta is None and pheno var is None just do positions
    else:
        # print("prep_group: meta NOT present and pheno_var NOT present")
        meta_info = alignment_posinfo(hom_id, panva_p, all_info_seq, phe_var=pheno_var, meta=meta, seqtype=seqtype,
                                      trimmed=trimmed)

    return meta_info
