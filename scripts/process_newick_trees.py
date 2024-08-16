#!/ubin/python

# Astrid van den Brandt
# Date start: 13-8-2024
# WUR Bioinformatics CPG

# Imports
import os
import re
import logging

# Use logger from config_logger.py
logger = logging.getLogger("export-to-panva")

# Functions
def clean_newick_labels(newick_string):
    def clean_label(label):
        # Assume labels are in the format "1_gene1" and we  want to extract "1"
        cleaned = label.split('_')[0]
        return cleaned

    # Use regex to find all labels and replace them with cleaned versions
    cleaned_newick = re.sub(r'(\d+_\S*?)(:)', lambda m: clean_label(m.group(1)) + m.group(2), newick_string)

    return cleaned_newick

def process_tree_file(tree_file_path, output_file_path):
    """
    Process a Newick tree file: clean the labels and save the cleaned tree.
    """
    if os.path.exists(tree_file_path):
        try:
            with open(tree_file_path, 'r') as file:
                newick_string = file.read()

            cleaned_newick = clean_newick_labels(newick_string)

            with open(output_file_path, 'w') as file:
                file.write(cleaned_newick)

            logger.info(f"Processed and cleaned tree file: {output_file_path}")

        except Exception as e:
            logger.error(f"An error occurred while processing the file {tree_file_path}: {e}")
    else:
        logger.warning(f"File not found: {tree_file_path}")

def worker(pangenome_path, panva_path, tree_mapping):
    try:
        tree_folder, tree_file, output_file = tree_mapping
        tree_folder_path = os.path.join(pangenome_path, tree_folder)
        tree_file_path = os.path.join(tree_folder_path, tree_file)
        out_path = os.path.join(panva_path, output_file)

        process_tree_file(tree_file_path, out_path)

    except Exception as e:
        logger.error(f"Error in worker: {e}")