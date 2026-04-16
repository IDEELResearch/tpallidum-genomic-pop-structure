import os

def get_genes_from_files(input_files):
    """
    Extracts unique gene names from a list of input BED files.
    
    Args:
        input_files (list): List of paths to BED files.
        
    Returns:
        list: A list of unique gene names.
    """
    genes = set()
    for file in input_files:
        with open(file, 'r') as f:
            for line in f:
                # Assuming the gene name is located in the fourth column of the BED file.
                tag = line.strip().split('\t')[3]
                gene = tag.split('_')[0]
                genes.add(gene)
    return list(genes)

def ensure_directory_exists(dir_path):
    """
    Ensures that a given directory exists, and if not, creates it.
    
    Args:
        dir_path (str): The path to the directory to check or create.
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

