import os
import shutil

def remove_analysis_directories(parent_directory):
    genome_directories = [name for name in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, name))]

    for genome_name in genome_directories:
        genome_directory = os.path.join(parent_directory, genome_name)
        analysis_directory = os.path.join(genome_directory, "analysis")

        if os.path.exists(analysis_directory) and os.path.isdir(analysis_directory):
            try:
                shutil.rmtree(analysis_directory)
                print(f"Removed 'analysis' directory in {genome_name}")
            except OSError as e:
                print(f"Error removing 'analysis' directory in {genome_name}: {e}")



# Example usage:
parent_directory = "/ResearchData/Microbiome/gblast"  # Replace this with the actual parent directory path

remove_analysis_directories(parent_directory)
