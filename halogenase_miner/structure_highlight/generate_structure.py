import sys
import subprocess

def generate_predicted_structure(fasta_file, output_dir):
    cmd = [
        "colabfold_batch",
        "--num-models", "1",
        "--num-recycle", "5",
        "--max-msa", "16:32",
        "--num-relax", "0",
        "--save-all",
        fasta_file,
        output_dir]
    esm_structure_prediction = subprocess.run(cmd)

    return f"The PDB file of the predicted structure is deposited in {output_dir}"
