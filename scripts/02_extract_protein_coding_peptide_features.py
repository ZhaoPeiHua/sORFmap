from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re

def get_codon_to_aa_map():
    """Return codon to amino acid mapping dictionary"""
    return {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP",
        "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

def translate_sequence(coding_sequence, codon_map):
    """
    Translate DNA sequence to protein sequence.
    
    Parameters:
    -----------
    coding_sequence : str
        DNA coding sequence
    codon_map : dict
        Codon to amino acid mapping
        
    Returns:
    --------
    str : Translated protein sequence
    """
    peptide = ""
    for i in range(0, len(coding_sequence), 3):
        if i + 2 < len(coding_sequence):
            codon = coding_sequence[i:i+3]
            peptide += codon_map[codon]
    return peptide

def extract_protein_coding_peptides_features(input_filepath, output_filepath):
    """
    Extract protein features from coding sequences.
    
    Parameters:
    -----------
    input_filepath : str
        Path to input file with coding sequences
    output_filepath : str
        Path to output file for protein features
    """
    codon_map = get_codon_to_aa_map()
    
    with open(input_filepath, 'r') as input_file, open(output_filepath, 'w') as output_file:
        # Write header
       # output_file.write("transcript_id\taromaticity\tinstability_index\tisoelectric_point\thydropathy\tpeptide_sequence\n")
        
        for line in input_file:
            line_components = line.strip().split("\t")
            coding_sequence = line_components[7]
            
            if re.search(r'^ATG', coding_sequence):
                peptide_sequence = translate_sequence(coding_sequence.strip(), codon_map)
                
                if len(peptide_sequence) > 9 and not re.search("(STOP).*", peptide_sequence):
                    # Create transcript identifier
                    transcript_id = '_'.join([
                        line_components[0],  # ENST
                        line_components[1],  # ENSG
                        line_components[2],  # tstart
                        line_components[3],  # tend
                        line_components[5]   # strand
                    ])
                    
                    # Calculate protein features
                    protein_analysis = ProteinAnalysis(peptide_sequence)
                    protein_features = [
                        transcript_id,
                        f"{protein_analysis.aromaticity():.4f}",
                        f"{protein_analysis.instability_index():.4f}",
                        f"{protein_analysis.isoelectric_point():.4f}",
                        f"{protein_analysis.gravy():.4f}",
                        peptide_sequence
                    ]
                    
                    output_file.write('\t'.join(protein_features) + '\n')

# Main execution
if __name__ == "__main__":
    input_filepath = "../data/input/PC.txt"
    output_filepath = "../data/intermediate/protein_coding_peptides_features.txt"
    extract_protein_coding_peptides_features(input_filepath, output_filepath)
