import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def extract_lncrna_peptides_features(input_filepath, output_filepath):
    """
    Extract protein features from translated lncRNA ORFs.
    
    Parameters:
    -----------
    input_filepath : str
        Path to input file with lncRNA ORF data
    output_filepath : str
        Path to output file for lncRNA protein features
    """
    with open(input_filepath, 'r') as input_file, open(output_filepath, 'w') as output_file:
        # Write header for output file
      #  output_file.write("lncrna_id\taromaticity\tinstability_index\tisoelectric_point\tgravy\tpeptide\n")
        
        # Skip header line if present
        next(input_file, None)
        
        for line in input_file:
            line_components = line.strip().split("\t")
            peptide_sequence = line_components[4]
            
            # Skip peptides with STOP or less than 10 amino acids
            if len(peptide_sequence) > 9 and not re.search("(STOP).*", peptide_sequence):
                protein_analysis = ProteinAnalysis(peptide_sequence)
                
                protein_features = [
                    line_components[0],
                    f"{protein_analysis.aromaticity():.4f}",
                    f"{protein_analysis.instability_index():.4f}",
                    f"{protein_analysis.isoelectric_point():.4f}",
                    f"{protein_analysis.gravy():.4f}",
                    peptide_sequence
                ]
                
                output_file.write('\t'.join(protein_features) + '\n')

# Main execution
input_filepath = '../data/intermediate/v37.LncRNA.translated.ORFs.txt'
output_filepath = '../data/intermediate/lncrna_peptides_features.txt'
extract_lncrna_peptides_features(input_filepath, output_filepath)
