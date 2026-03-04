import pandas as pd
from collections import defaultdict, Counter
import re

def parse_sample_from_header(header):
    """
    Extract sample ID from header string.
    Examples: 
    - 'Vfaba.Hedin2.R2.2g011475_F1_t30162' -> 'F1'
    - 'Vfaba.Hedin2.R2.2g011475_F246_t81320' -> 'F246'
    - 'Vfaba.Hedin2.R2.3g018667_F169_t44111' -> 'F169'
    """
    # Use regex to find pattern like _F followed by numbers
    match = re.search(r'_F(\d+)_', header)
    if match:
        return f'F{match.group(1)}'
    return None

def count_amino_acids(sequence):
    """
    Count each amino acid in the sequence.
    Returns a dictionary with amino acid counts.
    """
    # Standard 20 amino acids
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
                   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    
    # Count amino acids in sequence
    aa_count = Counter(sequence.upper())
    
    # Create dictionary with all amino acids (including zeros)
    result = {}
    for aa in amino_acids:
        result[aa] = aa_count.get(aa, 0)
    
    return result

def get_fileB_multipliers(data):
    """
    Process File B type data to get sequence ID multipliers.
    
    Args:
        data: List of dictionaries containing protein information
              Each dictionary should have keys: 'Sequence_IDs', 'Shared_Count'
    
    Returns:
        Dictionary mapping sequence_id to shared_count multiplier
    """
    sequence_multipliers = {}
    
    for entry in data:
        sequence_ids = entry['Sequence_IDs']
        shared_count = int(entry['Shared_Count'])
        
        # Split sequence IDs by comma and clean whitespace
        ids = [seq_id.strip() for seq_id in sequence_ids.split(',')]
        
        # Store multiplier for each sequence ID
        for seq_id in ids:
            sequence_multipliers[seq_id] = shared_count
    
    return sequence_multipliers

def process_fileA_with_multipliers(data, sequence_multipliers=None):
    """
    Process File A type data (with All_Headers column), applying multipliers from File B.
    
    Args:
        data: List of dictionaries containing protein information
              Each dictionary should have keys: 'All_Headers', 'Sequence'
        sequence_multipliers: Dictionary mapping sequence_id to shared_count multiplier
    
    Returns:
        Dictionary mapping sample_id to amino acid counts
    """
    sample_aa_counts = defaultdict(lambda: defaultdict(int))
    
    if sequence_multipliers is None:
        sequence_multipliers = {}
    
    for entry in data:
        sequence = entry['Sequence']
        all_headers = entry['All_Headers']
        
        # Split headers by semicolon and clean whitespace
        headers = [header.strip() for header in all_headers.split(';')]
        
        # Get amino acid counts for this sequence
        aa_counts = count_amino_acids(sequence)
        
        # Process each header to extract sample information
        for header in headers:
            sample_id = parse_sample_from_header(header)
            if sample_id:
                # Check if this sequence ID has a multiplier from File B
                multiplier = sequence_multipliers.get(header, 1)  # Default multiplier is 1
                
                # Add amino acid counts to the sample, multiplied by shared_count if available
                for aa, count in aa_counts.items():
                    sample_aa_counts[sample_id][aa] += count * multiplier
    
    return sample_aa_counts

def create_results_dataframe(sample_aa_counts):
    """
    Convert sample amino acid counts dictionary to DataFrame.
    
    Args:
        sample_aa_counts: Dictionary mapping sample_id to amino acid counts
    
    Returns:
        DataFrame with amino acid counts per sample
    """
    df_data = []
    for sample_id, aa_counts in sample_aa_counts.items():
        row = {'Sample': sample_id}
        row.update(aa_counts)
        df_data.append(row)
    
    # Create DataFrame and sort by sample name
    df = pd.DataFrame(df_data)
    if not df.empty:
        df = df.sort_values('Sample').reset_index(drop=True)
        
        # Ensure all amino acid columns are present (fill with 0 if missing)
        amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
                       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        for aa in amino_acids:
            if aa not in df.columns:
                df[aa] = 0
        
        # Reorder columns: Sample first, then amino acids in alphabetical order
        cols = ['Sample'] + sorted(amino_acids)
        df = df[cols]
    
    return df

def calculate_aa_counts_from_files(fileA_path=None, fileB_path=None):
    """
    Read data from CSV files and calculate amino acid counts.
    
    Args:
        fileA_path: Path to File A (with All_Headers column and Sequence)
        fileB_path: Path to File B (with Sequence_IDs and Shared_Count columns, no need for Sequence)
    
    Returns:
        DataFrame with amino acid counts per sample
    """
    sequence_multipliers = {}
    
    # Process File B to get multipliers
    if fileB_path:
        print(f"Processing File B for multipliers: {fileB_path}")
        dfB = pd.read_csv(fileB_path)
        dataB = dfB.to_dict('records')
        sequence_multipliers = get_fileB_multipliers(dataB)
        print(f"Found {len(sequence_multipliers)} sequence multipliers from File B")
    
    # Process File A with multipliers
    if fileA_path:
        print(f"Processing File A: {fileA_path}")
        dfA = pd.read_csv(fileA_path)
        dataA = dfA.to_dict('records')
        sample_counts = process_fileA_with_multipliers(dataA, sequence_multipliers)
    else:
        print("No File A provided!")
        return pd.DataFrame()
    
    # Create final DataFrame
    result_df = create_results_dataframe(sample_counts)
    return result_df

def example_usage():
    """
    Example of how to use the amino acid counter with sample data.
    """
    # Sample File A data (contains sequences)
    fileA_data = [
        {
            'Group_ID': 'Duplicate_Group_1',
            'Representative_Header': 'Vfaba.Hedin2.R2.2g011475_F1_t30162',
            'Sequence': 'RSDPENPFIFESNSFQTLFENDNGHIRLLQKFDERSKILENLQNYRLLEYKSKPRTIFLPQQTDADFILVVLSGKAILTVLKPNDRNSFNLERGDTIKLPAGTIAYLVNRDDNEDLRVLDLALPVNRPGQLQSFLLSGSQNQQSILSGFSKNILEASFNTGYEEIEKVLLEEHDKETQHRRSLRDKREHSQEEDVIVKLSTGQIEELSKNAKSSSKKSVSSESEPFNLRSRSPIYSNKFGRFFEITPEKNPQLQDLNVFVSSVEIKEGSLLLPHYNSRAIVIVAVNDGNGDFELVGQRNENQQGQRKEDDEEEEQGEEEINKQVKNYKAKLSRGDVFVIPAGYPVAIKASSNLDLLGFGINAKNNQRNFLAGEEDNVVIQIHQPVKELVFPGSAQEVDRLLENQKQSHFANAQPQHKERGSHETRDPLYSILDAF',
            'All_Headers': 'Vfaba.Hedin2.R2.2g011475_F1_t30162; Vfaba.Hedin2.R2.2g011475_F246_t81320; Vfaba.Hedin2.R2.2g011475_F318_t29235; Vfaba.Hedin2.R2.2g011475_F351_t49663; Vfaba.Hedin2.R2.2g011475_F402_t43182; Vfaba.Hedin2.R2.2g011475_F72_t57630'
        }
    ]
    
    # Sample File B data (only multipliers, no sequences needed)
    fileB_data = [
        {
            'Sequence_IDs': 'Vfaba.Hedin2.R2.2g011475_F1_t30162',
            'Is_Unique': 'Y',
            'Shared_Count': 2  # This sequence should be counted 2 times for F1
        },
        {
            'Sequence_IDs': 'Vfaba.Hedin2.R2.2g011475_F246_t81320,Vfaba.Hedin2.R2.2g011475_F318_t29235',
            'Is_Unique': 'N',
            'Shared_Count': 3  # These sequences should be counted 3 times each
        }
    ]
    
    print("Getting multipliers from File B...")
    sequence_multipliers = get_fileB_multipliers(fileB_data)
    print("Multipliers:", sequence_multipliers)
    
    print("\nProcessing File A with multipliers...")
    sample_counts = process_fileA_with_multipliers(fileA_data, sequence_multipliers)
    
    # Create final result
    result = create_results_dataframe(sample_counts)
    
    print("\nAmino acid counts per sample (with File B multipliers applied):")
    print("=" * 80)
    print(result.to_string(index=False))
    
    # Save to CSV
    result.to_csv('amino_acid_counts_by_sample.csv', index=False)
    print(f"\nResults saved to 'amino_acid_counts_by_sample.csv'")
    
    return result

# Main usage functions:
# For processing both files:
# result = calculate_aa_counts_from_files('fileA.csv', 'fileB.csv')

# For processing only File A:
# result = calculate_aa_counts_from_files(fileA_path='fileA.csv')

# Run example
#if __name__ == "__main__":
#    example_result = example_usage()

# 如果你想使用真实文件，取消下面的注释并修改文件名：
if __name__ == "__main__":
    #result = calculate_aa_counts_from_files('legumina_aligned_all_redundancy_analysis.csv', 'output_merged_sequences_comparison.csv')
    result = calculate_aa_counts_from_files('albumin_all_nonredundant_redundancy_analysis.csv', 'output_merged_sequences_comparison.csv')
    print(result)
    result.to_csv('albumin_amino_acid_results.csv', index=False)
    print("结果已保存到 final_results.csv")