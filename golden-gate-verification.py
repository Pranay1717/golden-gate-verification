from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import AlignIO
import tempfile
import os
import subprocess
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
from Bio.Align import MultipleSeqAlignment
import pandas as pd
import streamlit as st

def construct_file_paths(gb_files):
    """Construct full file paths for the GenBank files based on their corresponding folders, checking for prefixed variants."""
    base_folder = "data"  # Base folder for all data
    folders = [
        "A.__Backbone_-_2024-09-20",
        "B.__Regulated_gene_-_2024-09-20",
        "C.__Regulated_promoter_-_2024-09-20",
        "D._Expression_promoters_-_2024-09-20",
        "E.__Expression_genes_-_2024-09-20"
    ]
    
    full_paths = []
    st.write("Constructing plasmid from the following components: ")
    for i, file in enumerate(gb_files):
        folder = folders[i]
        full_path = check_file_exists(os.path.join(base_folder, folder), file, i + 1)  # Include "data" folder in the path
        full_paths.append(full_path)
    
    st.write()
    return full_paths


def organize_assembled_plasmids(file_names):
    """Modify each file name to be in the 'Assembled Plasmids' folder."""
    folder_name = "Assembled Plasmids"
    new_file_paths = [os.path.join(folder_name, file_name) for file_name in file_names]
    return new_file_paths

def organize_sequencing_results(file_names):
    """Modify each file name to be in the 'Assembled Plasmids' folder."""
    folder_name = "QRMK5Y_genbank-files"
    new_file_paths = [os.path.join(folder_name, file_name) for file_name in file_names]
    return new_file_paths

def check_file_exists(folder, filename, index):
    """Check if the file exists in the folder, accounting for possible prefixes."""
    original_path = os.path.join(folder, filename)
    if os.path.exists(original_path):
        return original_path
    prefixed_name_one_space = f"{index}. {filename}"
    prefixed_name_two_space = f"{index}.  {filename}"
    prefixed_name_three_space = f"{index}.   {filename}"
    prefixed_name_tab = f"{index}.\t{filename}"
    
    prefixed_path_one_space = os.path.join(folder, prefixed_name_one_space)
    if os.path.exists(prefixed_path_one_space):
        st.write(prefixed_path_one_space)
        return prefixed_path_one_space
    
    prefixed_path_two_space = os.path.join(folder, prefixed_name_two_space)
    if os.path.exists(prefixed_path_two_space):
        st.write(prefixed_path_two_space)
        return prefixed_path_two_space
    
    prefixed_path_three_space = os.path.join(folder, prefixed_name_three_space)
    if os.path.exists(prefixed_path_three_space):
        st.write(prefixed_path_three_space)
        return prefixed_path_three_space
    
    prefixed_path_tab = os.path.join(folder, prefixed_name_tab)
    if os.path.exists(prefixed_path_tab):
        st.write(prefixed_path_tab)
        return prefixed_path_tab
    
    return original_path

def extract_sequence_from_gb(file_path):
    """Extract the sequence from a GenBank file."""
    record = SeqIO.read(file_path, "genbank")
    return record.seq

def complement_base(base):
    """Return the complement of a single base."""
    base = base.upper()
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement.get(base, base)  # Return the complement, or the base itself if not found.

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return ''.join(complement_base(base) for base in reversed(seq))

def match_overhangs(seq1, seq2, overhang_length):
    """Check if the overhangs between seq1 and seq2 match, ignoring case."""
    return seq1[-overhang_length:].upper() == seq2[:overhang_length].upper()

def find_single_occurrence(sequence, subseq):
    """Find a single occurrence of subseq in sequence, raise error if not exactly one."""
    positions = [i for i in range(len(sequence)) if sequence[i:i+len(subseq)] == subseq]
    if len(positions) != 1:
        raise ValueError(f"There must be exactly one occurrence of {subseq}. Found {len(positions)}.")
    return positions[0]

def bsa_digest(dna_sequence, overhang_length):
    """Perform a BsaI-style digest on a circular DNA sequence."""
    dna_sequence = dna_sequence.upper()

    gagacc = "GAGACC"
    ggtctc = "GGTCTC"

    gagacc_pos = find_single_occurrence(dna_sequence, gagacc)
    ggtctc_pos = find_single_occurrence(dna_sequence, ggtctc)

    start_pos = (gagacc_pos - 1) % len(dna_sequence) #1 - overhang_length) % len(dna_sequence)
    end_pos = (ggtctc_pos + len(ggtctc) + 2) % len(dna_sequence)
    
    if start_pos > end_pos: #make sure not to include the Bsa1 recognition sites
        extracted_seq = dna_sequence[end_pos-1:start_pos]
    else:
        extracted_seq = dna_sequence[end_pos-1:] + dna_sequence[:start_pos]

    return extracted_seq

def reverse_complement_if_needed(seq1, seq2, overhang_length):
    """Reverse complement seq2 if necessary based on overhangs."""
    rev_comp = False
    if match_overhangs(seq1, seq2, overhang_length):
        return rev_comp, seq2
    elif match_overhangs(seq1, reverse_complement(seq2), overhang_length):
        return rev_comp, reverse_complement(seq2)
    elif match_overhangs(reverse_complement(seq1), seq2, overhang_length):
        rev_comp = True
        return rev_comp, seq2
    elif match_overhangs(reverse_complement(seq1), reverse_complement(seq2), overhang_length):
        rev_comp = True
        return rev_comp, reverse_complement(seq2)
    else:
        raise ValueError("Overhangs do not match. Manual assembly required.")

def assemble_plasmid(gb_files, overhang_length=2):
    """Assemble a plasmid from a list of GenBank files."""
    sequences = [bsa_digest(extract_sequence_from_gb(file), overhang_length) for file in gb_files]
    
    plasmid_sequence = sequences[0]
    
    st.write("Backbone: " + str(gb_files[0]))
    st.write("Plasmid is now length " + str(len(plasmid_sequence)))
    
    for i in range(1, len(sequences)):
        st.write("Ligating " + str(gb_files[i]))
        rev_comp, next_seq = reverse_complement_if_needed(plasmid_sequence, sequences[i], overhang_length)
        if rev_comp:
            plasmid_sequence = reverse_complement(plasmid_sequence).copy()
        plasmid_sequence += next_seq[overhang_length:]
        st.write("Success! Plasmid is now length " + str(len(plasmid_sequence)))
    
    return plasmid_sequence

def save_plasmid_to_fasta(plasmid_sequence, output_file):
    """Save the plasmid sequence to a FASTA file and provide a download link in Streamlit."""
    output_directory = "Assembled Plasmids"
    
    # Construct the full path for the output file
    full_output_path = os.path.join(output_directory, output_file)
    #st.write(output_file)
    #st.write(full_output_path)
    # Save the FASTA content to the output directory
    record = SeqRecord(Seq(plasmid_sequence), id="plasmid", description="Assembled plasmid")
    SeqIO.write(record, output_file, "fasta")  # Save directly to full_output_path

    # Provide a download button in Streamlit
    with open(output_file, "rb") as file:
        st.download_button(
            label="Download Assembled Plasmid (FASTA)",
            data=file,
            file_name=output_file,
            mime="text/fasta"
        )
    
    st.write("Assembled plasmid saved to: " + str(output_file))

def update_fasta_with_reverse_complement(fasta_seq):
    """Update the sequence in a SeqRecord to its reverse complement."""
    fasta_seq.seq = fasta_seq.seq.reverse_complement()

def visualize_alignment(alignment):
    """Visualize the alignment using matplotlib."""
    seq1 = str(alignment[0].seq)
    seq2 = str(alignment[1].seq)    
    mismatch_positions = []
    for i, (base1, base2) in enumerate(zip(seq1, seq2)):
        if base1 != base2:
            mismatch_positions.append(i)

    plt.figure(figsize=(20, 2))
    
    for pos in mismatch_positions:
        plt.plot([pos, pos], [0, 1], color='red', linewidth=0.5, alpha=1)  # Vertical line for mismatch

    plt.xlim(0, len(seq1) - 1)
    plt.ylim(0, 1)
    plt.xlabel('Position in Alignment')
    plt.title('Mismatch Visualization: Red = Mismatch')
    plt.xticks([])
    plt.yticks([])
    st.pyplot(plt)
    plt.clf()
    #plt.show()
    
def read_fasta(file_path):
    """Read a FASTA file and return the sequence as a string."""
    sequence = ""
    
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            # Skip header lines (lines starting with '>')
            if not line.startswith(">"):
                sequence += line  # Concatenate sequence lines
    
    return sequence
    
def read_fasta_from_upload(uploaded_file):
    """Read a FASTA file (from an UploadedFile) and return the sequence as a string."""
    sequence = ""
    
    # Use the file-like object directly without opening it
    for line in uploaded_file:
        line = line.decode("utf-8").strip()  # Decode bytes to string
        # Skip header lines (lines starting with '>')
        if not line.startswith(">"):
            sequence += line  # Concatenate sequence lines
    
    return sequence
    
def run_mafft(fasta_file, gb_file, rc=0):
    """Run MAFFT alignment between sequences from a FASTA file and a GenBank file."""
    
    # Read sequences from the provided files
    fasta_seq = read_fasta(fasta_file) #SeqIO.read(fasta_file, "fasta")
    gb_seq = read_fasta_from_upload(gb_file) #SeqIO.read(gb_file, "fasta")
    
    # Update fasta_seq if reverse complementation is needed
    if rc:
        update_fasta_with_reverse_complement(fasta_seq)
    
    # Display lengths of the sequences
    st.write("Template (from .fasta) Length: " + str(len(fasta_seq)))
    st.write("Sequence (from .fasta) Length: " + str(len(gb_seq)))

    # Optionally, extend the sequences (ensure this logic is as intended)
    fasta_seq = fasta_seq + fasta_seq
    gb_seq = gb_seq + gb_seq

    # Create a temporary file for input to MAFFT
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as temp_input:
        # Convert fasta_seq string to SeqRecord
        fasta_seq_record = SeqRecord(Seq(fasta_seq), id="fasta_seq", description="FASTA sequence")

        # Ensure gb_seq is a SeqRecord if it's a string
        if isinstance(gb_seq, str):
            gb_seq_record = SeqRecord(Seq(gb_seq), id="gb_seq", description="GenBank sequence")
        else:
            gb_seq_record = gb_seq  # Assuming gb_seq is already a SeqRecord

        # Write the SeqRecord objects to the temporary file
        SeqIO.write([fasta_seq_record, gb_seq_record], temp_input, "fasta")
        input_file = temp_input.name

    # Prepare the output file name
    output_file = input_file + "_aligned.fasta"

    try:
        # Run the MAFFT command
        subprocess.run(
            ["mafft", "--auto", "--adjustdirection", input_file],
            stdout=open(output_file, 'w'),
            stderr=subprocess.PIPE,
            check=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"MAFFT failed: {e.stderr.decode()}")
    
    # Read the resulting alignment
    alignment = AlignIO.read(output_file, "fasta")
    
    # Initialize variables to find the best alignment
    aligned_seq1 = str(alignment[0].seq)
    aligned_seq2 = str(alignment[1].seq)

    original_length = len(fasta_seq) // 2
    best_start = 0
    best_matches = 0

    # Find the best matching region
    for i in range(len(aligned_seq1) - original_length + 1):
        window_seq1 = aligned_seq1[i:i + original_length]
        window_seq2 = aligned_seq2[i:i + original_length]
        matches = sum(1 for base1, base2 in zip(window_seq1, window_seq2) if base1 == base2)
        
        if matches > best_matches:
            best_matches = matches
            best_start = i

    # Extract the best aligned sequences
    best_aligned_seq1 = aligned_seq1[best_start:best_start + original_length]
    best_aligned_seq2 = aligned_seq2[best_start:best_start + original_length]

    st.write("Determined optimal reverse complementation + rotation for alignment.")
    
    # Create a new MultipleSeqAlignment object for the best alignment
    best_alignment = MultipleSeqAlignment([
        SeqRecord(Seq(best_aligned_seq1), id="fasta_seq"),
        SeqRecord(Seq(best_aligned_seq2), id="gb_seq")
    ])

    # Visualize the best alignment
    visualize_alignment(best_alignment)
    
    # Calculate similarity score and differences
    sim_score, differences = calculate_similarity_and_differences(best_alignment)
    
    # Clean up temporary files
    os.remove(input_file)
    os.remove(output_file)
    
    return sim_score, differences

def calculate_similarity_and_differences(alignment):
    """Calculate similarity score and identify locations of differences in the alignment."""
    aligned_seq1 = alignment[0].seq
    aligned_seq2 = alignment[1].seq
    
    total_bases = len(aligned_seq1)
    matches = 0
    differences = []
    
    for i in range(total_bases):
        base1 = aligned_seq1[i]
        base2 = aligned_seq2[i]
        if base1 == base2:
            matches += 1
        else:
            # Record the position (1-based) and the differing bases
            differences.append((i + 1, base1, base2))
    
    similarity_score = (matches / total_bases) * 100
    
    return similarity_score, differences

def get_alignment(fasta_file, gb_file):
    
    st.write("\nRunning MAFFT Alignment")
    
    sim_score, differences = run_mafft(fasta_file, gb_file)
    
    st.write(f"Similarity score: {sim_score:.2f}%")
    st.write("Differences:")
    if len(differences) < 5:
        for diff in differences:
            st.write(f"Position {diff[0]}: {diff[1]} vs {diff[2]}")
    else:
        st.write("High quantity of differences (omitted)")
        
    if sim_score < 50:
        st.write("Low similarity score, trying reverse complement.")
        sim_score, differences = run_mafft(fasta_file, gb_file, rc=1)
        st.write(f"Similarity score: {sim_score:.2f}%")
        st.write("Differences:")
        if len(differences) < 5:
            for diff in differences:
                st.write(f"Position {diff[0]}: {diff[1]} vs {diff[2]}")
        else:
            st.write("High quantity of differences (omitted)")

st.title("Plasmid Assembly and Alignment for Golden Gate Verification")

#input1 = st.text_input("Sequencing Results")
input0 = st.file_uploader("Upload Sequencing Results .fasta File", type=["fasta"])
sequencing_results = [input0]

input1 = st.text_input("Sequencing Results Filename")
input2 = st.text_input("Part A")
input3 = st.text_input("Part B")
input4 = st.text_input("Part C")
input5 = st.text_input("Part D")
input6 = st.text_input("Part E")

if st.button("Submit"):
    os.makedirs("Assembled Plasmids", exist_ok=True)
    row = [input1,input2,input3,input4,input5,input6]

    #sequencing_results = [row[0].replace('.fasta', '.gbk')]
    gb_files = row[1:6]
    plasmid_names = [input1.split('_', 1)[1]]

    plasmid_names = organize_assembled_plasmids(plasmid_names)
    #sequencing_results = organize_sequencing_results(sequencing_results)
    gb_files = construct_file_paths(gb_files)

    assembled_plasmid = assemble_plasmid(gb_files, overhang_length=4)
    save_plasmid_to_fasta(assembled_plasmid, plasmid_names[0])

    for i in range(len(plasmid_names)):
        get_alignment(plasmid_names[i], sequencing_results[i])

    st.success("Processing completed!")
