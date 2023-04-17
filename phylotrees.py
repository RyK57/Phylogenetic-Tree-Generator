import streamlit as st
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

def read_fasta(file):
    # Read the contents of the file
    contents = file.getvalue().decode("utf-8")
    # Split the contents into lines
    lines = contents.split("\n")
    # Remove any empty lines
    lines = [line for line in lines if line.strip()]
    # Parse the FASTA sequences
    sequences = {}
    max_length = 0
    for line in lines:
        if line.startswith(">"):
            current_sequence = line[1:]
            sequences[current_sequence] = ""
        else:
            sequences[current_sequence] += line

            # Keep track of the maximum sequence length
            if len(sequences[current_sequence]) > max_length:
                max_length = len(sequences[current_sequence])

    # Add padding to ensure all sequences are the same length
    for sequence_name in sequences:
        sequence = sequences[sequence_name]
        padding = max_length - len(sequence)
        sequences[sequence_name] = sequence + "-" * padding

    # Convert the sequences dictionary to a list of SeqRecord objects
    seq_records = [SeqRecord(Seq(seq), id=name) for name, seq in sequences.items()]
    # Align the sequences and return the result
    return align_sequences(seq_records)


def align_sequences(seq_records):
    # Use the first sequence as the template for the alignment
    template_seq = seq_records[0].seq
    # Create a list of sequences to be aligned
    seqs_to_align = [str(seq_record.seq) for seq_record in seq_records]
    # Align the sequences
    alignment = MultipleSeqAlignment(seq_records)
    # Return the alignment
    return alignment

def build_phylogenetic_tree(alignment):
    # Create a distance calculator object
    calculator = DistanceCalculator('identity')
    # Calculate the distances between sequences in the alignment
    distance_matrix = calculator.get_distance(alignment)
    # Create a distance-based tree constructor object
    constructor = DistanceTreeConstructor()
    # Build the tree
    tree = constructor.upgma(distance_matrix)
    return tree

# Set page title
st.set_page_config(page_title="FASTA File Upload")

# Create upload button
uploaded_file = st.file_uploader(
    "Upload a text file containing FASTA sequences", type=["txt"]
)

file_url = 'https://drive.google.com/file/d/1vwHW4NcsU6tfJPbKVr1CGiZLtz5NpgtY/view?usp=sharing'

# Create a link to download the file
st.markdown(f'<a href="https://drive.google.com/file/d/1vwHW4NcsU6tfJPbKVr1CGiZLtz5NpgtY/view?usp=sharing" download>Example Input file (M19961 Accession Number NCBI)</a>', unsafe_allow_html=True)

# If file is uploaded, read and display alignment
if uploaded_file is not None:
    alignment = read_fasta(uploaded_file)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    st.write(f"Consensus sequence:\n{consensus}")
    st.write("Aligned sequences:")
    # Set maximum alignment length to display
    MAX_ALIGNMENT_LENGTH = 1000

    # Print pairwise alignment for each sequence
    for seq_record in alignment:
        st.write(seq_record.id)
        st.write(seq_record.seq)
        st.write("")

        # Align the sequence with the consensus sequence
        pairwise_alignments = align_sequences([seq_record, SeqRecord(Seq(consensus), id="Consensus")])

        # Print pairwise alignment for each pair of sequences
        for a in pairwise_alignments:
            # Check if alignment is shorter than the maximum length
            if len(a[0]) <= MAX_ALIGNMENT_LENGTH and len(a[1]) <= MAX_ALIGNMENT_LENGTH and len(a[2]) <= MAX_ALIGNMENT_LENGTH:
                alignments = pairwise2.align.globalxx(a[0], a[1])
                alignment_str = format_alignment(*alignments[0], full_sequences=False)
                st.write(alignment_str)

    # Calculate distance matrix
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    st.write("Distance matrix:")
    st.write(distance_matrix)

    # Construct the tree
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.build_tree(alignment)
    st.write("Phylogenetic tree:")
    st.write(tree)

    # Draw the tree
    Phylo.draw(tree)






## Reads and Aligns sequences codes

# import streamlit as st
# from Bio.Align import MultipleSeqAlignment, AlignInfo
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
# from Bio import pairwise2
# from Bio.pairwise2 import format_alignment

# def read_fasta(file):
#     # Read the contents of the file
#     contents = file.getvalue().decode("utf-8")
#     # Split the contents into lines
#     lines = contents.split("\n")
#     # Remove any empty lines
#     lines = [line for line in lines if line.strip()]
#     # Parse the FASTA sequences
#     sequences = {}
#     max_length = 0
#     for line in lines:
#         if line.startswith(">"):
#             current_sequence = line[1:]
#             sequences[current_sequence] = ""
#         else:
#             sequences[current_sequence] += line

#             # Keep track of the maximum sequence length
#             if len(sequences[current_sequence]) > max_length:
#                 max_length = len(sequences[current_sequence])

#     # Add padding to ensure all sequences are the same length
#     for sequence_name in sequences:
#         sequence = sequences[sequence_name]
#         padding = max_length - len(sequence)
#         sequences[sequence_name] = sequence + "-" * padding

#     # Convert the sequences dictionary to a list of SeqRecord objects
#     seq_records = [SeqRecord(Seq(seq), id=name) for name, seq in sequences.items()]
#     # Align the sequences and return the result
#     return align_sequences(seq_records)


# def align_sequences(seq_records):
#     # Use the first sequence as the template for the alignment
#     template_seq = seq_records[0].seq
#     # Create a list of sequences to be aligned
#     seqs_to_align = [str(seq_record.seq) for seq_record in seq_records]
#     # Align the sequences
#     alignment = MultipleSeqAlignment(seq_records)
#     # Return the alignment
#     return alignment


# # Set page title
# st.set_page_config(page_title="FASTA File Upload")

# # Create upload button
# uploaded_file = st.file_uploader(
#     "Upload a text file containing FASTA sequences", type=["txt"]
# )

# file_url = 'https://drive.google.com/file/d/1vwHW4NcsU6tfJPbKVr1CGiZLtz5NpgtY/view?usp=sharing'

# # Create a link to download the file
# st.markdown(f'<a href="https://drive.google.com/file/d/1vwHW4NcsU6tfJPbKVr1CGiZLtz5NpgtY/view?usp=sharing" download>Example Input file (M19961 Accession Number NCBI)</a>', unsafe_allow_html=True)

# # If file is uploaded, read and display alignment
# if uploaded_file is not None:
#     alignment = read_fasta(uploaded_file)
#     summary_align = AlignInfo.SummaryInfo(alignment)
#     consensus = summary_align.dumb_consensus()
#     st.write(f"Consensus sequence:\n{consensus}")
#     st.write("Aligned sequences:")
#    # Set maximum alignment length to display
#     MAX_ALIGNMENT_LENGTH = 1000

#         # Print pairwise alignment for each sequence
#     for seq_record in alignment:
#         st.write(seq_record.id)
#         st.write(seq_record.seq)
#         st.write("")

#         # Align the sequence with the consensus sequence
#         pairwise_alignments = align_sequences([seq_record, SeqRecord(Seq(consensus), id="Consensus")])

#         # Print pairwise alignment for each pair of sequences
#         for a in pairwise_alignments:
#             # Check if alignment is shorter than the maximum length
#             if len(a[0]) <= MAX_ALIGNMENT_LENGTH and len(a[1]) <= MAX_ALIGNMENT_LENGTH and len(a[2]) <= MAX_ALIGNMENT_LENGTH:
#                 alignments = pairwise2.align.globalxx(a[0], a[1])
#                 alignment_str = format_alignment(*alignments[0], full_sequences=False)
#                 st.write(alignment_str)
