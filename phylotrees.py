import streamlit as st
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import matplotlib


matplotlib.use('Agg')

# Primary accent for interactive elements
primaryColor = '#7792E3'

# Background color for the main content area
backgroundColor = '#273346'

# Background color for sidebar and most interactive widgets
secondaryBackgroundColor = '#B9F1C0'

# Color used for almost all text
textColor = '#FFFFFF'

# Font family for all text in the app, except code blocks
# Accepted values (serif | sans serif | monospace) 
# Default: "sans serif"
font = "sans serif"

def read_fasta(file):
    # Reading the contents of the file
    contents = file.getvalue().decode("utf-8")
    # Splitting the contents into lines
    lines = contents.split("\n")
    # Removing any empty lines
    lines = [line for line in lines if line.strip()]
    # Parsing the FASTA sequences
    sequences = {}
    max_length = 0
    for line in lines:
        if line.startswith(">"):
            current_sequence = line[1:]
            sequences[current_sequence] = ""
        else:
            sequences[current_sequence] += line

            # Keeping track of the maximum sequence length
            if len(sequences[current_sequence]) > max_length:
                max_length = len(sequences[current_sequence])

    # Adding padding to ensure all sequences are the same length
    for sequence_name in sequences:
        sequence = sequences[sequence_name]
        padding = max_length - len(sequence)
        sequences[sequence_name] = sequence + "-" * padding

    # Converting the sequences dictionary to a list of SeqRecord objects
    seq_records = [SeqRecord(Seq(seq), id=name) for name, seq in sequences.items()]
    # Aligning the sequences and return the result
    return align_sequences(seq_records)


def align_sequences(seq_records):
    # Using the first sequence as the template for the alignment
    template_seq = seq_records[0].seq
    # Creating a list of sequences to be aligned
    seqs_to_align = [str(seq_record.seq) for seq_record in seq_records]
    # Aligning the sequences
    alignment = MultipleSeqAlignment(seq_records)
    # Returning the alignment
    return alignment

def build_phylogenetic_tree(alignment):
    # Creating a distance calculator object
    calculator = DistanceCalculator('identity')
    # Calculating the distances between sequences in the alignment
    distance_matrix = calculator.get_distance(alignment)
    # Creating a distance-based tree constructor object
    constructor = DistanceTreeConstructor()
    # Buildng the tree
    tree = constructor.upgma(distance_matrix)
    return tree

# # Setting page title
# st.set_page_config(page_title="FASTA File Upload")

# Upload button
uploaded_file = st.file_uploader(
    "Upload a text file containing FASTA sequences", type=["txt"]
)

file_url = 'https://drive.google.com/file/d/1vwHW4NcsU6tfJPbKVr1CGiZLtz5NpgtY/view?usp=sharing'

# link to download the example sequence file
st.markdown(f'<a href="https://drive.google.com/file/d/1coCSpNrDI599WICcsCj5cs6C8U9Hb46g/view?usp=sharing" download>Example Input file (Carpodacus Mexicanus)</a>', unsafe_allow_html=True)

 st.title("Phylogenetic Tree Generator")
 st.write("-Rithvik Sabnekar")

    

# If file is uploaded, reading and displaying alignment
if uploaded_file is not None:
    alignment = read_fasta(uploaded_file)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    st.write(f"Consensus sequence:\n{consensus}")
    st.write("Aligned sequences:")
    # Setting maximum alignment length to display
    MAX_ALIGNMENT_LENGTH = 1000

    # Printing pairwise alignment for each sequence
    for seq_record in alignment:
        st.write(seq_record.id)
        st.write(seq_record.seq)
        st.write("")

        # Aligning the sequence with the consensus sequence
        pairwise_alignments = align_sequences([seq_record, SeqRecord(Seq(str(consensus)), id="Consensus")])


        # Printing pairwise alignment for each pair of sequences
        for a in pairwise_alignments:
            # Check if alignment is shorter than the maximum length
            if len(a[0]) <= MAX_ALIGNMENT_LENGTH and len(a[1]) <= MAX_ALIGNMENT_LENGTH and len(a[2]) <= MAX_ALIGNMENT_LENGTH:
                alignments = pairwise2.align.globalxx(a[0], a[1])
                alignment_str = format_alignment(*alignments[0], full_sequences=False)
                st.write(alignment_str)

    # Calculating distance matrix
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    st.write("Distance matrix:")
    st.write(distance_matrix)

    # Constructing the tree
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.build_tree(alignment)
    st.write("Phylogenetic tree:")
    st.write(tree)

    # Drawing the tree
    Phylo.draw(tree)
    
    plot = matplotlib.pyplot.gcf()
    # Displaying the plot in Streamlit ap
    st.pyplot(plot)
