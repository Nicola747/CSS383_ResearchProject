from species import Species
from variables import Variables

# Build the set/array of species for this study
organisms = [
    Species(name="zig-zag eel", latin_name="Mastacembelus armatus", file_name="zig_zag_eel"),
    Species(name="Rio pearlfish", latin_name="Nematolebias whitei", file_name="rio_pearlfish"),
    Species(name="pikepearch", latin_name="Sander lucioperca", file_name="pikeperch"),
    Species(name="jeweled blenny", latin_name="Salarias fasciatus", file_name="jeweled_blenny"),
    Species(name="Atlantic salmon", latin_name="Salmo salar", file_name="atlantic_salmon"),
    Species(name="channel catfish", latin_name="Ictalurus punctatus", file_name="channel_catfish"),
    Species(name="goldfish", latin_name="Carassius auratus", file_name="goldfish"),
    Species(name="Indian glassy fish", latin_name="Parambassis ranga", file_name="indian_glassy_fish"),
    Species(name="mummichog", latin_name="Fundulus heteroclitus", file_name="mummichog"),
    Species(name="Nile tilapia", latin_name="Oreochromis niloticus", file_name="nile_tilapia"),
    Species(name="orbiculate cardinalfish", latin_name="Sphaeramia orbicularis", file_name="orbiculate_cardinalfish"),
    Species(name="Japanese medaka", latin_name="Oryzias latipes", file_name="japanese_medaka"),
    Species(name="live sharksucker", latin_name="Echeneis naucrates", file_name="live_sharksucker"),
    Species(name="rainbow trout", latin_name="Oncorhynchus mykiss", file_name="rainbow_trout"),
    Species(name="pinecone soldierfish", latin_name="Myripristis murdjan", file_name="pinecone_soldierfish"),
    Species(name="tongue sole", latin_name="Cynoglossus semilaevis", file_name="tongue_sole"),
    Species(name="wolf-eel", latin_name="Anarrhichthys ocellatus", file_name="wolf_eel"),
    Species(name="electric eel", latin_name="Electrophorus electricus", file_name="electric_eel"),
    Species(name="sockeye salmon", latin_name="Oncorhynchus nerka", file_name="sockeye_salmon"),
    Species(name="guppy", latin_name="Poecilia reticulata", file_name="guppy"),
    Species(name="clown anemonefish", latin_name="Amphiprion ocellaris", file_name="clown_anemonefish"),
    Species(name="tiger tale seahorse", latin_name="Hippocampus comes", file_name="tiger_tail_seahorse")
]

#    - Compute DNA and Protein alignments for every pair of consecutive species and 
#        Store percentages in two arrays: DNA_IDENTITY and PROT_IDENTITY (to represent Y-axis values)
#    - Compute ABS(CUB1 - CUB2) for every pair of two consecutive species
#        Store CAI/CUB differences in array: CUB_DIFF (to represent X-axis values)

variables = []

from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mode = "global"

for index in range(0, len(organisms), 2):

    cub_diff = abs(organisms[index].CUB - organisms[index+1].CUB)

    org1 = organisms[index]
    org2 = organisms[index+1]

    dna_seq1 = org1.DNASequence
    dna_seq2 = org2.DNASequence
    alignments = aligner.align(dna_seq1, dna_seq2)

    dna = alignments[0].score
    dna_percentage = dna/len(dna_seq1) * 100.0

    protein_seq1 = org1.ProteinSequence
    protein_seq2 = org2.ProteinSequence
    alignments = aligner.align(protein_seq1, protein_seq2)

    prot = alignments[0].score
    prot_percentage = prot / len(protein_seq1) * 100.0
 
    variables.append(Variables(org1.Name, org2.Name, cub_diff, dna_percentage, prot_percentage))

# Sort collection of variables by their X (CUB/CAI) values for plotting purposes
variables.sort(key=lambda x: x.cub_diff, reverse=False)

CUB_DIFF = []
DNA_IDENTITY = []
PROT_IDENTITY = []

print("")

for index in range(0, len(variables), 1):
    var = variables[index]

    CUB_DIFF.append(var.cub_diff)
    DNA_IDENTITY.append(var.dna_identity)
    PROT_IDENTITY.append(var.protein_identity)

    print("Species: " + var.species1 + " and " + var.species2)
    print("CUB Difference: " + "{:.2f}".format(var.cub_diff))
    print("DNA Identity: " + "{:.2f}%".format(var.dna_identity) + " PROTEIN Identity: " + "{:.2f}%".format(var.protein_identity))
    print("")

# plot the curves
import matplotlib.pyplot as graph

graph.plot(CUB_DIFF, PROT_IDENTITY, 'ro-', label='Protein Alignment')
graph.plot(CUB_DIFF, DNA_IDENTITY, 'bs:', label='DNA Alignment')

graph.title('DNA/Protein Alignment Identity Percentages vs. CUB Differences Evolution')
graph.xlabel("Differences of CUB/CAI Values per Pair of Species")
graph.ylabel("Identity Percentages per Pair of Species")
graph.legend()

graph.show()