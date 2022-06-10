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
    Species(name="tiger tale seahorse", latin_name="Hippocampus comes", file_name="tiger_tail_seahorse"),
    
    Species(name="spotted gar", latin_name="Lepisosteus oculatus", file_name="spotted_gar"),
    Species(name="northern pike", latin_name="Esox lucius", file_name="northern_pike"),
    Species(name="Amazon molly", latin_name="Poecilia formosa", file_name="amazon_molly"),
    Species(name="Asian bonytongue", latin_name="Scleropages formosus", file_name="asian_bonytongue"),
    Species(name="Atlantic cod", latin_name="Gadus morhua", file_name="atlantic_cod"),
    Species(name="ballan wrasse", latin_name="Labrus bergylta", file_name="ballan_wrasse"),
    Species(name="barramundi perch", latin_name="Lates calcarifer", file_name="barramundi_perch"),
    Species(name="bicolor damselfish", latin_name="Stegastes partitus", file_name="bicolor_damselfish"),
    Species(name="Burton's mouthbrooder", latin_name="Haplochromis burtoni", file_name="burtons_mouthbrooder"),
    Species(name="climbing perch", latin_name="Anabas testudineus", file_name="climbing_perch"),
    Species(name="common carp", latin_name="Cyprinus carpio", file_name="common_carp"),
    Species(name="denticle herring", latin_name="Denticeps clupeoides", file_name="denticle_herring"),
    Species(name="eastern happy", latin_name="Astatotilapia calliptera", file_name="eastern_happy"),
    Species(name="flier cichlid", latin_name="Archocentrus centrarchus", file_name="flier_cichlid"),
    Species(name="greater amberjack", latin_name="Seriola dumerili", file_name="greater_amberjack"),
    Species(name="green swordtail", latin_name="Xiphophorus hellerii", file_name="green_swordtail"),
    Species(name="Japanese flounder", latin_name="Paralichthys olivaceus", file_name="japanese_flounder"),
    Species(name="Mexican tetra", latin_name="Astyanax mexicanus", file_name="mexican_tetra"),
    Species(name="Monterrey platyfish", latin_name="Xiphophorus couchianus", file_name="monterrey_platyfish"),
    Species(name="orangethroat darter", latin_name="Etheostoma spectabile", file_name="orangethroat_darter"),
    Species(name="red-bellied piranha", latin_name="Pygocentrus nattereri", file_name="red_bellied_piranha"),
    Species(name="reedfish", latin_name="Erpetoichthys calabaricus", file_name="reedfish"),
    Species(name="sailfin molly", latin_name="Poecilia latipinna", file_name="sailfin_molly"),
    Species(name="Siamese fighting fish", latin_name="Betta splendens", file_name="siamese_fighting_fish"),
    Species(name="southern platyfish", latin_name="Xiphophorus maculatus", file_name="southern_platyfish"),
    Species(name="spiny chromis", latin_name="Acanthochromis polyacanthus", file_name="spiny_chromis"),
    Species(name="torafugu", latin_name="Takifugu rubripes", file_name="torafugu"),
    Species(name="turqoise killfish", latin_name="Nothobranchius furzeri", file_name="turqoise_killfish"),
    Species(name="yellow catfish", latin_name="Tachysurus fulvidraco", file_name="yellow_catfish"),
    Species(name="zebra mbuna", latin_name="Maylandia zebra", file_name="zebra_mbuna")
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

#
# CORRELATION ANALYSIS
#

# scatter plot the DNA_IDENTITY and PROT_IDENTITY data sets in order to explore their correlation
graph.scatter(DNA_IDENTITY, PROT_IDENTITY)
graph.show()

# compute a couple of Correlation Coefficients to support conclusions from scatter plot above

# calculate Pearson's correlation
from scipy.stats import pearsonr
pearson_correlation, _ = pearsonr(DNA_IDENTITY, PROT_IDENTITY)
print('Pearson\'s correlation: %.3f' % pearson_correlation)

# calculate Spearman's correlation coefficient
from scipy.stats import spearmanr
spearman_correlation, _ = spearmanr(DNA_IDENTITY, PROT_IDENTITY)
print('Spearman\'s correlation: %.3f' % spearman_correlation)
print("")