from fastafileparser import FASTAFileParser

class Species():
    # Collection of properties that define a Species

    def __init__(self, name, latin_name, file_name):

        self.Name = name
        self.LatinName = latin_name

        import os
        dirname = os.path.dirname(__file__)

        self.GeneFilePath = os.path.join(dirname, "files\gene_" + file_name + ".fna")
        self.RNAFilePath = os.path.join(dirname, "files\\rna_" + file_name + ".fna")
        self.ProteinFilePath = os.path.join(dirname, "files\protein_" + file_name + ".faa")

        self.GeneFileParser = FASTAFileParser(self.GeneFilePath)
        self.ProteinFileParser = FASTAFileParser(self.ProteinFilePath)

        self.set_codon_usage_bias()
        self.set_dna_sequence()
        self.set_protein_sequence()

    def set_codon_usage_bias(self):

        """ Compute Codon Usage Bias / Codon Adaptation Index """
        from Bio.SeqUtils import CodonUsage as CU

        myIndex = CU.CodonAdaptationIndex()
        myIndex.generate_index(self.RNAFilePath)

        self.CUB = sum (myIndex.index.values())

    def set_dna_sequence(self):
        """ Load DNA sequence from DNA FASTA file """
        self.DNASequence = self.GeneFileParser.get_sequence()

    def set_protein_sequence(self):
        """ Load Protein sequence from DNA FASTA file """
        self.ProteinSequence = self.ProteinFileParser.get_sequence()