class Variables():
    # The X, Y, Z variables to be plotted
    # X: CUB DIFFERENCE
    # Y: DNA IDENTITY PERCENTAGE
    # Z: PROTEIN IDENTITY PERCENTAGE

    def __init__(self, species1, species2, cub_diff, dna_identity, protein_identity):
        self.species1 = species1
        self.species2 = species2
        self.cub_diff = cub_diff
        self.dna_identity = dna_identity
        self.protein_identity = protein_identity