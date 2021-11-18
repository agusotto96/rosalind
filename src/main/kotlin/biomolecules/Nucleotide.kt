package biomolecules

sealed interface Nucleotide : Monomer

sealed interface DNANucleotide : Nucleotide

sealed interface RNANucleotide : Nucleotide

object Adenine : DNANucleotide, RNANucleotide {
    override val symbol = 'A'
}

object Cytosine : DNANucleotide, RNANucleotide {
    override val symbol = 'C'
}

object Guanine : DNANucleotide, RNANucleotide {
    override val symbol = 'G'
}

object Thymine : DNANucleotide {
    override val symbol = 'T'
}

object Uracil : RNANucleotide {
    override val symbol = 'U'
}

internal val purines = setOf(Adenine, Guanine)

internal val pyrimidines = setOf(Cytosine, Thymine, Uracil)