package biomolecules

class Protein internal constructor(aminoAcids: List<AminoAcid>) : Polymer<Protein, AminoAcid>(aminoAcids) {

    override fun create(monomers: List<AminoAcid>) = Protein(monomers)

    fun mass() = monomers.sumOf(AminoAcid::mass)

}

fun consensus(sequences: List<Protein>) = consensus(sequences, ::Protein)

fun Protein(symbols: String): Protein {
    val nucleotides = symbols.map {
        AMINO_ACIDS_BY_SYMBOL[it] ?: throw IllegalArgumentException("$it is not a valid Protein monomer symbol")
    }
    return Protein(nucleotides)
}

private val AMINO_ACIDS_BY_SYMBOL = AminoAcid.values()
    .associateBy(AminoAcid::symbol)
