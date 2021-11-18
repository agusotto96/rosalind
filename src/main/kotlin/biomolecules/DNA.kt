package biomolecules

class DNA internal constructor(nucleotides: List<DNANucleotide>) : NucleicAcid<DNA, DNANucleotide>(nucleotides) {

    fun transcript() = RNA(monomers.map(::findTranscript))

    override fun create(monomers: List<DNANucleotide>) = DNA(monomers)

    override fun findComplement(nucleotide: DNANucleotide) = when (nucleotide) {
        Adenine -> Thymine
        Cytosine -> Guanine
        Guanine -> Cytosine
        Thymine -> Adenine
    }

    private fun findTranscript(nucleotide: DNANucleotide) = when (nucleotide) {
        Adenine -> Adenine
        Cytosine -> Cytosine
        Guanine -> Guanine
        Thymine -> Uracil
    }

}

fun consensus(sequences: List<DNA>) = consensus(sequences, ::DNA)

fun DNA(symbols: String): DNA {
    val nucleotides = symbols.map {
        DNA_NUCLEOTIDES_BY_SYMBOL[it] ?: throw IllegalArgumentException("$it is not a valid DNA monomer symbol")
    }
    return DNA(nucleotides)
}

private val DNA_NUCLEOTIDES_BY_SYMBOL = DNANucleotide::class.sealedSubclasses
    .mapNotNull { it.objectInstance }
    .associateBy(DNANucleotide::symbol)