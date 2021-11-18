package biomolecules

class RNA internal constructor(nucleotides: List<RNANucleotide>) : NucleicAcid<RNA, RNANucleotide>(nucleotides) {

    fun reverseTranscript() = DNA(monomers.map(::findReverseTranscript))

    fun translate(): Set<Protein> {
        val candidates = mutableSetOf<MutableList<AminoAcid>>()
        val proteins = mutableSetOf<Protein>()
        val triplets = monomers.chunked(3)
        for (triplet in triplets) {
            val aminoAcid = codons[triplet]
            if (aminoAcid != null) {
                if (aminoAcid == AminoAcid.Methionine) {
                    candidates.add(mutableListOf())
                }
                candidates.forEach { it.add(aminoAcid) }
            } else {
                proteins.addAll(candidates.map(::Protein))
                candidates.clear()
            }
        }
        return proteins
    }

    override fun create(monomers: List<RNANucleotide>) = RNA(monomers)

    override fun findComplement(nucleotide: RNANucleotide) = when (nucleotide) {
        Adenine -> Uracil
        Cytosine -> Guanine
        Guanine -> Cytosine
        Uracil -> Adenine
    }

    private fun findReverseTranscript(nucleotide: RNANucleotide) = when (nucleotide) {
        Adenine -> Adenine
        Cytosine -> Cytosine
        Guanine -> Guanine
        Uracil -> Thymine
    }

}

fun consensus(sequences: List<RNA>) = consensus(sequences, ::RNA)

fun RNA(symbols: String): RNA {
    val nucleotides = symbols.map {
        RNA_NUCLEOTIDES_BY_SYMBOL[it] ?: throw IllegalArgumentException("$it is not a valid RNA monomer symbol")
    }
    return RNA(nucleotides)
}

private val RNA_NUCLEOTIDES_BY_SYMBOL = RNANucleotide::class.sealedSubclasses
    .mapNotNull { it.objectInstance }
    .associateBy(RNANucleotide::symbol)