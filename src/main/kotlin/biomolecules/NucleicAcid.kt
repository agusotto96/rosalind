package biomolecules

sealed class NucleicAcid<N : NucleicAcid<N, M>, M : Nucleotide>(nucleotides: List<M>) : Polymer<N, M>(nucleotides) {

    fun reverseComplement() = create(monomers.map(::findComplement).asReversed())

    fun gcContent() = monomers.count(::isGc) / (count() / 100.0)

    fun isPalindrome() = this == reverseComplement()

    protected abstract fun findComplement(nucleotide: M): M

    private fun isGc(nucleotide: Nucleotide) = nucleotide in setOf(Guanine, Cytosine)

    fun readingFrames(): Set<N> {
        val readingFrames = mutableSetOf<List<M>>()
        val reverseComplement = reverseComplement()
        readingFrames.add(this.monomers)
        readingFrames.add(this.monomers.drop(1))
        readingFrames.add(this.monomers.drop(2))
        readingFrames.add(reverseComplement.monomers)
        readingFrames.add(reverseComplement.monomers.drop(1))
        readingFrames.add(reverseComplement.monomers.drop(2))
        readingFrames.remove(emptyList())
        return readingFrames.map(::create).toSet()
    }

    fun transitionTransversionRate(other: N): Double {
        var transitions = 0.0
        var transversions = 0.0
        for (nucleotide in this.monomers.zip(other.monomers)) {
            if (nucleotide.first != nucleotide.second) {
                if (nucleotide.first in purines && nucleotide.second in purines) {
                    transitions++
                } else if (nucleotide.first in pyrimidines && nucleotide.second in pyrimidines) {
                    transitions++
                } else {
                    transversions++
                }
            }
        }
        return transitions / transversions
    }

    fun palindromes(min: Int, max: Int): Map<N, List<Int>> {
        val palindromes = mutableMapOf<N, MutableList<Int>>()
        for (size in min..max) {
            val subsequences = monomers.windowed(size)
            for (subsequence in subsequences.map(::create).withIndex()) {
                if (subsequence.value.isPalindrome()) {
                    val indexes = palindromes[subsequence.value]
                    if (indexes == null) {
                        palindromes[subsequence.value] = mutableListOf(subsequence.index)
                    } else {
                        indexes.add(subsequence.index)
                    }
                }
            }
        }
        return palindromes
    }

}