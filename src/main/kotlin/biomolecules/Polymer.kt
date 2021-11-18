package biomolecules

sealed class Polymer<P : Polymer<P, M>, M : Monomer>(protected val monomers: List<M>) {

    init {
        if (monomers.isEmpty()) throw IllegalArgumentException("Polymer cannot be empty.")
    }

    fun subsequence(from: Int, to: Int) = create(monomers.subList(from, to + 1))

    fun countMismatches(other: P) = this.monomers.zip(other.monomers).count { it.first != it.second }

    fun monomersByCount() = monomers.groupingBy { it }.eachCount()

    fun count() = monomers.count()

    fun windowed(size: Int) = monomers.windowed(size).map(::create)

    fun locate(subsequence: P) = windowed(subsequence.count())
        .mapIndexedNotNull { index, window -> if (window == subsequence) index else null }

    fun contains(subsequence: P) = windowed(subsequence.count()).contains(subsequence)

    override fun toString() = monomers.map(Monomer::symbol).joinToString("")

    protected abstract fun create(monomers: List<M>): P

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        other as Polymer<*, *>
        if (monomers != other.monomers) return false
        return true
    }

    override fun hashCode() = monomers.hashCode()

    fun withIndex() = monomers.withIndex()

    fun indices() = monomers.indices

}

fun <P : Polymer<P, M>, M : Monomer> sharedSubsequence(sequences: List<P>): P? {
    val shortestSequence = sequences.minByOrNull { it.count() } ?: return null
    for (i in shortestSequence.count() downTo 1) {
        val subsequences = shortestSequence.windowed(i)
        for (subsequence in subsequences) {
            if (sequences.all { it.contains(subsequence) }) {
                return subsequence
            }
        }
    }
    return null
}

fun <P : Polymer<P, M>, M : Monomer> profile(sequences: List<P>): Map<M, List<Int>> {
    val profile = mutableMapOf<M, MutableList<Int>>()
    val longestSequence = sequences.maxByOrNull { it.count() } ?: return profile
    for (sequence in sequences) {
        for (symbol in sequence.withIndex()) {
            val count = profile[symbol.value]
            if (count == null) {
                val initialCount = MutableList(longestSequence.count()) { 0 }
                initialCount[symbol.index] = 1
                profile[symbol.value] = initialCount
            } else {
                count[symbol.index]++
            }
        }
    }
    return profile
}

internal fun <P : Polymer<P, M>, M : Monomer> consensus(sequences: List<P>, create: (List<M>) -> P): P? {
    val longestSequence = sequences.maxByOrNull { it.count() } ?: return null
    val consensus = mutableListOf<M>()
    val profile = profile(sequences)
    for (i in longestSequence.indices()) {
        var maxCount = 0
        var mostCommonSymbol: M? = null
        for (row in profile) {
            if (row.value[i] > maxCount) {
                maxCount = row.value[i]
                mostCommonSymbol = row.key
            }
        }
        if (mostCommonSymbol == null) return null
        consensus.add(mostCommonSymbol)
    }
    return create(consensus)
}
