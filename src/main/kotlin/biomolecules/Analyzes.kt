package biomolecules

fun String(polymer: Polymer): String {
    return polymer.map(::symbol).joinToString("")
}

fun <M : Monomer> countMonomers(polymer: List<M>): Map<M, Int> {
    return polymer.groupingBy { it }.eachCount()
}

fun <M : Monomer> subPolymer(polymer: List<M>, from: Int, to: Int): List<M> {
    return polymer.subList(from, to + 1)
}

fun countMismatches(a: Polymer, b: Polymer): Int {
    return a.zip(b).count { it.first != it.second }
}

fun locate(polymer: Polymer, subPolymer: Polymer): List<Int> {
    return polymer
        .windowed(subPolymer.count())
        .mapIndexedNotNull { index, window -> if (window == subPolymer) index else null }
}

fun contains(polymer: Polymer, subPolymer: Polymer): Boolean {
    return polymer.windowed(subPolymer.count()).contains(subPolymer)
}

fun allContain(polymers: List<Polymer>, subPolymer: Polymer): Boolean {
    return polymers.all { contains(it, subPolymer) }
}

fun <M : Monomer> shortestPolymer(polymers: List<List<M>>): List<M>? {
    return polymers.minByOrNull(Polymer::count)
}

fun <M : Monomer> longestPolymer(polymers: List<List<M>>): List<M>? {
    return polymers.maxByOrNull(Polymer::count)
}

fun <M : Monomer> sharedSubPolymer(polymers: List<List<M>>): List<M>? {
    val shortestPolymer = shortestPolymer(polymers) ?: return null
    for (count in shortestPolymer.count() downTo 1) {
        for (subPolymer in shortestPolymer.windowed(count)) {
            if (allContain(polymers, subPolymer)) return subPolymer
        }
    }
    return null
}

fun <M : Monomer> profile(polymers: List<List<M>>): Map<M, List<Int>> {
    val profile = mutableMapOf<M, MutableList<Int>>()
    val longestSequence = polymers.maxByOrNull { it.count() } ?: return profile
    for (sequence in polymers) {
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

fun <M : Monomer> consensus(polymers: List<List<M>>): List<M>? {
    val longestPolymer = longestPolymer(polymers) ?: return null
    val consensus = mutableListOf<M>()
    val profile = profile(polymers)
    for (i in longestPolymer.indices) {
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
    return consensus
}

fun reverseComplement(nucleicAcid: NucleicAcid): NucleicAcid {
    return nucleicAcid.map(::complement).asReversed()
}

fun isPalindrome(nucleicAcid: NucleicAcid): Boolean {
    return nucleicAcid == reverseComplement(nucleicAcid)
}

fun readingFrames(nucleicAcid: NucleicAcid): Set<NucleicAcid> {
    val readingFrames = mutableSetOf<NucleicAcid>()
    val reverseComplement = reverseComplement(nucleicAcid)
    readingFrames.add(nucleicAcid)
    readingFrames.add(nucleicAcid.drop(1))
    readingFrames.add(nucleicAcid.drop(2))
    readingFrames.add(reverseComplement)
    readingFrames.add(reverseComplement.drop(1))
    readingFrames.add(reverseComplement.drop(2))
    return readingFrames
}

fun gcContent(nucleicAcid: NucleicAcid): Double {
    return nucleicAcid.count(::isGc) / (nucleicAcid.count() / 100.0)
}

fun <N : Nucleotide> palindromes(nucleicAcid: List<N>, min: Int, max: Int): Map<List<N>, List<Int>> {
    val palindromes = mutableMapOf<List<N>, MutableList<Int>>()
    for (size in min..max) {
        val subsequences = nucleicAcid.windowed(size)
        for (subsequence in subsequences.withIndex()) {
            if (isPalindrome(subsequence.value)) {
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

fun transitionTransversionRate(a: NucleicAcid, b: NucleicAcid): Double {
    var transitions = 0.0
    var transversions = 0.0
    for (nucleotide in a.zip(b)) {
        if (nucleotide.first != nucleotide.second) {
            if (isPurine(nucleotide.first) && isPurine(nucleotide.second)) {
                transitions++
            } else if (isPyrimidine(nucleotide.first) && isPyrimidine(nucleotide.second)) {
                transitions++
            } else {
                transversions++
            }
        }
    }
    return transitions / transversions
}

fun transcript(dna: DNA): RNA {
    return dna.map(::transcript)
}

fun reverseTranscript(rna: RNA): DNA {
    return rna.map(::reverseTranscript)
}

fun triplets(rna: RNA): List<Triple<RNANucleotide, RNANucleotide, RNANucleotide>> {
    return rna.chunked(3).map { Triple(it[0], it[1], it[2]) }
}

fun translate(rna: RNA): Set<Protein> {
    val candidates = mutableSetOf<MutableList<AminoAcid>>()
    val proteins = mutableSetOf<Protein>()
    for (triplet in triplets(rna)) {
        val aminoAcid = aminoAcid(triplet)
        if (aminoAcid != null) {
            if (aminoAcid == AminoAcid.Methionine) {
                candidates.add(mutableListOf())
            }
            candidates.forEach { it.add(aminoAcid) }
        } else {
            proteins.addAll(candidates)
            candidates.clear()
        }
    }
    return proteins
}

fun mass(protein: Protein): Double {
    return protein.sumOf(::mass)
}