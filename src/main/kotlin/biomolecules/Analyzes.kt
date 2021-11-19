package biomolecules

fun String(polymer: Polymer) = polymer.map(::symbol).joinToString("")

fun countMonomers(polymer: Polymer) = polymer.groupingBy { it }.eachCount()

fun subPolymer(polymer: Polymer, from: Int, to: Int) = polymer.subList(from, to + 1)

fun countMismatches(a: Polymer, b: Polymer) = a.zip(b).count { it.first != it.second }

fun locate(polymer: Polymer, subPolymer: Polymer) = polymer
    .windowed(subPolymer.count())
    .mapIndexedNotNull { index, window -> if (window == subPolymer) index else null }

fun contains(polymer: Polymer, subPolymer: Polymer) = polymer.windowed(subPolymer.count()).contains(subPolymer)

fun allContain(polymers: List<Polymer>, subPolymer: Polymer) = polymers.all { contains(it, subPolymer) }

fun shortestPolymer(polymers: List<Polymer>) = polymers.minByOrNull(Polymer::count)

fun longestPolymer(polymers: List<Polymer>) = polymers.maxByOrNull(Polymer::count)

fun sharedSubPolymer(polymers: List<Polymer>): Polymer? {
    val shortestPolymer = shortestPolymer(polymers) ?: return null
    for (count in shortestPolymer.count() downTo 1) {
        for (subPolymer in shortestPolymer.windowed(count)) {
            if (allContain(polymers, subPolymer)) return subPolymer
        }
    }
    return null
}

fun profile(polymers: List<Polymer>): Map<Monomer, List<Int>> {
    val profile = mutableMapOf<Monomer, MutableList<Int>>()
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

fun consensus(polymers: List<Polymer>): Polymer? {
    val longestPolymer = longestPolymer(polymers) ?: return null
    val consensus = mutableListOf<Monomer>()
    val profile = profile(polymers)
    for (i in longestPolymer.indices) {
        var maxCount = 0
        var mostCommonSymbol: Monomer? = null
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

fun reverseComplement(nucleicAcid: NucleicAcid) = nucleicAcid.map(::complement).asReversed()

fun isPalindrome(nucleicAcid: NucleicAcid) = nucleicAcid == reverseComplement(nucleicAcid)

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

fun gcContent(nucleicAcid: NucleicAcid) = nucleicAcid.count(::isGc) / (nucleicAcid.count() / 100.0)

fun palindromes(nucleicAcid: NucleicAcid, min: Int, max: Int): Map<NucleicAcid, List<Int>> {
    val palindromes = mutableMapOf<NucleicAcid, MutableList<Int>>()
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

fun transcript(dna: DNA): RNA = dna.map(::transcript)

fun reverseTranscript(rna: RNA): DNA = rna.map(::reverseTranscript)

fun triplets(rna: RNA) = rna.chunked(3).map { Triple(it[0], it[1], it[2]) }

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

fun mass(protein: Protein): Double = protein.sumOf(::mass)