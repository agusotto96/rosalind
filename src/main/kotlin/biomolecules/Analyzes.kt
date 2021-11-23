package biomolecules

fun <M : Monomer<M>> String(polymer: List<M>): String {
    return polymer.map(Monomer<M>::symbol).joinToString("")
}

fun <M : Monomer<M>> countMonomers(polymer: List<M>): Map<M, Int> {
    return polymer.groupingBy { it }.eachCount()
}

fun <M : Monomer<M>> subPolymer(polymer: List<M>, from: Int, to: Int): List<M> {
    return polymer.subList(from, to + 1)
}

fun <M : Monomer<M>> countMismatches(a: List<M>, b: List<M>): Int {
    return a.zip(b).count { it.first != it.second }
}

fun <M : Monomer<M>> locate(polymer: List<M>, subPolymer: List<M>): List<Int> {
    return polymer
        .windowed(subPolymer.count())
        .mapIndexedNotNull { index, window -> if (window == subPolymer) index else null }
}

fun <M : Monomer<M>> contains(polymer: List<M>, subPolymer: List<M>): Boolean {
    return polymer.windowed(subPolymer.count()).contains(subPolymer)
}

fun <M : Monomer<M>> allContain(polymers: List<List<M>>, subPolymer: List<M>): Boolean {
    return polymers.all { contains(it, subPolymer) }
}

fun <M : Monomer<M>> shortestPolymer(polymers: List<List<M>>): List<M>? {
    return polymers.minByOrNull(List<M>::count)
}

fun <M : Monomer<M>> longestPolymer(polymers: List<List<M>>): List<M>? {
    return polymers.maxByOrNull(List<M>::count)
}

fun <M : Monomer<M>> sharedSubPolymer(polymers: List<List<M>>): List<M>? {
    val shortestPolymer = shortestPolymer(polymers) ?: return null
    for (count in shortestPolymer.count() downTo 1) {
        for (subPolymer in shortestPolymer.windowed(count)) {
            if (allContain(polymers, subPolymer)) return subPolymer
        }
    }
    return null
}

fun <M : Monomer<M>> profile(polymers: List<List<M>>): Map<M, List<Int>> {
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

fun <M : Monomer<M>> consensus(polymers: List<List<M>>): List<M>? {
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

fun <N : Nucleotide<N>> reverseComplement(nucleicAcid: List<N>): List<N> {
    return nucleicAcid.map(Nucleotide<N>::complement).asReversed()
}

fun <N : Nucleotide<N>> isPalindrome(nucleicAcid: List<N>): Boolean {
    return nucleicAcid == reverseComplement(nucleicAcid)
}

fun <N : Nucleotide<N>> readingFrames(nucleicAcid: List<N>): List<List<N>> {
    val readingFrames = mutableListOf<List<N>>()
    val reverseComplement = reverseComplement(nucleicAcid)
    readingFrames.add(nucleicAcid)
    readingFrames.add(nucleicAcid.drop(1))
    readingFrames.add(nucleicAcid.drop(2))
    readingFrames.add(reverseComplement)
    readingFrames.add(reverseComplement.drop(1))
    readingFrames.add(reverseComplement.drop(2))
    return readingFrames
}

fun <N : Nucleotide<N>> gcContent(nucleicAcid: List<N>): Double {
    return nucleicAcid.count(Nucleotide<N>::isGc) / (nucleicAcid.count() / 100.0)
}

fun <N : Nucleotide<N>> palindromes(nucleicAcid: List<N>, min: Int, max: Int): Map<List<N>, List<Int>> {
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

fun <N : Nucleotide<N>> transitionTransversionRate(a: List<N>, b: List<N>): Double {
    var transitions = 0.0
    var transversions = 0.0
    for (nucleotide in a.zip(b)) {
        if (nucleotide.first != nucleotide.second) {
            if (nucleotide.first.isPurine() && nucleotide.second.isPurine()) {
                transitions++
            } else if (nucleotide.first.isPyrimidine() && nucleotide.second.isPyrimidine()) {
                transitions++
            } else {
                transversions++
            }
        }
    }
    return transitions / transversions
}

fun transcript(dna: DNA): RNA {
    return dna.map(DNANucleotide::transcript)
}

fun reverseTranscript(rna: RNA): DNA {
    return rna.map(RNANucleotide::reverseTranscript)
}

fun triplets(rna: RNA): List<Triple<RNANucleotide, RNANucleotide, RNANucleotide>> {
    return rna.chunked(3).map { Triple(it[0], it[1], it[2]) }
}

fun translate(
    rna: RNA,
    aminoAcidsByTriplet: Map<Triple<RNANucleotide, RNANucleotide, RNANucleotide>, AminoAcid?> = defaultAminoAcidsByTriplet,
    startingAminoAcids: Set<AminoAcid> = defaultStartingAminoAcids
): Set<Protein> {
    val candidates = mutableSetOf<MutableList<AminoAcid>>()
    val proteins = mutableSetOf<Protein>()
    for (triplet in triplets(rna)) {
        val aminoAcid = aminoAcidsByTriplet[triplet]
        if (aminoAcid != null) {
            if (startingAminoAcids.contains(aminoAcid)) {
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
    return protein.sumOf(AminoAcid::mass)
}

private val defaultStartingAminoAcids = setOf(AminoAcid.Methionine)

private val defaultAminoAcidsByTriplet = mapOf(
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Uracil) to AminoAcid.Phenylalanine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Uracil) to AminoAcid.Leucine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Uracil) to AminoAcid.Isoleucine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Uracil) to AminoAcid.Valine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Cytosine) to AminoAcid.Phenylalanine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Cytosine) to AminoAcid.Leucine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Cytosine) to AminoAcid.Isoleucine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Cytosine) to AminoAcid.Valine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Adenine) to AminoAcid.Leucine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Adenine) to AminoAcid.Leucine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Adenine) to AminoAcid.Isoleucine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Adenine) to AminoAcid.Valine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Guanine) to AminoAcid.Leucine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Guanine) to AminoAcid.Leucine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Guanine) to AminoAcid.Valine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Uracil) to AminoAcid.Serine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Uracil) to AminoAcid.Proline,
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Uracil) to AminoAcid.Threonine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Uracil) to AminoAcid.Alanine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Cytosine) to AminoAcid.Serine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Cytosine) to AminoAcid.Proline,
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Cytosine) to AminoAcid.Threonine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Cytosine) to AminoAcid.Alanine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Adenine) to AminoAcid.Serine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Adenine) to AminoAcid.Proline,
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Adenine) to AminoAcid.Threonine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Adenine) to AminoAcid.Alanine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Guanine) to AminoAcid.Serine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Guanine) to AminoAcid.Proline,
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Guanine) to AminoAcid.Threonine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Guanine) to AminoAcid.Alanine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Adenine, RNANucleotide.Uracil) to AminoAcid.Tyrosine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Uracil) to AminoAcid.Histidine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Uracil) to AminoAcid.Asparagine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Uracil) to AminoAcid.AsparticAcid,
    Triple(RNANucleotide.Uracil, RNANucleotide.Adenine, RNANucleotide.Cytosine) to AminoAcid.Tyrosine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Cytosine) to AminoAcid.Histidine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Cytosine) to AminoAcid.Asparagine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Cytosine) to AminoAcid.AsparticAcid,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Adenine) to AminoAcid.Glutamine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Adenine) to AminoAcid.Lysine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Adenine) to AminoAcid.GlutamicAcid,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Guanine) to AminoAcid.Glutamine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Guanine) to AminoAcid.Lysine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Guanine) to AminoAcid.GlutamicAcid,
    Triple(RNANucleotide.Uracil, RNANucleotide.Guanine, RNANucleotide.Uracil) to AminoAcid.Cysteine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Uracil) to AminoAcid.Arginine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Uracil) to AminoAcid.Serine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Uracil) to AminoAcid.Glycine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Guanine, RNANucleotide.Cytosine) to AminoAcid.Cysteine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Cytosine) to AminoAcid.Arginine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Cytosine) to AminoAcid.Serine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Cytosine) to AminoAcid.Glycine,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Adenine) to AminoAcid.Arginine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Adenine) to AminoAcid.Arginine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Adenine) to AminoAcid.Glycine,
    Triple(RNANucleotide.Uracil, RNANucleotide.Guanine, RNANucleotide.Guanine) to AminoAcid.Tryptophan,
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Guanine) to AminoAcid.Arginine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Guanine) to AminoAcid.Arginine,
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Guanine) to AminoAcid.Glycine,
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Guanine) to AminoAcid.Methionine
)
