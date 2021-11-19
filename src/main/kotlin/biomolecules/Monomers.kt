package biomolecules

sealed interface Monomer

sealed interface Nucleotide : Monomer

enum class DNANucleotide : Nucleotide {
    Adenine,
    Cytosine,
    Guanine,
    Thymine
}

enum class RNANucleotide : Nucleotide {
    Adenine,
    Cytosine,
    Guanine,
    Uracil
}

enum class AminoAcid : Monomer {
    Alanine,
    Cysteine,
    AsparticAcid,
    GlutamicAcid,
    Phenylalanine,
    Glycine,
    Histidine,
    Isoleucine,
    Lysine,
    Leucine,
    Methionine,
    Asparagine,
    Proline,
    Glutamine,
    Arginine,
    Serine,
    Threonine,
    Valine,
    Tryptophan,
    Tyrosine
}

fun symbol(monomer: Monomer): Char = when (monomer) {
    DNANucleotide.Adenine -> 'A'
    DNANucleotide.Cytosine -> 'C'
    DNANucleotide.Guanine -> 'G'
    DNANucleotide.Thymine -> 'T'
    RNANucleotide.Adenine -> 'A'
    RNANucleotide.Cytosine -> 'C'
    RNANucleotide.Guanine -> 'G'
    RNANucleotide.Uracil -> 'U'
    AminoAcid.Alanine -> 'A'
    AminoAcid.Cysteine -> 'C'
    AminoAcid.AsparticAcid -> 'D'
    AminoAcid.GlutamicAcid -> 'E'
    AminoAcid.Phenylalanine -> 'F'
    AminoAcid.Glycine -> 'G'
    AminoAcid.Histidine -> 'H'
    AminoAcid.Isoleucine -> 'I'
    AminoAcid.Lysine -> 'K'
    AminoAcid.Leucine -> 'L'
    AminoAcid.Methionine -> 'M'
    AminoAcid.Asparagine -> 'N'
    AminoAcid.Proline -> 'P'
    AminoAcid.Glutamine -> 'Q'
    AminoAcid.Arginine -> 'R'
    AminoAcid.Serine -> 'S'
    AminoAcid.Threonine -> 'T'
    AminoAcid.Valine -> 'V'
    AminoAcid.Tryptophan -> 'W'
    AminoAcid.Tyrosine -> 'Y'
}

fun complement(nucleotide: Nucleotide): Nucleotide = when (nucleotide) {
    DNANucleotide.Adenine -> DNANucleotide.Thymine
    DNANucleotide.Cytosine -> DNANucleotide.Guanine
    DNANucleotide.Guanine -> DNANucleotide.Cytosine
    DNANucleotide.Thymine -> DNANucleotide.Adenine
    RNANucleotide.Adenine -> RNANucleotide.Uracil
    RNANucleotide.Cytosine -> RNANucleotide.Guanine
    RNANucleotide.Guanine -> RNANucleotide.Cytosine
    RNANucleotide.Uracil -> RNANucleotide.Adenine
}

fun isGc(nucleotide: Nucleotide): Boolean = when (nucleotide) {
    DNANucleotide.Guanine,
    DNANucleotide.Cytosine,
    RNANucleotide.Guanine,
    RNANucleotide.Cytosine -> true
    else -> false
}

fun isPurine(nucleotide: Nucleotide): Boolean = when (nucleotide) {
    DNANucleotide.Adenine,
    DNANucleotide.Guanine,
    RNANucleotide.Adenine,
    RNANucleotide.Guanine -> true
    else -> false
}

fun isPyrimidine(nucleotide: Nucleotide): Boolean = !isPurine(nucleotide)

fun mass(aminoAcid: AminoAcid): Double = when (aminoAcid) {
    AminoAcid.Alanine -> 071.03711
    AminoAcid.Cysteine -> 103.00919
    AminoAcid.AsparticAcid -> 115.02694
    AminoAcid.GlutamicAcid -> 129.04259
    AminoAcid.Phenylalanine -> 147.06841
    AminoAcid.Glycine -> 057.02146
    AminoAcid.Histidine -> 137.05891
    AminoAcid.Isoleucine -> 113.08406
    AminoAcid.Lysine -> 128.09496
    AminoAcid.Leucine -> 113.08406
    AminoAcid.Methionine -> 131.04049
    AminoAcid.Asparagine -> 114.04293
    AminoAcid.Proline -> 097.05276
    AminoAcid.Glutamine -> 128.05858
    AminoAcid.Arginine -> 156.10111
    AminoAcid.Serine -> 087.03203
    AminoAcid.Threonine -> 101.04768
    AminoAcid.Valine -> 099.06841
    AminoAcid.Tryptophan -> 186.07931
    AminoAcid.Tyrosine -> 163.06333
}

fun transcript(nucleotide: DNANucleotide): RNANucleotide = when (nucleotide) {
    DNANucleotide.Adenine -> RNANucleotide.Adenine
    DNANucleotide.Cytosine -> RNANucleotide.Cytosine
    DNANucleotide.Guanine -> RNANucleotide.Guanine
    DNANucleotide.Thymine -> RNANucleotide.Uracil
}

fun reverseTranscript(nucleotide: RNANucleotide): DNANucleotide = when (nucleotide) {
    RNANucleotide.Adenine -> DNANucleotide.Adenine
    RNANucleotide.Cytosine -> DNANucleotide.Cytosine
    RNANucleotide.Guanine -> DNANucleotide.Guanine
    RNANucleotide.Uracil -> DNANucleotide.Thymine
}

fun aminoAcid(triplet: Triple<RNANucleotide, RNANucleotide, RNANucleotide>): AminoAcid? = when (triplet) {
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Uracil) -> AminoAcid.Phenylalanine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Uracil) -> AminoAcid.Leucine
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Uracil) -> AminoAcid.Isoleucine
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Uracil) -> AminoAcid.Valine
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Cytosine) -> AminoAcid.Phenylalanine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Cytosine) -> AminoAcid.Leucine
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Cytosine) -> AminoAcid.Isoleucine
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Cytosine) -> AminoAcid.Valine
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Adenine) -> AminoAcid.Leucine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Adenine) -> AminoAcid.Leucine
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Adenine) -> AminoAcid.Isoleucine
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Adenine) -> AminoAcid.Valine
    Triple(RNANucleotide.Uracil, RNANucleotide.Uracil, RNANucleotide.Guanine) -> AminoAcid.Leucine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Uracil, RNANucleotide.Guanine) -> AminoAcid.Leucine
    Triple(RNANucleotide.Guanine, RNANucleotide.Uracil, RNANucleotide.Guanine) -> AminoAcid.Valine
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Uracil) -> AminoAcid.Serine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Uracil) -> AminoAcid.Proline
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Uracil) -> AminoAcid.Threonine
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Uracil) -> AminoAcid.Alanine
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Cytosine) -> AminoAcid.Serine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Cytosine) -> AminoAcid.Proline
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Cytosine) -> AminoAcid.Threonine
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Cytosine) -> AminoAcid.Alanine
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Adenine) -> AminoAcid.Serine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Adenine) -> AminoAcid.Proline
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Adenine) -> AminoAcid.Threonine
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Adenine) -> AminoAcid.Alanine
    Triple(RNANucleotide.Uracil, RNANucleotide.Cytosine, RNANucleotide.Guanine) -> AminoAcid.Serine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Cytosine, RNANucleotide.Guanine) -> AminoAcid.Proline
    Triple(RNANucleotide.Adenine, RNANucleotide.Cytosine, RNANucleotide.Guanine) -> AminoAcid.Threonine
    Triple(RNANucleotide.Guanine, RNANucleotide.Cytosine, RNANucleotide.Guanine) -> AminoAcid.Alanine
    Triple(RNANucleotide.Uracil, RNANucleotide.Adenine, RNANucleotide.Uracil) -> AminoAcid.Tyrosine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Uracil) -> AminoAcid.Histidine
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Uracil) -> AminoAcid.Asparagine
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Uracil) -> AminoAcid.AsparticAcid
    Triple(RNANucleotide.Uracil, RNANucleotide.Adenine, RNANucleotide.Cytosine) -> AminoAcid.Tyrosine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Cytosine) -> AminoAcid.Histidine
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Cytosine) -> AminoAcid.Asparagine
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Cytosine) -> AminoAcid.AsparticAcid
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Adenine) -> AminoAcid.Glutamine
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Adenine) -> AminoAcid.Lysine
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Adenine) -> AminoAcid.GlutamicAcid
    Triple(RNANucleotide.Cytosine, RNANucleotide.Adenine, RNANucleotide.Guanine) -> AminoAcid.Glutamine
    Triple(RNANucleotide.Adenine, RNANucleotide.Adenine, RNANucleotide.Guanine) -> AminoAcid.Lysine
    Triple(RNANucleotide.Guanine, RNANucleotide.Adenine, RNANucleotide.Guanine) -> AminoAcid.GlutamicAcid
    Triple(RNANucleotide.Uracil, RNANucleotide.Guanine, RNANucleotide.Uracil) -> AminoAcid.Cysteine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Uracil) -> AminoAcid.Arginine
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Uracil) -> AminoAcid.Serine
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Uracil) -> AminoAcid.Glycine
    Triple(RNANucleotide.Uracil, RNANucleotide.Guanine, RNANucleotide.Cytosine) -> AminoAcid.Cysteine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Cytosine) -> AminoAcid.Arginine
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Cytosine) -> AminoAcid.Serine
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Cytosine) -> AminoAcid.Glycine
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Adenine) -> AminoAcid.Arginine
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Adenine) -> AminoAcid.Arginine
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Adenine) -> AminoAcid.Glycine
    Triple(RNANucleotide.Uracil, RNANucleotide.Guanine, RNANucleotide.Guanine) -> AminoAcid.Tryptophan
    Triple(RNANucleotide.Cytosine, RNANucleotide.Guanine, RNANucleotide.Guanine) -> AminoAcid.Arginine
    Triple(RNANucleotide.Adenine, RNANucleotide.Guanine, RNANucleotide.Guanine) -> AminoAcid.Arginine
    Triple(RNANucleotide.Guanine, RNANucleotide.Guanine, RNANucleotide.Guanine) -> AminoAcid.Glycine
    Triple(RNANucleotide.Adenine, RNANucleotide.Uracil, RNANucleotide.Guanine) -> AminoAcid.Methionine
    else -> null
}

val DNANucleotidesBySymbol = DNANucleotide.values().associateBy(::symbol)

val RNANucleotidesBySymbol = RNANucleotide.values().associateBy(::symbol)

val aminoAcidsBySymbol = AminoAcid.values().associateBy(::symbol)