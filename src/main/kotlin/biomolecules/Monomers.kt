package biomolecules

internal sealed interface Monomer<Self : Monomer<Self>> {
    fun symbol(): Char
}

internal sealed interface Nucleotide<Self : Nucleotide<Self>> : Monomer<Self> {
    fun complement(): Self
    fun isGc(): Boolean
    fun isPurine(): Boolean
    fun isPyrimidine(): Boolean
}

enum class DNANucleotide : Nucleotide<DNANucleotide> {

    Adenine, Cytosine, Guanine, Thymine;

    override fun complement() = when (this) {
        Adenine -> Thymine
        Cytosine -> Guanine
        Guanine -> Cytosine
        Thymine -> Adenine
    }

    override fun symbol() = when (this) {
        Adenine -> 'A'
        Cytosine -> 'C'
        Guanine -> 'G'
        Thymine -> 'T'
    }

    override fun isGc() = when (this) {
        Guanine,
        Cytosine -> true
        else -> false
    }

    override fun isPurine() = when (this) {
        Adenine,
        Guanine -> true
        else -> false
    }

    override fun isPyrimidine() = when (this) {
        Cytosine,
        Thymine -> true
        else -> false
    }

    fun transcript() = when (this) {
        Adenine -> RNANucleotide.Adenine
        Cytosine -> RNANucleotide.Cytosine
        Guanine -> RNANucleotide.Guanine
        Thymine -> RNANucleotide.Uracil
    }

}

enum class RNANucleotide : Nucleotide<RNANucleotide> {

    Adenine, Cytosine, Guanine, Uracil;

    override fun complement() = when (this) {
        Adenine -> Uracil
        Cytosine -> Guanine
        Guanine -> Cytosine
        Uracil -> Adenine
    }

    override fun symbol() = when (this) {
        Adenine -> 'A'
        Cytosine -> 'C'
        Guanine -> 'G'
        Uracil -> 'U'
    }

    override fun isGc() = when (this) {
        Guanine,
        Cytosine -> true
        else -> false
    }

    override fun isPurine() = when (this) {
        Adenine,
        Guanine -> true
        else -> false
    }

    override fun isPyrimidine() = when (this) {
        Cytosine,
        Uracil -> true
        else -> false
    }

    fun reverseTranscript() = when (this) {
        Adenine -> DNANucleotide.Adenine
        Cytosine -> DNANucleotide.Cytosine
        Guanine -> DNANucleotide.Guanine
        Uracil -> DNANucleotide.Thymine
    }

}

enum class AminoAcid : Monomer<AminoAcid> {

    Alanine, Cysteine, AsparticAcid, GlutamicAcid, Phenylalanine, Glycine, Histidine, Isoleucine, Lysine, Leucine,
    Methionine, Asparagine, Proline, Glutamine, Arginine, Serine, Threonine, Valine, Tryptophan, Tyrosine;

    override fun symbol() = when (this) {
        Alanine -> 'A'
        Cysteine -> 'C'
        AsparticAcid -> 'D'
        GlutamicAcid -> 'E'
        Phenylalanine -> 'F'
        Glycine -> 'G'
        Histidine -> 'H'
        Isoleucine -> 'I'
        Lysine -> 'K'
        Leucine -> 'L'
        Methionine -> 'M'
        Asparagine -> 'N'
        Proline -> 'P'
        Glutamine -> 'Q'
        Arginine -> 'R'
        Serine -> 'S'
        Threonine -> 'T'
        Valine -> 'V'
        Tryptophan -> 'W'
        Tyrosine -> 'Y'
    }

    fun mass() = when (this) {
        Alanine -> 071.03711
        Cysteine -> 103.00919
        AsparticAcid -> 115.02694
        GlutamicAcid -> 129.04259
        Phenylalanine -> 147.06841
        Glycine -> 057.02146
        Histidine -> 137.05891
        Isoleucine -> 113.08406
        Lysine -> 128.09496
        Leucine -> 113.08406
        Methionine -> 131.04049
        Asparagine -> 114.04293
        Proline -> 097.05276
        Glutamine -> 128.05858
        Arginine -> 156.10111
        Serine -> 087.03203
        Threonine -> 101.04768
        Valine -> 099.06841
        Tryptophan -> 186.07931
        Tyrosine -> 163.06333
    }

}
