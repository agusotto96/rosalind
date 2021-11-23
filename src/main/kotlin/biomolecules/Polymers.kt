package biomolecules

typealias DNA = List<DNANucleotide>

typealias RNA = List<RNANucleotide>

typealias Protein = List<AminoAcid>

fun DNA(symbols: String): DNA = symbols.mapNotNull(DNANucleotidesBySymbol::get)

fun RNA(symbols: String): RNA = symbols.mapNotNull(RNANucleotidesBySymbol::get)

fun Protein(symbols: String): Protein = symbols.mapNotNull(aminoAcidsBySymbol::get)

private val DNANucleotidesBySymbol = DNANucleotide.values().associateBy(Monomer<*>::symbol)

private val RNANucleotidesBySymbol = RNANucleotide.values().associateBy(Monomer<*>::symbol)

private val aminoAcidsBySymbol = AminoAcid.values().associateBy(Monomer<*>::symbol)