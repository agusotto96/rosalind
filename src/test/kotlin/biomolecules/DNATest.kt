package biomolecules

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test

internal class DNATest {

    @Test
    fun transcript() {
        val expected = "GAUGGAACUUGACUACGUAAAUU"
            .let(::RNA)
        val actual = "GATGGAACTTGACTACGTAAATT"
            .let(::DNA)
            .let(DNA::transcript)
        assertEquals(expected, actual)
    }

    @Test
    fun gcContent() {
        val expected = 60.91954
        val actual = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
            .let(::DNA)
            .let(DNA::gcContent)
        assertEquals(expected, actual, 0.00001)
    }

    @Test
    fun sharedSubsequence() {
        val expected = "CGTA"
            .let(::DNA)
        val actual = listOf("ACGTACGT", "AACCGTATA")
            .map(::DNA)
            .let(::sharedSubsequence)
        assertEquals(expected, actual)
    }

    @Test
    fun reverseComplement() {
        val expected = "ACCGGGTTTT"
            .let(::DNA)
        val actual = "AAAACCCGGT"
            .let(::DNA)
            .let(DNA::reverseComplement)
        assertEquals(expected, actual)
    }

    @Test
    fun profile() {
        val expected = mapOf(
            Adenine to listOf(5, 1, 0, 0, 5, 5, 0, 0),
            Cytosine to listOf(0, 0, 1, 4, 2, 0, 6, 1),
            Guanine to listOf(1, 1, 6, 3, 0, 1, 0, 0),
            Thymine to listOf(1, 5, 0, 0, 0, 1, 1, 6)
        )
        val actual = listOf("ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC", "TTGGAACT", "ATGCCATT", "ATGGCACT")
            .map(::DNA)
            .let(::profile)
        assertEquals(expected, actual)
    }

    @Test
    fun consensus() {
        val expected = DNA("ATGCAACT")
        val actual = listOf("ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC", "TTGGAACT", "ATGCCATT", "ATGGCACT")
            .map(::DNA)
            .let(::consensus)
        assertEquals(expected, actual)
    }

    @Test
    fun transitionTransversionRate() {
        val a = "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT"
            .let(::DNA)
        val b = "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"
            .let(::DNA)
        val expected = 1.214
        val actual = a.transitionTransversionRate(b)
        assertEquals(expected, actual, 0.001)
    }

}