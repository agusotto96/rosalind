package biomolecules

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test

internal class RNATest {

    @Test
    fun translate() {
        val expected = listOf("MAMAPRTEINSTRING", "MAPRTEINSTRING")
            .map(::Protein)
            .toSet()
        val actual = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
            .let(::RNA)
            .let(RNA::translate)
        assertEquals(expected, actual)
    }

}