package biomolecules

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test

internal class ProteinTest {

    @Test
    fun mass() {
        val expected = 821.392
        val actual = "SKADYEK"
            .let(::Protein)
            .let(Protein::mass)
        assertEquals(expected, actual, 0.001)
    }

}