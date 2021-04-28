from unittest import TestCase

import dot_dna.primer
from dot_dna import parser
import os

test_file = os.path.join("tests", "test.dna")


class TestSnapGene(TestCase):
    def setUp(self) -> None:
        self.snap = parser.SnapGene(test_file)
        self.snap.parse()

    def test_parse_primers(self):
        print(self.snap.primers)
        if not self.snap.primers:
            self.fail()

    def test_parse_features(self):
        if not self.snap.features:
            self.fail()

    def test_parse_notes(self):
        if not self.snap.notes_content:
            self.fail()

    def test_parse_seq_properties(self):
        if not self.snap.seq_properties:
            self.fail()

    def test_get_translated(self):
        test1 = self.snap.get_translated(2)
        if not test1[0] == "M" or not test1[1] == 1:
            self.fail()

        test2 = self.snap.get_translated(2, 1)
        print(test2)


class TestPrimer(TestCase):
    def test_from_string(self):
        p = dot_dna.primer.Primer()
        p.from_string("TCTAGTTGTTCCAGAGATATTCCATACC")

