from unittest import TestCase
from dot_dna import parser
import os

test_file = os.path.join("tests", "test.dna")


class TestSnapGene(TestCase):
    def setUp(self) -> None:
        self.snap = parser.SnapGene(test_file)
        self.snap.parse()

    def test_parse_primers(self):
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
