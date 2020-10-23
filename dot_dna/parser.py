import struct
import xml.etree.ElementTree as ET

from dot_dna.primer import Primer


class SnapGene:
    def __init__(self, file_path):
        self.file_path = file_path
        self.primers = []
        self.features = {}
        self.notes_content = {}
        self.seq_properties = {}
        self.meta = {}
        self.translation = ""

    def parse(self):
        with open(self.file_path, "rb") as infile:
            if infile.read(1) == b"\t":
                document_length = struct.unpack(">I", infile.read(4))
                if document_length[0] == 14:
                    title = struct.unpack("8s", infile.read(8))
                    if title[0] == b"SnapGene":
                        is_dna = False
                        if struct.unpack(">H", infile.read(2))[0] == 1:
                            is_dna = True
                        self.meta = dict(
                            is_dna=is_dna,
                            export_version=struct.unpack(">H", infile.read(2))[0],
                            import_version=struct.unpack(">H", infile.read(2))[0]
                        )
                        while True:
                            block = infile.read(1)
                            if block == b"":
                                break
                            block_size = struct.unpack(">I", infile.read(4))
                            # print(block, ord(block))
                            if ord(block) == 0:
                                self.parse_seq_properties(block_size, infile)
                            elif ord(block) == 6:
                                root = self.get_xml(block_size, infile)
                                self.parse_notes(root)
                            elif ord(block) == 10:
                                root = self.get_xml(block_size, infile)
                                self.parse_features(root)

                            elif ord(block) == 5:
                                root = self.get_xml(block_size, infile)
                                self.parse_primers(root)
                            else:
                                infile.read(block_size[0])
                                pass
                    else:
                        raise ValueError("File is in an incorrect format.")
                else:
                    raise ValueError("File is in an incorrect format.")
            else:
                raise ValueError("File is in an incorrect format.")

    def parse_primers(self, root):
        for primer in root:
            if primer.tag == "Primer":
                pri = Primer()
                pri.from_element(primer)
                self.primers.append(pri)
            # if primer.tag == "Primer":
            #     print(primer)
            #     pr = dict(primer.attrib)
            #     pr["data"] = {"Sites": []}
            #     for b in primer:
            #         site = dict(b.attrib)
            #         site["components"] = []
            #         for c in b:
            #             site["components"].append(c.attrib)
            #         pr["data"]["Sites"].append(site)
            #    primers.append(pr)

    def parse_features(self, root):
        for f in root:
            feature = dict(f.attrib)
            feature["data"] = {"Q": {}}
            if f.tag not in self.features:
                self.features[f.tag] = []
            for i in f:
                if i.tag == "Q":
                    for v in i:
                        for a in v.attrib:
                            feature["data"]["Q"][i.attrib["name"]] = v.attrib[a]
                            if i.attrib["name"] == "translation":
                                self.translation = v.attrib[a][:]
                else:
                    feature["data"][i.tag] = i.attrib
            self.features[f.tag].append(feature)

    def parse_notes(self, root):
        for i in root:
            self.notes_content[i.tag] = i.text

    def get_xml(self, block_size, infile):
        xml_content = infile.read(block_size[0])
        root = ET.fromstring(struct.unpack("%ss" % len(xml_content), xml_content)[0])
        return root

    def parse_seq_properties(self, block_size, infile):
        p = struct.unpack(">b", infile.read(1))
        # print(p)
        if p[0] & 0x01:
            self.seq_properties["topology"] = "circular"
        else:
            self.seq_properties["topology"] = "linear"
        if p[0] & 0x02 > 0:
            self.seq_properties["stranded"] = "double"
        else:
            self.seq_properties["stranded"] = "single"
        self.seq_properties["a_methylated"] = p[0] & 0x04 > 0
        self.seq_properties["c_methylated"] = p[0] & 0x08 > 0
        self.seq_properties["ki_methylated"] = p[0] & 0x10 > 0
        self.seq_properties["length"] = block_size[0] - 1
        self.seq_properties["seq"] = struct.unpack("%ss" % self.seq_properties["length"],
                                                   infile.read(self.seq_properties["length"]))

    def get_translated(self, position, strand=0):
        if strand == 1:
            position = self.seq_properties["length"] - position
        aa_position = int((position - position % 3)/3)
        return self.translation[aa_position], aa_position+1


