from xml.etree.ElementTree import Element


class Primer:
    def __init__(self):
        self.strand = "0"
        self.start = 0
        self.stop = 0
        self.seq = ""
        self.temperature = 0
        self.name = ""
        self.composition = {}
        self.length = 0
    def from_element(self, element: Element):
        if element.tag == "Primer":
            for k in element.attrib:
                if k == "name":
                    self.name = element.attrib[k]
                elif k == "sequence":
                    self.set_seq(element.attrib[k])
            for bindingSite in element:
                if "simplified" not in bindingSite.attrib:
                    location = bindingSite.attrib["location"].split("-")
                    self.start = int(location[0])
                    self.stop = int(location[1])
                self.temperature = int(bindingSite.attrib["meltingTemperature"])
                self.strand = bindingSite.attrib["boundStrand"]

    def from_string(self, seq):
        self.set_seq(seq)
        self.temperature = self.calculate_melting_temp()

    def calculate_melting_temp(self):
        if len(self.seq) > 13:
            return 64.9 + 41*(self.composition["G"]+self.composition["C"]-16.4)/(self.composition["G"]+self.composition["C"]+self.composition["A"]+self.composition["T"])
        else:
            return 2*(self.composition["A"]+self.composition["T"]) + 4*(self.composition["G"]+self.composition["C"])

    def to_dict(self):
        return {"Position": "{}-{}".format(self.start, self.stop), "Name": self.name, "Direction": "Forward" if self.strand == "0" else "Reverse", "Sequence": self.seq, "Temperature": self.temperature}

    def __repr__(self):
        l = "{} {} {}-{} T:{:.2f}".format(self.strand, self.name, self.start, self.stop, self.temperature)
        for i in "ATGC":
            l += " {}:{:.2f}%".format(i, self.composition[i]/self.length*100)
        return l
    def set_seq(self, seq):
        self.seq = seq
        self.composition = {}
        self.length = len(seq)
        for i in "ATGC":
            self.composition[i] = self.seq.count(i)