# dotDNA
 SnapGene .dna file parser.


# Examples
Parse .dna file and output primers into a csv:

    from dot_dna import parser
    from csv import DictWriter
    
    column = ["Name", "Position", "Direction", "Temperature", "Sequence"]
    file_path = "tests/test.dna"
    
    snapgene = parser.SnapGene(file_path)
    snapgene.parse()
    with open(output_path, "wt", newline="") as output:
        w = DictWriter(output, fieldnames=column)
        w.writeheader()
        for p in snapgene.primers:
            w.writerow(p.to_dict())

