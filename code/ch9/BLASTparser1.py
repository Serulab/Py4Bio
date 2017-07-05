from Bio.Blast import NCBIXML
with open('../../samples/sampleXblast.xml') as xmlfh:
    for record in NCBIXML.parse(xmlfh):
        for align in record.alignments:
            print(align.title)
