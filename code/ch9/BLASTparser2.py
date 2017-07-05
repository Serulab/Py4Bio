from Bio.Blast import NCBIXML
threshold = 0.0001
with open('../../samples/other.xml') as xmlfh:
    blast_record = next(NCBIXML.parse(xmlfh))
for align in blast_record.alignments:
    if align.hsps[0].expect < threshold:
        print(align.accession)
