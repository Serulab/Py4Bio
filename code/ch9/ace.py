from Bio.Sequencing import Ace

fn = '../../samples/contig1.ace'
acefilerecord = Ace.read(open(fn))

# For each contig:
for ctg in acefilerecord.contigs:
    print('==========================================')
    print('Contig name: %s'%ctg.name)
    print('Bases: %s'%ctg.nbases)
    print('Reads: %s'%ctg.nreads)
    print('Segments: %s'%ctg.nsegments)
    print('Sequence: %s'%ctg.sequence)
    print('Quality: %s'%ctg.quality)
    # For each read in contig:
    for read in ctg.reads:
        print('Read name: %s'%read.rd.name)
        print('Align start: %s'%read.qa.align_clipping_start)
        print('Align end: %s'%read.qa.align_clipping_end)
        print('Qual start: %s'%read.qa.qual_clipping_start)
        print('Qual end: %s'%read.qa.qual_clipping_end)
        print('Read sequence: %s'%read.rd.sequence)
        print('==========================================')