import pprint
from Bio.Sequencing import Phd
fn = '../../samples/phd1'
with open(fn) as fh:
    rp = Phd.read(fh)
# All the comments are in a dictionary
pprint.pprint(rp.comments)
# Sequence information
print('Sequence: %s' % rp.seq)
# Quality information for each base
print('Quality: %s' % rp.sites)
