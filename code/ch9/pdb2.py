import gzip
import io
from Bio.PDB.PDBParser import PDBParser

def disorder(structure):
   for chain in structure[0].get_list():
       for residue in chain.get_list():
           for atom in residue.get_list():
               if atom.is_disordered():
                   print(residue, atom)
   return None

pdbfn = '../../samples/pdb1apk.ent.gz'
handle = gzip.GzipFile(pdbfn)
handle = io.StringIO(handle.read().decode('utf-8'))
parser = PDBParser()
structure = parser.get_structure('test', handle)
disorder(structure)
