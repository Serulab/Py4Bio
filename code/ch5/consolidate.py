from path import Path
d = Path('../../samples/bioinfo/seqs/')
with open('outfile.fasta', 'w') as f_out:
    for file_name in d.walk('*.fasta'):
        with open(file_name) as f_in:
            data = f_in.read()
            f_out.write(data)
