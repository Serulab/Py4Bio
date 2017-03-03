INPUT_FILE = 'fasta22_out.fas'
OUTPUT_FILE = 'fasta33.fas'

with open(INPUT_FILE) as in_fh:
    with open(OUTPUT_FILE, 'w') as out_fh:
        for line in in_fh:
            if line.startswith('>'):
                line = line.replace('\n', '[Rattus norvegicus]\n')
            out_fh.write(line)
