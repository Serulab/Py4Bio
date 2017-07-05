from Bio import Entrez
my_em = 'user@example.com'
db = "pubmed"
# Search de Entrez website using esearch from eUtils
# esearch returns a handle (called h_search)
h_search = Entrez.esearch(db=db, email=my_em,
                         term='python and bioinformatics')
# Parse the result with Entrez.read()
record = Entrez.read(h_search)
# Get the list of Ids returned by previous search
res_ids = record["IdList"]
# For each id in the list
for r_id in res_ids:
   # Get summary information for each id
   h_summ = Entrez.esummary(db=db, id=r_id, email=my_em)
   # Parse the result with Entrez.read()
   summ = Entrez.read(h_summ)
   print(summ[0]['Title'])
   print(summ[0]['DOI'])
   print('==============================================')
