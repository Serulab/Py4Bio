from Bio import Entrez
my_em = 'user@example.com'
db = "gene"
term = 'cobalamin synthase homo sapiens'
h_search = Entrez.esearch(db=db, email=my_em, term=term)
record = Entrez.read(h_search)
res_ids = record["IdList"]
for r_id in res_ids:
   h_summ = Entrez.esummary(db=db, id=r_id, email=my_em)
   s = Entrez.read(h_summ)
   print(r_id)
   name = s['DocumentSummarySet']['DocumentSummary'][0]['Name']
   print(name)
   su = s['DocumentSummarySet']['DocumentSummary'][0]['Summary']
   print(su)
   print('==============================================')
