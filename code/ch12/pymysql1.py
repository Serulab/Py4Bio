import pymysql
db = pymysql.connect(host="", user="", passwd="", db="")
cursor = db.cursor()
recs = cursor.execute('SELECT * FROM Students')
for rec in range(recs):
    print(cursor.fetchone())
