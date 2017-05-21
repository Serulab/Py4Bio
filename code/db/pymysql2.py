import pymysql
db = pymysql.connect(host="", user="", passwd="", db="")
cursor = db.cursor()
cursor.execute('SELECT * FROM Students')
for rec in cursor:
    print(rec)
