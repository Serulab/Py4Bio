import pymysql
db = pymysql.connect(host='localhost',
        user='root', passwd='secret', db='PythonU')
cursor = db.cursor()
recs = cursor.execute('SELECT * FROM Students')
for x in range(recs):
    print(cursor.fetchone())
