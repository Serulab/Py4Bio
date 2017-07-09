import pymysql
db = pymysql.connect(host='localhost',
        user='root', passwd='secret', db='PythonU')
cursor = db.cursor()
cursor.execute('SELECT * FROM Students')
for row in cursor:
    print(row)
