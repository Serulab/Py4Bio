import sqlite3
db = sqlite3.connect('../../samples/PythonU.db')
cursor = db.cursor()
cursor.execute('Select * from Students')
for row in cursor:
    print(row)
