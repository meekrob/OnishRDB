#!/usr/bin/env python3
# conda install mysql-connector-python
import sys
import mysql.connector as mariadb

connection = mariadb.connect(
    host="129.82.125.11",
    port="3307",
    user="worm",
    password="",
    database="information_schema"
    )



cursor = connection.cursor()
tablename = "PromoterPeakOverlap"
if len(sys.argv) > 1:
    tablename = sys.argv[1]

cursor.execute("select DATA_TYPE, COLUMN_NAME, COLUMN_KEY, COLUMN_COMMENT from information_schema.columns where table_name = '%s'; " % tablename)

print('```mermaid')
print()
print('erDiagram')
print()

print(tablename, "{")
for data_type, column_name, column_key, column_comment in cursor:
    if column_key == 'PRI': 
        column_key = 'PK'
    elif column_key == 'MUL':
        column_key = ''

    column_comment = column_comment.replace('"','')

    print("%s %s %s \'%s\'" % (data_type, column_name, column_key, column_comment))

print("}")
