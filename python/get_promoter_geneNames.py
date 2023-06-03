#!/usr/bin/env python3
# conda install mysql-connector-python
import sys, pickle
import mysql.connector as mariadb



connection = mariadb.connect(
    host="129.82.125.11",
    port="3307",
    user="worm",
    database="NishimuraLab"
    )

cursor = connection.cursor()
stmt = "select geneName from promoters"
cursor.execute(stmt)

geneNames = [ x[0] for x in cursor.fetchall() ]

connection.close()
with open("geneNames.pickle", "wb") as pickout:
    pickle.dump(geneNames, pickout)
with open("geneNames.pickle", "rb") as pickin:
    GN = pickle.load(pickin)
    print("elt-2" in GN)
