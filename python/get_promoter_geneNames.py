#!/usr/bin/env python3
# conda install mysql-connector-python
from sys import stderr
import pickle
import mysql.connector as mariadb

with mariadb.connect(
    host="129.82.125.11",
    port="3307",
    user="worm",
    database="NishimuraLab"
    ) as connection:

    cursor = connection.cursor()
    stmt = "select geneName,WBID from promoters"
    print("executing", '`' + stmt + '`', file=stderr, end="...")
    cursor.execute(stmt)
    print(" done.", file=stderr)

    geneNames = {}
    for fields in cursor:
        geneNames[ fields[0] ] = fields[1]

    print("writing pickle 'geneNames.pickle'", file=stderr, end="...")
    with open("geneNames.pickle", "wb") as pickout:
        pickle.dump(geneNames, pickout)
    print(" done.", file=stderr)

    # test opening pickle
    with open("geneNames.pickle", "rb") as pickin:
        GN = pickle.load(pickin)
        print('testing: is elt-2 from pickle?', file=stderr, end="...? ")
        print("elt-2" in GN, file=stderr)

        print("GN['elt-2']:", GN['elt-2'])
