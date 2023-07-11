#!/usr/bin/env python3
# conda install mysql-connector-python
import sys
import mysql.connector as mariadb

with mariadb.connect(
    host="129.82.125.11",
    port="3307",
    user="worm",
    database="williams2023"
    ) as connection:

    cursor = connection.cursor()
    stmt = "select distinct WBID from log2FoldChangeWide WHERE outcome_01 in ('equal','enriched')"
    print(stmt, file=sys.stderr)
    cursor.execute(stmt)

    output = cursor.fetchall()
    print("got", len(output), "records", file=sys.stderr)
    for row in output: print(row)

