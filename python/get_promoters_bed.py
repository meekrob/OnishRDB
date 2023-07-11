#!/usr/bin/env python3
# conda install mysql-connector-python
import sys
import mysql.connector as mariadb

with mariadb.connect(
    host="129.82.125.11",
    port="3307",
    user="worm",
    database="NishimuraLab"
    ) as connection:

    cursor = connection.cursor()
    stmt = "select chrom, start, end, CONCAT(geneName, ':', WBID), 1000, strand from promoters order by chrom, start"
    print("executing", '`' + stmt + '`', file=sys.stderr, end="...")
    cursor.execute(stmt)
    print("done", file=sys.stderr)

    print("writing promoters.bed", file=sys.stderr, end="...")
    with open("promoters.bed", "w") as promoter_out:
        for fields in cursor:
            print(*fields, sep="\t", file=promoter_out)
    print("done", file=sys.stderr)
