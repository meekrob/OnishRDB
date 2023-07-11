#!/usr/bin/env python3
# conda install mysql-connector-python
import sys, getpass
import mysql.connector as mariadb

password = getpass.getpass("MariaDB password for root> ")

with mariadb.connect(
    host="129.82.125.11",
    port="3307",
    user="root",
    password=password,
    database="NishimuraLab"
    ) as connection:

    stmt = "show engine innodb status;"
    cursor = connection.cursor()
    cursor.execute(stmt)


    for fields in cursor:
        print(*fields)
