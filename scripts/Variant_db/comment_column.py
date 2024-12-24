#!/usr/bin/python3
import sqlite3
import sys

def add_comments_column(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    cursor.execute('ALTER TABLE Variants ADD COLUMN COMMENTS TEXT;')
    print("COMMENTS column added")

    conn.commit()
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: {} <database_name>".format(sys.argv[0]))
        sys.exit(1)

    db_name = sys.argv[1]
    add_comments_column(db_name)

