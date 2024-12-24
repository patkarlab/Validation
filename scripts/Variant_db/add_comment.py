#!/usr/bin/python3
import sqlite3
import sys

def add_comment(db_name, chrom, pos, ref, alt, comments):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    cursor.execute('SELECT COMMENTS FROM Variants WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (chrom, pos, ref, alt))
    existing_comments = cursor.fetchone()

    if existing_comments[0] is not None:
        old_comment = existing_comments[0]
        comments = f"{old_comment}, {comments}"
        cursor.execute('UPDATE Variants SET COMMENTS=? WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (comments, chrom, pos, ref, alt))
        print("Comment updated", comments, existing_comments)
    else:
        comments = f"{comments}"
        cursor.execute('UPDATE Variants SET COMMENTS=? WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (comments, chrom, pos, ref, alt))
        print("Comment added", comments, existing_comments)

    conn.commit()
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: {} <db_name> <chrom> <pos> <ref> <alt> <comments>".format(sys.argv[0]))
        sys.exit(1)

    db_name = sys.argv[1]
    chrom = sys.argv[2]
    pos = int(sys.argv[3])
    ref = sys.argv[4]
    alt = sys.argv[5]
    comments = sys.argv[6]

    add_comment(db_name, chrom, pos, ref, alt, comments)
 
