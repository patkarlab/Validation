#!/usr/bin/python3
import sqlite3
import sys

def remove_comment(db_name, chrom, pos, ref, alt, comment_to_remove):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    cursor.execute('SELECT COMMENTS FROM Variants WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (chrom, pos, ref, alt))
    existing_comments = cursor.fetchone()

    if existing_comments is not None and existing_comments[0] is not None:
        existing_comments = existing_comments[0]
        
        comments_list = existing_comments.split(', ')
        
        if comment_to_remove.lower() == 'all':
            comments_list = []
            print("All comments removed from variant") 
        elif comment_to_remove in comments_list:
            comments_list.remove(comment_to_remove)
            print(f"Comment '{comment_to_remove}' removed from variant")
        else:
            print(f"Comment '{comment_to_remove}' not found in variant")

        updated_comments = ', '.join(filter(None, comments_list))

        cursor.execute('UPDATE Variants SET COMMENTS=? WHERE CHROM=? AND POS=? AND REF=? AND ALT=?;', (updated_comments, chrom, pos, ref, alt))

    else:
        print("Variant not found or no comments to remove")

    conn.commit()
    conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: {} <db_name> <chrom> <pos> <ref> <alt> <comment_to_remove>".format(sys.argv[0]))
        sys.exit(1)

    db_name = sys.argv[1]
    chrom = sys.argv[2]
    pos = int(sys.argv[3])
    ref = sys.argv[4]
    alt = sys.argv[5]
    comment_to_remove = sys.argv[6]

    remove_comment(db_name, chrom, pos, ref, alt, comment_to_remove)

