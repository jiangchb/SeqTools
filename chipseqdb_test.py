#
# A test harness for chipseqdb.py
#

from chipseqdb import *

con = build_db("test.db")
print dir(con)
print type(con)
print con