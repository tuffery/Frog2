import _pysssr
sssr = _pysssr.sssr

##print 1
##try:
##    sssr(1, 1)
##except TypeError, msg:
##    print "caught type error:", msg

##print 2
##try:
##    sssr(1, [1])
##except TypeError, msg:
##    print "caught type error:", msg

##print 3
##try:
##    sssr(1, [(1,2,3)])
##except TypeError, msg:
##    print "caught type error:", msg

import sys
while 1:
    i = sssr(4, [(1,2), (2,3), (3, 1)])
#    print sys.getrefcount(i) - 1
#    asdf
#print sssr(1, [])
#print sssr(4, [(1,2), (2,3), (3, 1)])

import sys
i = sssr(4, [(1,2), (2,3), (3, 1)])
print sys.getrefcount(i) - 1
print i
l = []
print sys.getrefcount(l) - 1
