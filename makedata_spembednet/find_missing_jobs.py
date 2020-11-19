import os,sys
import numpy as np

f = open("test",'r')
flines = f.readlines()

stride = 20
narrayids=100
talley = np.zeros(narrayids, dtype=np.int )

for ll in flines:
    ll = ll.strip()
    jobid = int(ll.split("jobid")[1].split("_")[0])
    #print ll,jobid
    arrayid = jobid/stride
    talley[arrayid] += 1

#missing_arrayid=np.argwhere(talley!=stride)
missing_arrayid=np.argwhere(talley==0).reshape(-1)
print "talley: ",talley
print "missing: ",missing_arrayid.tolist()
print "total: ",len(missing_arrayid)
commalist = ""
for x in missing_arrayid.tolist():
    commalist+="%d,"%(x)
print commalist
