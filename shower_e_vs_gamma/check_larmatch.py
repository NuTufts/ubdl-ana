from __future__ import print_function
import os,sys

import ROOT as rt

inputfile = rt.TFile(sys.argv[1],"read")
larlite_id_tree = inputfile.Get("larlite_id_tree")
nentries =  0
try:
    nentries = larlite_id_tree.GetEntries()
except:
    nentries = 0

if nentries>0:
    print("1")
else:
    print("0")
