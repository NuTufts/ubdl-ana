import os,sys

plines = os.popen("find output_pi0_filtered | grep .root")
llines = plines.readlines()

