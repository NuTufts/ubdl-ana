from __future__ import print_function
import os,sys,json

#sample = "intrinsicnue"
sample = "bnbnu"

larbysdata = "/cluster/tufts/wongjiradlab/larbys/data/mcc9"

datajsonfile = {"intrinsicnue":larbysdata+"/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge/prodgenie_bnb_intrinsic_nue_overlay_run3b_ssnet_wc_v2_nocrtremerge_run3b_ssnet_merged_dlreco.json",
                "bnbnu":larbysdata+"/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge/prodgenie_bnb_overlay_run3b_ssnet_wc_v2_nocrtremerge_run3b_ssnet_merged_dlreco.json"}
filelist = {"intrinsicnue":"intrinsic_v1.txt",
            "bnbnu":"bnb_nu_overlay_v1.txt"}

jfile = open(datajsonfile[sample],'r')
data = json.load(jfile)

mcpot = {}

for entry in data:
    fname = entry['file_name']
    h = fname.split("_")[-1].split(".")[0]
    mcpot[str(h)] = entry['mc.pot']

print("Loaded JSON")

flist = open(filelist[sample],"r")
flines = flist.readlines()

nfiles = 0
nmissing = 0
totpot = 0.0
for f in flines:
    f = f.strip()
    h = os.path.basename(f).split("_")[2]
    #print h
    if h in mcpot:
        totpot += float(mcpot[h])
        nfiles += 1
    else:
        nmissing += 1

print("Total MC POT: ",totpot)
print("nfiles: ",nfiles)
print("nmissing: ",nmissing)
