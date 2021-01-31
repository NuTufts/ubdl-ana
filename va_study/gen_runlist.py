from __future__ import print_function
import os,sys,re


samplename = "mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE"
inputlist="inputlists/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE_MRCNN_INPUTS_LIST.txt"

outfolder="output/%s/vareco/larlite/"%(samplename)
larmatch_outfolder="/cluster/tufts/wongjiradlab/nutufts/data/v0/mcc9_v29e_dl_run3b_intrinsic_nue_LowE/larmatch/"

# get list of finished vareco files
cmd = "find %s -name vareco_*.root -size +1k | sort" % (outfolder)
print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()

finished = []

for f in flist:
    f = f.strip()
    #print(f)
    base = os.path.basename(f)
    #print(base.split("fileid")[-1])
    x = re.split("[_-]+",base.split("fileid")[-1])
    #print(x)
    try:
        jobid = int(x[0])        
    except:        
        print("error parsing file: ",f," :: ",x)
        sys.exit(-1)
    finished.append(jobid)
    #print(jobid," ",base)    

finished.sort()
print("Number of finished files: ",len(finished))

# need list of larmatch files
cmd = "find %s -name larmatch_kps*larlite.root -size +100k | sort"%(larmatch_outfolder)
print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()
lm_finished = {}
for f in flist:
    f = f.strip()
    #print(f)
    h = f.split("larmatch_kps_")[-1].split("-jobid")[0]
    print(h)
    lm_finished[h] = f

pnjobs = os.popen("cat %s | wc -l"%(inputlist))
njobs = int(pnjobs.readlines()[0])

print("Number of files in input list: ",njobs)
missing = []
for i in range(njobs):
    if i not in finished:
        missing.append(i)
missing.sort()        
print("Number of missing files: ",len(missing))

foutname="runlist_vareco_%s.txt"%(samplename)
print("Making file, %s, with file IDs to run"%(foutname))
fout = open(foutname,'w')
for fileid in missing:
    # get dlmerged input
    pout = os.popen( "sed -n %dp %s"%(fileid+1,inputlist))
    dlmerged = pout.readlines()[0].strip()
    h = dlmerged.split("merged_dlreco_")[-1].split(".root")[0]
    if h in lm_finished:
        lmfile = lm_finished[h]
        print("%d %s %s"%(fileid,dlmerged,lmfile),file=fout)
    else:
        print("No larmatch file for %s"%(os.path.basename(dlmerged)))
fout.close()


