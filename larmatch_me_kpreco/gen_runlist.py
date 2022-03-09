from __future__ import print_function
import os,sys,re

# BNB NU RUN3
#samplename = "mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
#inputlist="/cluster/tufts/wongjiradlab//nutufts/dlgen2prod/run3inputlists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list"
#stem="merged_dlreco"
#larflow_outdir="/cluster/tufts/wongjiradlabnu/nutufts/data/larmatch/larmatch_me_v2/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge/"
#ana_outdir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_me_kpreco/output/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"

# BNB NUE CC RUN3
#samplename = "mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge"
#inputlist="/cluster/tufts/wongjiradlab//nutufts/dlgen2prod/maskrcnn_input_filelists/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlreco"
#larflow_outdir="/cluster/tufts/wongjiradlabnu/nutufts/data/larmatch/larmatch_me_v2/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge/"
#ana_outdir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_me_kpreco/output/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge"

# BNB NUE CC RUN3
samplename = "mcc9_v29e_dl_run3_G1_extbnb_dlana"
inputlist="/cluster/tufts/wongjiradlab//nutufts/dlgen2prod/maskrcnn_input_filelists/mcc9_v29e_dl_run3_G1_extbnb_dlana_MRCNN_INPUTS_LIST.txt"
stem="merged_dlana"
larflow_outdir="/cluster/tufts/wongjiradlabnu/nutufts/data/larmatch/larmatch_me_v2/mcc9_v29e_dl_run3_G1_extbnb_dlana/"
ana_outdir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_me_kpreco/output/mcc9_v29e_dl_run3_G1_extbnb_dlana/"

# get list of finished larmatch files
cmd = "find %s -name larmatchme_*.root -size +1k | sort" % (larflow_outdir)
print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()

finished_larmatch = {}

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
    finished_larmatch[jobid] = f
    #print(jobid," ",base)    

finished = finished_larmatch.keys()
print("Number of finished larmatch files: ",len(finished))
for n,x in enumerate(finished):
    print(" [",x,"]: ",finished_larmatch[x])
    if n>=9:
        break
print("[enter] to continue")
input()

# get ana output files
ana_finished = []


pnjobs = os.popen("cat %s | wc -l"%(inputlist))
njobs = int(pnjobs.readlines()[0])

print("Number of files in input list: ",njobs)
missing = []
for i in range(njobs):
    if i not in finished_larmatch:
        continue
    if i not in ana_finished:
        missing.append(i)
missing.sort()        
print("Number of missing files: ",len(missing))

foutname="runlist_kprecoana_%s.txt"%(samplename)
print("Making file, %s, with file IDs to run"%(foutname))
fout = open(foutname,'w')
for fileid in missing:
    # get dlmerged input
    pout = os.popen( "sed -n %dp %s"%(fileid+1,inputlist))
    dlmerged = pout.readlines()[0].strip()
    h = dlmerged.split(stem+"_")[-1].split(".root")[0]
    lmfile = finished_larmatch[fileid]
    print("%d %s %s"%(fileid,dlmerged,lmfile),file=fout)
fout.close()


