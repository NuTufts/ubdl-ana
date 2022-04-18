from __future__ import print_function
import os,sys,re

"""
We need to make a text file that pairs [fileid] [input dlmerged file] [larmatch file]
"""

docheck = False

# MC SAMPLE
#samplename = "mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
#inputlist="pi0select_mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list"
#dl_folder="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/output_pi0_filtered"
#lm_folder="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/output_larmatch"
#reco_folder="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/output_reco/ana"
#input_stem="pi0select"

# DATA SAMPLE
samplename = "mcc9_v29e_dl_run3_pi0_lowBDT_sideband"
inputlist="mcc9_v29e_dl_run3_pi0_lowBDT_sideband.list"
lm_folder="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/output_larmatch/data/"
reco_folder="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/output_reco/data/ana"
input_stem="dlfilter_allsamples1"

# get list of finished reco files
cmd = "find %s -name larflowreco_file*.root -size +50k | sort" % (reco_folder)
print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()

finished = []
fmap = {}
# loop through output larmatch files, get fileid
for f in flist:
    f = f.strip()
    print(f)
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
    fmap[jobid] = f
    #print(jobid," ",base)

finished.sort()
print("Number of finished files: ",len(finished))
print("[enter] to continue")
input()

## Check if files are ok
fok = []
for f in finished:
    fname = fmap[f]
    isok = False
    if not docheck:
        isok = True
    else:
        #cmd = "python3 check_larmatch.py %s"%(fname)    
        cmd = "root -x -q \"check_larmatch_files.C(\\\"%s\\\")\""%(fname)
        print(cmd)
        pout = os.popen(cmd)
        lout = pout.readlines()
        nevents = int(lout[-1].strip().split()[-1])
        if nevents>0:
            isok = True
            
    if isok>0:
        fok.append(f)
        
print("checked ok",len(fok))
print("[enter] to continue")
#input()

fnjobs = open(inputlist,'r')
pinputlist=fnjobs.readlines()
njobs = len(pinputlist)

print("Number of files in input list: ",njobs)
missing = []
for i,inputfile in enumerate(pinputlist):
    inputfile = inputfile.strip()
    if i not in fok:
        missing.append((i,inputfile))
missing.sort()        
print("Number of missing files: ",len(missing))

foutname="larflowreco_runlist.txt"
print("Making file, %s, with file IDs to run"%(foutname))
fout = open(foutname,'w')
for (fileid,inputfile) in missing:
    lmfile = os.path.basename(inputfile).replace(input_stem,"larmatchme_fileid%04d"%(fileid)).replace(".root","_larlite.root")
    lmpath = "%s/%03d/%s"%(lm_folder,int(fileid/100),lmfile)
    if os.path.exists(lmpath):
        print("%d\t%s\t%s"%(fileid,inputfile,lmpath),file=fout)
    else:
        print("larmatch file does not exist for ",os.path.basename(inputfile))
        print(" lmpath=",lmpath)
        input()
fout.close()


