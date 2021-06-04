from __future__ import print_function
import os,sys,re


#samplename = "mcc9_v29e_dl_run3b_intrinsic_nue_LowE"
#inputlist="../maskrcnn_input_filelists/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlreco"

#samplename = "mcc9_v29e_dl_run3_G1_extbnb_dlana"
#inputlist="../maskrcnn_input_filelists/mcc9_v29e_dl_run3_G1_extbnb_dlana_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlana"

samplename = "mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
inputlist="/cluster/tufts/wongjiradlab/nutufts/dlgen2prod/run3inputlists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list"
stem="merged_dlreco"

#samplename = "mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge"
#inputlist="/cluster/tufts/wongjiradlab/nutufts/dlgen2prod/maskrcnn_input_filelists/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlreco"

#samplename = "mcc9jan_run1_bnb5e19"
#inputlist="../run1inputlists/mcc9_v28_wctagger_bnb5e19_filelist.txt"
#stem="merged_dlreco"

outfolder="/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/output/%s/showerdqdx/"%(samplename)
larmatch_outfolder="/cluster/tufts/wongjiradlab/nutufts/data/v0/%s/larmatch/"%(samplename)
kpsana_outfolder="/cluster/tufts/wongjiradlab/nutufts/data/v1/%s/larflowreco/ana/"%(samplename)

# get list of showerdqdx files (the ones we are trying to make here)
cmd = "find %s -name showerdqdx_*.root -size +1k | sort" % (outfolder)
plist = os.popen(cmd)
flist = plist.readlines()
finished_showerdqdx = []
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
    finished_showerdqdx.append(jobid)
    #print(jobid," ",base)    

finished_showerdqdx.sort()
print("Number of finished SHOWERDQDX files: ",len(finished_showerdqdx))
#input()


# get list of finished larflow reco files
cmd = "find %s -name larflowreco_*.root -size +1k | sort" % (kpsana_outfolder)
print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()
finished_kpsreco = []
kpsreco_dict = {}
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
    finished_kpsreco.append(jobid)
    kpsreco_dict[jobid] = f
    #print(jobid," ",base)    
finished_kpsreco.sort()
print("Number of finished LARFLOW KPSRECO files: ",len(finished_kpsreco))
input()

# need list of larmatch files
cmd = "find %s -name larmatch_kps*larlite.root -size +100k | sort"%(larmatch_outfolder)
print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()
lm_finished = {} # holds hash to filename
for f in flist:
    f = f.strip()
    #print(f)
    f1 = f.replace("_larlite.root","")    
    if samplename in ["mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge",
                      "mcc9_v29e_dl_run3_G1_extbnb_dlana",
                      "mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge",
                      "mcc9jan_run1_bnb5e19"]:
        h = f1.split("larmatch_kps_fileid")[-1].split("_")[-1]
    elif samplename in ["mcc9_v29e_dl_run3b_intrinsic_nue_LowE"]:
        h = f.split("larmatch_kps_")[-1].split("-jobid")[0]
    else:
        raise ValueError("sample, %s, not setup for parsing larmatch files"%(samplename))    
    #print(h,": ",os.path.basename(f))
    lm_finished[h] = f    
print("larmatch files with hashes: ",len(lm_finished))
input("[ENTER] to continue")
             

# get list of dlmerged input files
pnjobs = os.popen("cat %s | wc -l"%(inputlist))
njobs = int(pnjobs.readlines()[0])

print("Number of files in input list: ",njobs)
missing = []
for i in range(njobs):
    if i not in finished_showerdqdx:
        missing.append(i)
missing.sort()        
print("Number of missing files: ",len(missing))

foutname="runlist_reco_%s.txt"%(samplename)
print("Making file, %s, with file IDs to run"%(foutname))
fout = open(foutname,'w')
for n,fileid in enumerate(missing):
    # get dlmerged input
    if n%500==0:
        print("processing %d of %d"%(n,len(missing)))
    pout = os.popen( "sed -n %dp %s"%(fileid+1,inputlist))
    dlmerged = pout.readlines()[0].strip()
    h = dlmerged.split(stem+"_")[-1].split(".root")[0]
    #base = os.path.basename( dlmerged ).split(
    #h = re.split("[_-]+",base.split("fileid")[-1])    
    #print(h)
    if h in lm_finished and fileid in kpsreco_dict:
        lmfile = lm_finished[h]
        kpsreco = kpsreco_dict[fileid]        
        print("%d %s %s"%(fileid,dlmerged,kpsreco),file=fout)
    else:
        print("No LARMATCH file with hash %s"%(h))
fout.close()


