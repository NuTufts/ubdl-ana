from __future__ import print_function
import os,sys,re

DLPRODDIR="/cluster/tufts/wongjiradlab/nutufts/dlgen2prod"

#samplename = "mcc9_v29e_dl_run3b_intrinsic_nue_LowE"
#inputlist="../maskrcnn_input_filelists/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlreco"

#samplename = "mcc9_v29e_dl_run3_G1_extbnb_dlana"
#inputlist="../maskrcnn_input_filelists/mcc9_v29e_dl_run3_G1_extbnb_dlana_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlana"

samplename = "mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
inputlist=DLPRODDIR+"/run3inputlists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list"
stem="merged_dlreco"

#samplename = "mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge"
#inputlist=DLPRODDIR+"/maskrcnn_input_filelists/mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge_MRCNN_INPUTS_LIST.txt"
#stem="merged_dlreco"

#samplename = "mcc9jan_run1_bnb5e19"
#inputlist="../run1inputlists/mcc9_v28_wctagger_bnb5e19_filelist.txt"
#stem="merged_dlreco"

outfolder="/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/extract_passing_events/output/%s/filtered_events/"%(samplename)
larmatch_outfolder="/cluster/tufts/wongjiradlab/nutufts/data/v0/%s/larmatch/"%(samplename)
reco_outfolder="/cluster/tufts/wongjiradlab/nutufts/data/v0/%s/larflowreco/ana/"%(samplename)

# GET INPUT LIST AND COMPILE HASH LIST
pnjobs = os.popen("cat %s | wc -l"%(inputlist)) # NUMBER OF SAMPLE FILES
njobs = int(pnjobs.readlines()[0])

fsourcelist = open(inputlist,'r')
lsourcelist = fsourcelist.readlines()
source_dict = {} # fileid:(hash,path)
fileid = 0
for fsource in lsourcelist:
    dlmerged = fsource.strip()
    # get dlmerged input
    h = dlmerged.split(stem+"_")[-1].split(".root")[0]
    source_dict[fileid] = (h,dlmerged)
    fileid += 1
    #print(h,os.path.basename(dlmerged))
print("Number of input files: ",fileid)
          
# NEED LIST OF LARFLOW RECO FILES
cmd = "find %s -name larflowreco_*.root -size +1k | sort" % (reco_outfolder)
#print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()

larflowreco_dict = {} # fileid:(hash, kpsreco file, larlite file)

for f in flist:
    f = f.strip()
    #print(f)
    base = os.path.basename(f)
    #print(base.split("fileid")[-1])
    x = re.split("[_-]+",base.split("fileid")[-1])
    #print(x)
    try:
        fileid = int(x[0])        
    except:        
        print("error parsing file: ",f," :: ",x)
        sys.exit(-1)
    
    larflowreco_dict[fileid] = (f,f.replace("kpsrecomanagerana","larlite").replace("/ana/","/larlite/"))
    #print(larflowreco_dict[fileid])
    #print(jobid," ",base)    

print("Number of larflow reco entries: ",len(larflowreco_dict))

# COMPILE LIST OF FINISHED FILES
filtered_dict = {} #fileid:(hash,filepath)
cmd = "find %s -name dlgen2filtered_*.root -size +1k | sort" % (outfolder)
#print(cmd)
plist = os.popen(cmd)
flist = plist.readlines()
for f in flist:
    f = f.strip()
    #print(f)
    base = os.path.basename(f)
    #print(base.split("fileid")[-1])
    x = re.split("[_-]+",base.split("fileid")[-1])
    #print(x)
    try:
        fileid = int(x[0])        
    except:        
        print("error parsing file: ",f," :: ",x)
        sys.exit(-1)
    
    filtered_dict[fileid] = (f)
    #print(filtered_dict[fileid])
    #print(jobid," ",base)    


print("Number of filterd files finished: ",len(filtered_dict))

# FOR EACH LARFLOW RECO FILE, CHECK IF OUTPUT EXISTS, IF NOT PREPARE ENTRY FOR RUNNING
NFINISHED = 0
runlist = []
for fileid,data in larflowreco_dict.items():
    if fileid in filtered_dict:
        continue # made it
    # compile inputs
    # kpsreco ana file
    # kpsreco larlite file
    # dlmerged file
    if fileid not in source_dict:
        print("Did not find FILEID[%d]"%(fileid))
        NFINISHED += 1
        continue
    kpsreco_ana = data[0]
    kpsreco_larlite = data[1]
    sourcedata = source_dict[fileid]

    # infer larmatch file
    larmatch = larmatch_outfolder + "/%03d/"%(fileid/100) + os.path.basename(kpsreco_larlite).replace("larflowreco","larmatch_kps")
    
    # check hash
    h = sourcedata[0]
    if h not in kpsreco_ana or h not in kpsreco_larlite:
        print("Entry mismatch!")
        print("  FILEID[%d]"%(fileid))
        print("  kpsreco-ana:     ",os.path.basename(kpsreco_ana))
        print("  kpsreco-larlite: ",os.path.basename(kpsreco_larlite))
        print("  dlmerged:        ",os.path.basename(sourcedata[1]))
        print("  larmatch:        ",os.path.basename(larmatch))
        sys.exit(0)
    runlist.append( (fileid,sourcedata[1],kpsreco_ana,kpsreco_larlite,larmatch) )

print("Number of jobs FINISHED: ",NFINISHED)
print("Number of jobs to run: ",len(runlist))

foutname="runlist_dlgen2filter_%s.txt"%(samplename)
print("Making file, %s, with file IDs to run"%(foutname))
fout = open(foutname,'w')
for (fileid,dlmerged,kpsana,kpslarlite,larmatch) in runlist:
    print("%d %s %s %s %s"%(fileid,kpsana,dlmerged,kpslarlite,larmatch),file=fout)
fout.close()


