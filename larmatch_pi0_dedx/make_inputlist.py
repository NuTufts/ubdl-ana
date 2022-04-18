import os,sys

DLMERGED_INPUTLIST = "mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list"

plines = os.popen("find $PWD/output_pi0_filtered | grep .root")
llines = plines.readlines()
pairs = []
for l in llines:
    l = l.strip()
    print(l)
    hashtag = os.path.basename(l).replace("pi0select_","").replace("_merged.root","")
    p = os.popen("grep %s %s"%(hashtag,DLMERGED_INPUTLIST))
    lp = p.readlines()
    if len(lp)==1:
        x = lp[0].strip()
        pairs.append( (l,x) )
        
out = open( "paired_"+DLMERGED_INPUTLIST, 'w' )
for (l,x) in pairs:
    print(l," ",x,file=out)
out.close()
