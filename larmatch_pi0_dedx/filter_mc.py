import os,sys,argparse
import ROOT as rt
from larcv import larcv
from larlite import larlite

parser = argparse.ArgumentParser("Filter MC")
parser.add_argument("--dlmerged",type=str,required=True)
parser.add_argument("--filelist",type=str,required=True)
args = parser.parse_args()

flist = open(args.filelist,'r')
flines = flist.readlines()
entry_list = []
for l in flines:
    l = l.strip()
    info = l.split(",")
    r = int(info[0])
    s = int(info[1])
    e = int(info[2])
    entry_list.append( (r,s,e) )
entry_list.sort()
#print(entry_list)

froot = rt.TFile( args.dlmerged, "open" )
larlite_id_tree = froot.Get("larlite_id_tree")
nentries = larlite_id_tree.GetEntries()

save_entry = []
for i in range(nentries):
    larlite_id_tree.GetEntry(i)
    r = larlite_id_tree._run_id
    s = larlite_id_tree._subrun_id
    e = larlite_id_tree._event_id

    rse = (r,s,e)
    #print(rse)
    if rse in entry_list:
        # save this event!
        save_entry.append(i)
froot.Close()

#save_entry = [3,8]

print("Entries to save: ",save_entry)
if len(save_entry)==0:
    sys.exit(0)


# we have an entry!
# open with larcv iomanager and larlite storage_manager
# save key trees

llout  = os.path.basename(args.dlmerged).replace("merged_dlreco","pi0select").replace(".root","_larlite.root")
lcvout = os.path.basename(args.dlmerged).replace("merged_dlreco","pi0select").replace(".root","_larcv.root")

ioll = larlite.storage_manager( larlite.storage_manager.kBOTH )
ioll.add_in_filename( args.dlmerged )
ioll.set_out_filename( llout )
ioll.set_data_to_read( larlite.data.kLArFlow3DHit, "larmatch" )
ioll.set_data_to_read( larlite.data.kMCTrack,  "mcreco" )
ioll.set_data_to_read( larlite.data.kMCShower, "mcreco" )
ioll.set_data_to_read( larlite.data.kMCTruth,  "generator" )
ioll.set_data_to_read( larlite.data.kOpFlash,  "simpleFlashBeam" )
ioll.set_data_to_read( larlite.data.kOpFlash,  "simpleFlashCosmic" )
ioll.open()

iolcv = larcv.IOManager( larcv.IOManager.kBOTH, "larcv", larcv.IOManager.kTickBackward )
iolcv.add_in_file(   args.dlmerged )
iolcv.set_out_file(  lcvout )
iolcv.specify_data_read( larcv.kProductImage2D, "wire" );
iolcv.specify_data_read( larcv.kProductImage2D, "thrumu" );
iolcv.specify_data_read( larcv.kProductImage2D, "ancestor" );
iolcv.specify_data_read( larcv.kProductImage2D, "segment" );
iolcv.specify_data_read( larcv.kProductImage2D, "instance" );
iolcv.specify_data_read( larcv.kProductImage2D, "larflow" );
iolcv.specify_data_read( larcv.kProductChStatus, "wire" );
iolcv.specify_data_read( larcv.kProductImage2D, "ubspurn_plane0" )
iolcv.specify_data_read( larcv.kProductImage2D, "ubspurn_plane1" )
iolcv.specify_data_read( larcv.kProductImage2D, "ubspurn_plane2" )
iolcv.specify_data_read( larcv.kProductSparseImage, "sparseuresnetout" )
iolcv.initialize()

nentries = ioll.get_entries()
for entry in save_entry:
    ioll.go_to(entry)
    iolcv.read_entry(entry)
    iolcv.save_entry()
ioll.next_event()    
ioll.close()
iolcv.finalize()

    

