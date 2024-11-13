import os,sys
import numpy as np

def isGoodRecoVertex( tree ):
    pass

def select_true_finalstates( true_photon_dict, true_vtx_info, true_protons, true_pions ):

    # in fiducial volume
    inFV = true_vtx_info["dwall"]>10.0

    # count number of detectable photons
    nphotons = 0    
    for i,data in enumerate(true_photon_dict):
        if data['detectable']==0:
            continue
        # how close to vertex, to be primary photon
        d = true_vtx_info["pos"][:3]-data["start"][:3]
        d = (d*d).sum()
        d = np.sqrt(d)
        print("detectable photon dist2vtx: ",d)
        if d<3.0:
            nphotons += 1
            
    # check for 1 proton above threshold
    nprotons = 0
    for i,data in enumerate(true_protons):
        KEproton = data["pmom"][3]-938.0
        if KEproton>=60.0:
            nprotons += 1

    # check for charged pions above threshold
    npions = 0
    for i,data in enumerate(true_pions):
        KEpion = data["pmom"][3]-139.0
        if KEpion>=60.0:
            npions += 1

    if inFV and npions==0:
        if nphotons==2 and nprotons>=1:    
            return "2gNp"
        elif nphotons==2 and nprotons==0:    
            return "2g0p"
        elif nphotons==1 and nprotons>=1:
            return "1gNp"
        elif nphotons==1 and nprotons==0:
            return "1g0p"

        
    return False
            





if __name__=="__main__":
    import ROOT as rt    

    #inputfile = "make_dlgen2_flat_ntuples_reco_v3lmshowerkp_gen2ntuple_evisvertex_output_0.root"
    inputfile = "make_dlgen2_flat_ntuples_reco_v3lmshowerkp_gen2ntuple_evisvertex_output_55.root"    
    rfile = rt.TFile(inputfile,"r")
    eventtree = rfile.Get("EventTree")
    pottree = rfile.Get("potTree")

    nentries = eventtree.GetEntries()

    N_true_2g1p = 0
    N_true_2g0p = 0
    N_true_1g1p = 0
    N_true_1g0p = 0

    final_states = ["2g0p","2gNp","1g0p","1gNp"]
    N_true_finalstates = {}
    N_reco_finalstates = {}    
    N_true_positives = {}
    N_false_positives = {}
    for fs in final_states:
        N_true_finalstates[fs] = 0
        N_reco_finalstates[fs] = 0        
        N_true_positives[fs] = 0
        N_false_positives[fs] = 0
    
    for ientry in range(nentries):

        print("======================================")
        print(" ENTRY[",ientry,"]")
        print("======================================")
              
              
        eventtree.GetEntry(ientry)

        trueVtxPos = np.array( (eventtree.trueVtxX, eventtree.trueVtxY, eventtree.trueVtxZ) )
        xdwall = min( np.abs(trueVtxPos[0]-0.0), np.abs(trueVtxPos[0]-256.0) )
        ydwall = min( np.abs(trueVtxPos[1]-106.0), np.abs(trueVtxPos[1]+107.0))
        zdwall = min( np.abs(trueVtxPos[1]-106.0), np.abs(trueVtxPos[1]+107.0))
        inTPC = False
        if ( trueVtxPos[0]>=0.0 and trueVtxPos[0]<=256.0
             and trueVtxPos[1]>-106.0 and trueVtxPos[1]<106.0
             and trueVtxPos[2]>0.0 and trueVtxPos[2]<1036.0):
            inTPC = True
        dwall_v = np.array( (xdwall, ydwall, zdwall) )
        dwall = np.min(dwall_v)
        true_vtx_info = {"pos":trueVtxPos,
                         "dwall_v":dwall_v,
                         "dwall":dwall}

        # true vertex
        print("True vertex:",true_vtx_info["pos"])
        print("in TPC: ",inTPC)
        print("dwall: ",dwall)


        # true photon information
        # confirm its NC.
        # get number of detectable photons and their start points
        true_det_photons = []
        true_protons = []
        true_pions = []
        true_above_threshold_tracks = []
        ndetectable_photons = 0
        for i in range(eventtree.nTrueSimParts):
            # common data
            true_part_data = {"trackid":eventtree.trueSimPartTID[i],
                              "motherid":eventtree.trueSimPartMID[i],
                              "pmom":np.array((eventtree.trueSimPartPx[i],
                                               eventtree.trueSimPartPy[i],
                                               eventtree.trueSimPartPz[i],
                                               eventtree.trueSimPartE[i])),
                              "start":np.array((eventtree.trueSimPartX[i],
                                                eventtree.trueSimPartY[i],
                                                eventtree.trueSimPartZ[i])),
                              "processid":eventtree.trueSimPartProcess[i],}
            
            if eventtree.trueSimPartPDG[i]==22:
                # a photon
                true_photon_data = true_part_data
                true_photon_data["planepixsum"] = np.array((eventtree.trueSimPartPixelSumUplane[i],
                                                            eventtree.trueSimPartPixelSumVplane[i],
                                                            eventtree.trueSimPartPixelSumYplane[i]))*0.0162
                true_photon_data["edeppos"] = np.array((eventtree.trueSimPartEDepX[i],
                                                        eventtree.trueSimPartEDepY[i],
                                                        eventtree.trueSimPartEDepZ[i]))
                                             
                
                """
                trueSimPartProcess: Integer indicating the process by which a detsim-tracked particle was created. 
                   0 for primary particles from neutrino interaction, 
                   1 for particles produced in a decay, 
                   2 for all other processes
                """
                
                nplanes_above = (true_photon_data["planepixsum"]>=15.0).sum()
                true_photon_data["nplanes_above"] = nplanes_above

                dist2vtx = (true_photon_data["edeppos"]-trueVtxPos)
                dist2vtx = np.sqrt( (dist2vtx*dist2vtx).sum() )
                true_photon_data["dist2vtx"] = dist2vtx
                
                if nplanes_above>=2:
                    true_photon_data["detectable"] = 1
                    ndetectable_photons += 1
                else:
                    true_photon_data["detectable"] = 0
                    
                true_det_photons.append(true_photon_data)
            elif eventtree.trueSimPartPDG[i]==2212:
                # check the vertex
                dist2vtx = (true_part_data["start"]-trueVtxPos)
                dist2vtx = np.sqrt( (dist2vtx*dist2vtx).sum() )
                true_part_data["dist2vtx"] = dist2vtx
                if dist2vtx<1.0:
                    true_protons.append( true_part_data )
            elif abs(eventtree.trueSimPartPDG[i])==211:
                # check the vertex
                dist2vtx = (true_part_data["start"]-trueVtxPos)
                dist2vtx = np.sqrt( (dist2vtx*dist2vtx).sum() )
                true_part_data["dist2vtx"] = dist2vtx
                if dist2vtx<1.0:
                    true_pions.append( true_part_data )
                
                
        print("Photon truth:")
        print("  number detectable: ",ndetectable_photons)
        print("")
        for i,data in enumerate(true_det_photons):
            print("  True Photon[",i,"]")
            for k,v in data.items():
                print("    ",k,": ",v)
        print("Proton truth:")
        for i,data in enumerate(true_protons):
            print("  True Proton[",i,"]")
            for k,v in data.items():
                print("    ",k,": ",v)
        print("Pion truth:")
        for i,data in enumerate(true_pions):
            print("  True Charged Pions[",i,"]")
            for k,v in data.items():
                print("    ",k,": ",v)

        # true selection
        true_final_state = select_true_finalstates( true_det_photons, true_vtx_info, true_protons, true_pions )
        if true_final_state is not False:
            if true_final_state=="2gNp":
                N_true_2g1p += 1
            elif true_final_state=="2g0p":
                N_true_2g0p += 1
            elif true_final_state=="1gNp":
                N_true_1g1p += 1
            elif true_final_state=="1g0p":
                N_true_1g0p += 1
        print("True final State: ",true_final_state)

        # ===============================================
        # ----------        
        #  Reco
        # ----------
        hasRecoVtx = eventtree.foundVertex==1
        vtxScore = eventtree.vtxScore
        inWCFV = eventtree.vtxIsFiducial
        recoVtxPos = np.array( (eventtree.vtxX,eventtree.vtxY,eventtree.vtxZ) )
        vtxdist = ((recoVtxPos-trueVtxPos)*(recoVtxPos-trueVtxPos)).sum()
        vtxDistToTrue = np.sqrt(vtxdist)        
        #kptype = eventtree.kpClusterType
        print("Reco Vertex:")
        print("  has reco vertex: ",hasRecoVtx)
        print("  vertex score (approx. vis energy): ",vtxScore," MeV")
        print("  inside WC FV: ",inWCFV)
        print("  Num Raw Tracks: ",eventtree.nTracks)
        print("  Num Raw ShowerS: ",eventtree.nShowers)
        print("  Vertex dist to true: ",vtxDistToTrue," cm")

        # Get Reco particle info
        reco_photons = []
        reco_electrons = []
        reco_protons = []
        reco_pions = []
        reco_muons = []
        
        # reco photon information
        print("Reco Shower info: ")
        for i in range(eventtree.nShowers):
            isSecondary = eventtree.showerIsSecondary[i]
            showerE = eventtree.showerCharge[i]*0.0162
            showerStartPos = np.array( (eventtree.showerStartPosX[i],eventtree.showerStartPosY[i],eventtree.showerStartPosZ[i]) )
            showerDir = np.array( (eventtree.showerStartDirX[i], eventtree.showerStartDirY[i], eventtree.showerStartDirZ[i]) )
            showerClassified = eventtree.showerClassified[i]
            showerLArPID = np.array( (eventtree.showerElScore[i], eventtree.showerPhScore[i], eventtree.showerMuScore[i], eventtree.showerPrScore[i], eventtree.showerPiScore[i]) )
            showerSoftMax = np.exp(showerLArPID)/np.sum(np.exp(showerLArPID))
            completeness = eventtree.showerComp[i]
            purity = eventtree.showerPurity[i]
            process = eventtree.showerProcess[i] # 0 primary, 1: secondary from neutral parent, 2: secondary from charged parent
            larpidprimary = eventtree.showerPrimaryScore[i]
            larpidFromNeutral = eventtree.showerFromNeutralScore[i]
            larpidFromCharged = eventtree.showerFromChargedScore[i]
            showerRecoE = eventtree.showerRecoE[i]
            
            dist2vtx = showerStartPos-recoVtxPos
            dist2vtx = np.sqrt( (dist2vtx*dist2vtx).sum() )

            # truth matching
            showerTrueTID = eventtree.showerTrueTID[i]
            showerTruePID = eventtree.showerTruePID[i]
            showerTrueCompleteness = eventtree.showerTrueComp[i]
            showerTruePurity       = eventtree.showerTruePurity[i]
            
            print("  RecoShower[",i,"]")            
            photon_data = {"E_MeV":showerRecoE,
                           "issecondary":isSecondary,
                           "isclassified":showerClassified,
                           "larpid_logits":showerLArPID,
                           "larpid":showerSoftMax,
                           "isprimary":larpidprimary,
                           "fromneutral":larpidFromNeutral,
                           "fromcharged":larpidFromCharged,
                           "comp":completeness,
                           "purity":purity,
                           "dir":showerDir,
                           "dist2vtx":dist2vtx,
                           "edeppos":showerStartPos,
                           "trueTID":showerTrueTID}
            for k,v in photon_data.items():
                print("    ",k,": ",v)
            print("   argmax(score):",np.argmax(photon_data["larpid"]))

            # keep as reco photon
            # basic definition for now
            if ( photon_data["E_MeV"]>15.0
                 and photon_data["purity"]>0.7
                 and photon_data["comp"]>0.7
                 and np.argmax(photon_data["larpid"])==1 ):
                print("    ** keep ** ")
                reco_photons.append( photon_data )
                

        #end of photon loop

        # analyze track data
        for i in range(eventtree.nTracks):
            track_data = {"start":np.array( (eventtree.trackStartPosX[i], eventtree.trackStartPosY[i], eventtree.trackStartPosZ[i]) ),
                          "dir":np.array( (eventtree.trackStartDirX[i], eventtree.trackStartDirY[i], eventtree.trackStartDirZ[i]) ),
                          "isClassified":eventtree.trackClassified[i],
                          "larpid_logits":np.array( (eventtree.trackElScore[i],
                                                     eventtree.trackPhScore[i],
                                                     eventtree.trackMuScore[i],
                                                     eventtree.trackPrScore[i],
                                                     eventtree.trackPiScore[i]) ),
                          "completeness":eventtree.trackComp[i],
                          "purity":eventtree.trackPurity[i],
                          "process":eventtree.trackProcess[i],
                          "isprimary_logit":eventtree.trackPrimaryScore[i],
                          "fromNeutral_logit":eventtree.trackFromNeutralScore[i],
                          "fromCharged_logit":eventtree.trackFromChargedScore[i],
                          "RecoE":eventtree.trackRecoE[i],
                          "truePID":eventtree.trackTruePID[i],
                          "trueTID":eventtree.trackTrueTID[i],
                          "trueE":eventtree.trackTrueE[i]}
            track_data["larpid"] = np.exp( track_data["larpid_logits"] )/np.sum( np.exp( track_data["larpid_logits"] ) )
            larpid = np.argmax( track_data["larpid"] )
            process_out = np.array( (eventtree.trackPrimaryScore[i], eventtree.trackFromNeutralScore[i], eventtree.trackFromChargedScore[i]) )
            process_normed = np.exp( process_out )/np.sum( np.exp(process_out) )
            track_data["isprimary"]   = process_normed[0]
            track_data["fromNeutral"] = process_normed[1]
            track_data["fromCharged"] = process_normed[2]
            track_data["pid_argmax"]  = larpid

            print("  Reco Track[",i,"]")
            for k,v in track_data.items():
                print("    ",k,": ",v)
            
            if track_data["isClassified"]==1:
                # defining observed  muons                
                if larpid==2 and track_data["process"]==0:
                    # muon candidate
                    KE = track_data["RecoE"]
                    if KE>60.0:
                        reco_muons.append(track_data)
                        print("    ** passing muon")
                # defining observed protons
                elif larpid==3 and track_data["process"]==0:
                    # proton candidate
                    KE = track_data["RecoE"]
                    if KE>60.0:
                        reco_protons.append(track_data)
                        print("    ** passing proton")
                # defining observed protons
                elif larpid==4 and track_data["process"]==0:
                    # pion candidate
                    KE = track_data["RecoE"]
                    if KE>60.0:
                        reco_pions.append(track_data)
                        print("    ** passing pion")

        # reco final state
        nphotons = len(reco_photons)
        nprotons = len(reco_protons)
        nmuons   = len(reco_muons)
        npions   = len(reco_pions)

        reco_final_state = False
        if nmuons==0 and npions==0:
            if nphotons==2 and nprotons==0:
                reco_final_state = "2g0p"
            elif nphotons==2 and nprotons>=1:
                reco_final_state = "2gNp"
            elif nphotons==1 and nprotons==0:
                reco_final_state = "1g0p"
            elif nphotons==1 and nprotons>=1:
                reco_final_state = "1gNp"

        print("TRUE FINAL STATE: ",true_final_state)
        print("RECO FINAL STATE: ",reco_final_state)
        print("  number of reco photons: ",nphotons)
        print("  number of reco protons: ",nprotons)
        print("  number of reco muons: ",nmuons)
        print("  number of reco pions: ",npions)                

        if true_final_state:
            N_true_finalstates[true_final_state] += 1
        if reco_final_state:
            N_reco_finalstates[reco_final_state] += 1
        if true_final_state and reco_final_state and true_final_state==reco_final_state:
            N_true_positives[reco_final_state] += 1
        if reco_final_state and reco_final_state!=true_final_state:
            N_false_positives[reco_final_state] += 1
            #input()
        
            

        # -------------------------------------------
        # evluation reco selection efficiency

        # -------------------------------------------
        # truth-evaluations

        # photon efficiency stats

        # photon purity stats

        # if Ng1p
        # proton efficiency stats
        # proton 
                
        if False and eventtree.nTracks>0:
            input()
            
        if False and ientry+1>=20:
            break

    print("Number of entries: ",nentries)
    print("Number of true 2g1p: ",N_true_2g1p)
    print("Number of true 2g0p: ",N_true_2g0p)
    print("Number of true 1g1p: ",N_true_1g1p)
    print("Number of true 1g0p: ",N_true_1g0p)

    for fs in final_states:
        print("FINAL STATE [",fs,"]")
        print("  Number of true ",fs,": ",N_true_finalstates[fs])
        print("  Number of reco ",fs,": ",N_reco_finalstates[fs])
        print("  Efficiency: ",float(N_true_positives[fs])/float(N_true_finalstates[fs]))
        print("  Purity: ",float(N_true_positives[fs])/float(N_reco_finalstates[fs]))
    
    print("FIN")
