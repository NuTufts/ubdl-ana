# Vertexing Analysis

The goal of these scripts are to evaluate the vertexing performance of the larflow reco. chain.

The quantities we are trying to measure:

1) efficiency of reconstructing neutrino vertices
2) rate of bad vertices on off-beam events

These quantities we want to condition on:

* neutrino energy
* hadronic visible energy
* vertex quality scores
* mode 
* final state topology


# Analysis Chain

The analysis chain as of 7/15/2020:

* larmatch, keypoint, affinity field networks applied
* KPSReco applied to output

To do for the reco chain:

* sparse_infill applied before
* clustering using Mask-RCNN
* 

# Workflow

1. define an input list
2. run the larmatch deploy script on each file (`submit_batch_kps_larmatch.sh`, `run_batch_kps_larmatch.sh`)
3. run the reco on each file (`submit_batch_kps_reco.sh`, `run_batch_kps_reco.sh`)
4. run the ana program on each file
5. condense the ana files into histograms

# Samples

Short description of samples found in `inputlists/`.

| inputlist | description |
| mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE | intrinsic nue interactions below 400 MeV |

# Scripts



