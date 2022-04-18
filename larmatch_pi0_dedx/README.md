dE/dx study using showers from Pi0

This attempts to use the dE/dx from reconstructed showers.

# Workflow: MC

## Filter MC events

We want to compare data to MC where the data was selected by the DL Gen1 lowBDT+pi0selection by Katie Mason.
Katie has run the selection MC and so we use a list of the passing MC events to filter events to reduce the
amount of processing we need.

MC event list (Run 3): `ListforAaroosh_pi0selection_bnb_overlay_run3.txt`

The following script is used to filter the events: `filter_mc.py`.

We use the following scripts to launch jobs that run the filter: `run_filter_mc.sh` and `submit_filter.sh`

The current list of output filtered files are: `pi0select_mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list`

## Run larmatch

To run larmatch minkowski network, we use: `run_batch_larmatchminkowski.sh`

To make/update the input run list use: `gen_larmatch_runlist.py` (make sure to update it for making the MC list)

## Run Reco

To run the reco, we need to give both the (filtered) dlmerged file and the larmatch file.

To make a list that pairs the inputs, use: `gen_reco_runlist.py`

Then to run the reco, we use: `run_batch_larflowreco_mc_cpu.sh` and `submit_batch_reco.sh`

# Workflow: DATA

The events are already filtered. The list of data files is: `mcc9_v29e_dl_run3_pi0_lowBDT_sideband.list`

## Run larmatch

We use the script to run larmatch: `run_batch_larmatchminkowski.sh` (the same as MC)

To generate the input list for running: `gen_larmatch_runlist.py` (the same as MC, needs to be updated)


## Run Reco

To run the reco, we need to give both the (filtered) dlmerged file and the larmatch file.

To make a list that pairs the inputs, use: `gen_reco_runlist.py`

Then to run the reco, we use: `run_batch_larflowreco_data_cpu.sh` and `submit_batch_reco.sh`
