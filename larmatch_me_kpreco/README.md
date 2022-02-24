# Larmatch (ME) keypoint reco analysis

Scripts for analyzing the performance of the larmatch (MinkowskiEngine) keypoint reconstruction.

We are interested in:
 * what is the precision and recall for finding different keypoints.
 * what is the distance between the true and closest reco keypoint (of the proper type)
 * the performance metrics as a function of truth metadata, e.g. visible energy, interaction mode, prox to dead region.
 * also interested in event-level metrics
 * want false positive rate per EXTBNB (i.e. no beam) event

Files
 * `make_inputlists.py`: script to collate files for grid jobs
 * `run_kpreco_ana.cxx`: routine that takes in larmatch network output and produces analysis metric ROOT tree
 * `plot_results.ipynb`: jupyter notebook for result plotting
 * `jobscript.sh`: runs the analysis code on a worker node
 * `submit.sh`: job submission script

## Dependencies

We build against the UBDL repository. [Link](https://github.com/larbys/ubdl)

