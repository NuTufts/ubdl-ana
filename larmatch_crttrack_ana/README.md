# Scripts to make CRT Track Files

## Workflow

* make input file list in `inputlists` directory. One way is to use `find`.  Example:

    find [FILE_DIR]/ -name [file-pattern]*root | sort > inputlists/[sample-name].list

* if one needs to, run larmatch on the files. you might check if larmatch was already on this sample somewhere.
  First, update `submit_batch_kps_larmatch.sh` by setting: `RUN_DLANA_DIR`, `SAMPLE_NAME`.
  The `array` parameter indicates which job instances to run. 
  If you have `N` files in the inputlist, the job instances are between `0-[N-1]`.
  It is good to set this to simple `--array=0` to try one test job first.

* then run the CRT track reco on the larmatch files. 
  update `submit_batch_crttrack.sh` in a similar fashion to the previous step.


The output of the larmatch jobs should be in the `outputdir` folder. 
The output of the CRT jobs should be in the `crttrack_outputdir`. 
There are corresponding logfiles to check the jobs/debug.
    