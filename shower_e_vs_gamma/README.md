#shower e vs gamma ana

Studies aimed at optiming e vs. gamma discrimination.

We use perfect reco showers to study variables, and then compare to reco showers.

Strategies to develop:

  1. better dedx measurement on trunk [done...at least additional iteration implemented in `larflow/Reco/ShowerdQdx.h/.cxx`]
  2. 3d-branching metrics
  3. missed 2nd shower pi0 metric
  4. shower likelihood

## Study of `larflow::reco::ShowerdQdx` class

* script to run dev jobs is `run_showerdqdx.sh`. It runs `larflow/Reco/test/dev_showerbilineardqdx.py`.
* script to launch jobs is `submit_showerdqdx.sh`
* script to analyze output is `plot_dqdx.py`



## Shower likelihood profile

Workflow:

  1. Make perfect reco file for perfect reco showers
  2. Run shower likelihood code to make profile histograms

