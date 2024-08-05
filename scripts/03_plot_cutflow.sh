#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --processes $processes
        --version $version
        --datasets $datasets
        --workflow local
        --yscale log
        --print-output 1
        # --remove-output -1
        --selector-steps "trigger,met_filter,b_veto,dilepton_veto,extra_lepton_veto,trigobj_prematch,trigobj_postmatch,One_higgs_cand_per_event,has_proper_tau_decay_products"  
    )
law run cf.PlotCutflow "${args[@]}"