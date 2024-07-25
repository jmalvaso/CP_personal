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
        --selector-steps "trigger,met_filter,b_veto,PreTrigObjMatch,PostTrigObjMatch,dilepton_veto,extra_lepton_veto,One_higgs_cand_per_event,has_proper_tau_decay_products'"  
    )
law run cf.PlotCutflow "${args[@]}"