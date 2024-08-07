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
        --selector-steps "trigger,met_filter,b_veto,dilepton_veto,Hcand_creation,Hcand_Trigger_Macthing,One_higgs_cand_per_event,extra_lepton_veto,has_proper_tau_decay_products"  
        "${@:2}"
    )
law run cf.PlotCutflow "${args[@]}"