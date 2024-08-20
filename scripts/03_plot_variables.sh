#!/bin/bash
source ./common.sh #to access set_common_vars() function
#The following function defines config, processes, version and datasets variables
set_common_vars "$1"
args=(
        --config $config
        --datasets $datasets
        --processes $processes
        --version $version
        --categories incl,mutau,etau,
        --variables hcand_mass,tau_1_phi,tau_1_eta,electron_1_pt,electron_1_phi,electron_1_eta,muon_1_pt,muon_1_phi,muon_1_eta
        --general-settings "cms-label=pw"
        "${@:2}"
    )
echo law run cf.PlotVariables1D "${args[@]}"
law run cf.PlotVariables1D "${args[@]}"

# tau_1_phi,tau_1_eta,electron_1_pt,electron_1_phi,electron_1_eta,muon_1_pt,muon_1_phi,muon_1_eta,