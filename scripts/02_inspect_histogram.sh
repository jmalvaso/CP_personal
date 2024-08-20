#!/bin/bash
args=(
    /eos/user/j/jmalvaso/higgs_cp_store/analysis_httcp/cf.SelectEvents/run3_2022_preEE_nano_tau_v12/wj_incl/nominal/calib__main/sel__main/dev/columns_0.parquet
)

cf_inspect "${args[@]}"

