#!/bin/bash

set_common_vars() {

version="test"
case $1 in
    "run2" )
        config=run2_UL2018_nano_tau_v10_limited
        datasets='data_ul2018_a_single_mu,data_ul2018_b_single_mu,'`
        `'data_ul2018_c_single_mu,data_ul2018_d_single_mu'
        processes="data"
    ;;
    "run2_ul17_test" )
        config=run2_UL2017_nano_tau_v10_limited
        datasets='data_mu_f,wj_incl,dy_incl'
        processes="data_mu,wj,dy_lep,"
    ;;
    "run3lim")
        config="run3_2022_preEE_nano_tau_v12_limited"
        datasets='dy_incl,wj_incl,data_mu_c,data_mu_d' 
        #'st_t_bbarq,st_tbar_bq,'`
        #`'st_t_wminus_to_lnu2q,st_t_wminus_to_2l2nu,'`
        #`'st_tbar_wplus_to_lnu2q,st_tbar_wplus_to_2l2nu'
        processes='dy_lep,wj,data' #,
    ;;
    "run3")
        config="run3_2022_preEE_nano_tau_v12"
        #Datasets to use
        data='data_mu_c,data_mu_d,data_mu_e,'
        bkg_ewk='wj_incl,ww,wz,zz,dy_incl'
        bkg_top='st_t_bbarq,st_tbar_bq,'`
        `'st_t_wminus_to_lnu2q,st_t_wminus_to_2l2nu,'`
        `'st_tbar_wplus_to_lnu2q,st_tbar_wplus_to_2l2nu,'
        bkg_ttbar='tt_sl,tt_dl,tt_fh'
        datasets="$data$bkg_ewk$bkg_top$bkg_ttbar"
        processes="dy_z2mumu,dy_z2tautau,vv,tt,st,wj,data"
    ;;
    "run3_data")
        config="run3_2022_preEE_nano_tau_v12"
        #Datasets to use
        datasets="data_mu_c,data_mu_d"
        processes="data"
    ;;
    "run3_DY")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="dy_incl"
    processes="dy_z2tautau,dy_z2mumu,dy_z2ee"
    ;;
    # "run3_data_dy_wj_vv_tt_st")
    # config="run3_2022_preEE_nano_tau_v12"
    # #Datasets to use
    # datasets="st_tchannel_t,st_tchannel_tbar,st_twchannel_t_sl,st_twchannel_tbar_sl,st_twchannel_t_dl,st_twchannel_tbar_dl,st_twchannel_tbar_fh,tt_sl,tt_dl,tt_fh,ww,zz,wz,dy_incl,wj_incl,data_mu_c,data_mu_d"
    # processes="st,tt,vv,dy_lep,wj,data"
    # ;;
    "run3_data_dy_wj_vv_tt_st")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="data_mu_c,data_mu_d,dy_incl,wj_incl,ww,zz,wz,tt_sl,tt_dl,tt_fh,st_tchannel_t,st_tchannel_tbar,st_twchannel_t_sl,st_twchannel_tbar_sl,st_twchannel_t_dl,st_twchannel_tbar_dl,st_twchannel_tbar_fh"
    processes="st,tt,vv,dy_lep,wj,data"
    ;; 
    "run3_wj")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="wj_incl"
    processes="wj"
    ;;
    "run3_dy_wj")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="wj_incl,dy_incl"
    processes="wj,dy_lep"
    ;;
    "run3_vv")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="ww,zz,wz"
    processes="ww,zz,wz"
    ;;
    "run3_ww")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="ww"
    processes="ww"
    ;;
    "run3_ttbar_sl")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="tt_sl"
    processes="tt_sl"
    ;;
    "run3_ttbar_dl")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="tt_dl"
    processes="tt_dl"
    ;;
    "run3_ttbar_fh")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="tt_fh"
    processes="tt_fh"
    ;;
    "run3_tt")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="tt_sl,tt_dl,tt_fh"
    processes="tt_sl,tt_dl,tt_fh"
    ;;
    "run3_tt_fh")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="tt_fh"
    processes="tt_fh"
    ;;
    "run3_st")
    config="run3_2022_preEE_nano_tau_v12"
    #Datasets to use
    datasets="st_tchannel_t,st_tchannel_tbar,st_twchannel_t_sl,st_twchannel_tbar_sl,st_twchannel_t_dl,st_twchannel_tbar_dl,st_twchannel_tbar_fh"
    processes="st_tchannel_t,st_tchannel_tbar,st_twchannel_t_sl,st_twchannel_tbar_sl,st_twchannel_t_dl,st_twchannel_tbar_dl,st_twchannel_tbar_fh"
    ;;                       
    *)
    echo "Unknown run argument! Choose from: [run2, run3, run3lim, run2_ul17_test]"
    exit
esac
}
