# coding: utf-8

"""
Configuration of the higgs_cp analysis.
"""

import functools

import law
import order as od
from scinum import Number

from columnflow.util import DotDict, maybe_import, dev_sandbox
from columnflow.config_util import (
    get_root_processes_from_campaign, 
    add_category,
    verify_config_processes,
)

ak = maybe_import("awkward")


def add_run3_preEE(ana: od.Analysis,
                      campaign: od.Campaign,
                      config_name           = None,
                      config_id             = None,
                      limit_dataset_files   = None,) -> od.Config :

    # get all root processes
    procs = get_root_processes_from_campaign(campaign)
    
    # create a config by passing the campaign, so id and name will be identical
    cfg = ana.add_config(campaign,
                        name  = config_name,
                        id    = config_id)

    # gather campaign data
    year = campaign.x.year
    
    # add processes we are interested in
    process_names = [
        "data",
        "data_mu",
        "data_e",
        #Drell-Yan
        "dy_lep",
        "dy_z2ee",
        "dy_z2mumu",
        # "dy_z2tautau",
        "dy_z2ll",
        # "dy_lep_m10to50",
        #W + jets
        "wj",
        #diboson
        "vv", #diboson inclusive
        "ww",
        "wz",
        "zz",
        #ttbar 
        "tt", # ttbar inclusive
        "tt_sl",
        "tt_dl",
        "tt_fh",
        #single top
        "st",
        #single top t-channel
        "st_tchannel_t",
        "st_tchannel_tbar",
        # single top tW channel
        "st_twchannel_t_sl",
        "st_twchannel_tbar_sl",
        "st_twchannel_t_dl",
        "st_twchannel_tbar_dl",
        # "st_twchannel_t_fh",
        #"st_twchannel_tbar_fh",
        # signal
        "h_ggf_tautau"
    ]
    for process_name in process_names:
        # add the process
        proc = cfg.add_process(procs.get(process_name))
        if proc.is_mc:
            # Updated color mapping to avoid repetition and ensure unique colors
            if proc.name == "st"                  : proc.color1 = (63, 144, 218)  # Sky Blue
            if proc.name == "st_tchannel_t"       : proc.color1 = (63, 144, 218)  # Sky Blue
            if proc.name == "st_tchannel_tbar"    : proc.color1 = (87, 144, 252)  # Dodger Blue
            if proc.name == "st_twchannel_t_sl"   : proc.color1 = (146, 218, 221) # Pale Turquoise
            if proc.name == "st_twchannel_tbar_sl": proc.color1 = (148, 164, 162) # Cadet Grey
            if proc.name == "st_twchannel_t_dl"   : proc.color1 = (169, 107, 89)  # Rosy Brown
            if proc.name == "st_twchannel_tbar_dl": proc.color1 = (200, 73, 169)  # Medium Violet Red
            if proc.name == "st_twchannel_tbar_fh": proc.color1 = (131, 45, 182)  # Amethyst
            if proc.name == "tt"                  : proc.color1 = (255, 169, 14)  # Orange
            if proc.name == "tt_sl"               : proc.color1 = (255, 169, 14)  # Orange
            if proc.name == "tt_dl"               : proc.color1 = (248, 156, 32)  # Dark Golden Rod
            if proc.name == "tt_fh"               : proc.color1 = (228, 37, 54)   # Crimson Red
            if proc.name == "vv"                  : proc.color1 = (101, 99, 100)  # Charcoal
            if proc.name == "ww"                  : proc.color1 = (101, 99, 100)  # Charcoal
            if proc.name == "zz"                  : proc.color1 = (185, 172, 112) # Olive Drab
            if proc.name == "wz"                  : proc.color1 = (122, 33, 221)  # Blue Violet
            if proc.name == "dy_lep"              : proc.color1 = (156, 156, 161) # Dark Gray
            if proc.name == "wj"                  : proc.color1 = (255, 94, 2)    # Orange Red


        # configuration of colors, labels, etc. can happen here
       

    # add datasets we need to study
    dataset_names = [
        #data
        "data_e_C",
        "data_e_D",
        "data_mu_C",
        "data_mu_D",
        #Drell-Yan
        "dy_incl",
        # "dy_lep_m10to50",
        #W+jets
        "wj_incl",
        #Diboson
        "ww",
        "wz",
        "zz",
        #ttbar
        "tt_sl",
        "tt_dl",
        "tt_fh",
        #single top t-channel
        "st_tchannel_t",
        "st_tchannel_tbar",
        # single top tW channel
        "st_twchannel_t_sl",
        "st_twchannel_tbar_sl",
        "st_twchannel_t_dl",
        "st_twchannel_tbar_dl",
        # "st_twchannel_t_fh",
        # "st_twchannel_tbar_fh",
        # signal
        "signal"
        ]
    
    for dataset_name in dataset_names:
        # add the dataset
        dataset = cfg.add_dataset(campaign.get_dataset(dataset_name))

        # for testing purposes, limit the number of files to 1
        for info in dataset.info.values():
            if limit_dataset_files:
                info.n_files = min(info.n_files, limit_dataset_files) #<<< REMOVE THIS FOR THE FULL DATASET

    # verify that the root process of all datasets is part of any of the registered processes
    verify_config_processes(cfg, warn=True)

  
    from httcp.config.triggers import add_triggers_run3_2022_preEE
    add_triggers_run3_2022_preEE(cfg)
    
    from httcp.config.met_filters import add_met_filters
    add_met_filters(cfg)

    # default objects, such as calibrator, selector, producer, ml model, inference model, etc
    cfg.x.default_calibrator = "main"
    cfg.x.default_selector = "main"
    cfg.x.default_producer = "main"
    cfg.x.default_ml_model = None
    cfg.x.default_inference_model = "example"
    cfg.x.default_categories = ("incl",)
    cfg.x.default_variables = ("event","channel_id")
    cfg.x.default_weight_producer = "all_weights"

    # process groups for conveniently looping over certain processs
    # (used in wrapper_factory and during plotting)
    cfg.x.process_groups = {
        "diboson": ["ww", "wz", "zz"],
        "tt" : ["tt_sl","tt_dl","tt_fh"],
        "st" : ["st_tchannel_t","st_tchannel_tbar","st_twchannel_t_sl","st_twchannel_tbar_sl","st_twchannel_t_dl","st_twchannel_tbar_dl"]
    }

    # dataset groups for conveniently looping over certain datasets
    # (used in wrapper_factory and during plotting)
    cfg.x.dataset_groups = {}

    # category groups for conveniently looping over certain categories
    # (used during plotting)
    cfg.x.category_groups = {}

    # variable groups for conveniently looping over certain variables
    # (used during plotting)
    cfg.x.variable_groups = {}

    # shift groups for conveniently looping over certain shifts
    # (used during plotting)
    cfg.x.shift_groups = {}

    # selector step groups for conveniently looping over certain steps
    # (used in cutflow tasks)
    cfg.x.selector_step_groups = {
        "default": ["json", "met_filter", "dl_res_veto", "trigger", "lepton", "jet"],
    }
    #  cfg.x.selector_step_labels = {"Initial":0, 
    #                                "Trigger": , "Muon"}
     
    # whether to validate the number of obtained LFNs in GetDatasetLFNs
    # (currently set to false because the number of files per dataset is truncated to 2)
    cfg.x.validate_dataset_lfns = False

    # lumi values in inverse pb
    # https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2?rev=2#Combination_and_correlations
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#DATA_AN2
    #Only F and G eras
    cfg.x.luminosity = Number(7980, {
        "lumi_13p6TeV_2022": 0.014j,
        
    })
    

    # names of muon correction sets and working points
    # (used in the muon producer)
    #cfg.x.muon_sf_names = ("NUM_TightRelIso_DEN_TightIDandIPCut", f"{year}")

    # register shifts
    cfg.add_shift(name="nominal", id=0)
  
    cfg.x.deep_tau = DotDict.wrap({
        "tagger": "DeepTau2018v2p5",
        "vs_e"          : "VVLoose",
        "vs_mu"         : "Tight",
        "vs_jet"        : "Medium",
        "vs_e_jet_wps"  : {'VVVLoose'   : 1,
                           'VVLoose'    : 2,
                           'VLoose'     : 3,
                           'Loose'      : 4,
                           'Medium'     : 5,
                           'Tight'      : 6,
                           'VTight'     : 7,
                           'VVTight'    : 8},
        "vs_mu_wps"     : {'VLoose' : 1,
                           'Loose'  : 2,
                           'Medium' : 3,
                           'Tight'  : 4}
        })

    cfg.x.btag_working_points = DotDict.wrap(
        {
            2016 : {
                "deepjet": { #TODO: make a link to this numbers
                    "loose": 0.0532,
                    "medium": 0.3040,
                    "tight": 0.7476,
                },
                "deepcsv": {
                    "loose": 0.1355,
                    "medium": 0.4506,
                    "tight": 0.7738,
                },
            },
            2022 : {
                "deepjet" : { #https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
                    "loose": 0.0583,
                    "medium": 0.3086,
                    "tight": 0.7183,
                }
            }
                    
                
        },
    )
    
    import os
    external_path = os.path.join(os.environ.get('HTTCP_BASE'), "httcp/data/corrections")
    # cfg.x.external_files = DotDict.wrap({
    #     # lumi files
    #     "lumi": {
    #         "golden": ("/afs/cern.ch/user/j/jmalvaso/public/Cert_Collisions2022_355100_362760_GoldenJSON.txt", "v1"),  # noqa
    #         "normtag": ("/afs/cern.ch/user/j/jmalvaso/public/normtag_PHYSICSJSON.txt", "v1"),
    #         # "golden": (f"{external_path}/Cert_Collisions2022_355100_362760_GoldenJSON.txt", "v1"),  # noqa
    #         # "normtag": (f"{external_path}/normtag_PHYSICS.json", "v1"),
    #     },
    #     "pileup":{
    #         #"json": ("/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/PileUp/EFG/pileup_JSON.txt", "v1")
    #         # "data" : "/afs/cern.ch/work/d/dwinterb/public/Run3_corrections/pu_data_2022_preEE.root",
    #         # "mc"   : "/afs/cern.ch/work/d/dwinterb/public/Run3_corrections/pu_mc_2022.root"         
    #         "data" : "/afs/cern.ch/user/j/jmalvaso/public/Data_PileUp_2022_preEE.root",
    #         "mc"   : "/afs/cern.ch/user/j/jmalvaso/public/MC_PileUp_2022.root",
    #     },
    
    corr_dir = "/afs/cern.ch/user/j/jmalvaso/CP_personal/httcp/corrections"
    cfg.x.external_files = DotDict.wrap({
        # lumi files
        "lumi": {
            "golden": (f"{corr_dir}/Cert_Collisions2022_355100_362760_Golden.json", "v1"),  # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
            "normtag": (f"{corr_dir}/normtag_PHYSICS.json", "v1"), #/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags
        },
        "pu_sf": (f"{corr_dir}/jsonpog-integration/POG/LUM/2022_Summer22/puWeights.json.gz", "v1"),
        "muon_correction" : f"{corr_dir}/jsonpog-integration/POG/MUO/2022_Summer22/muon_Z.json.gz",
        "tau_correction" : f"{corr_dir}/tau_DeepTau2018v2p5_2022_preEE.json.gz"

    })


    # target file size after MergeReducedEvents in MB
    cfg.x.reduced_file_size = 512.0
    
    from httcp.config.variables import keep_columns
    keep_columns(cfg)

    # event weight columns as keys in an OrderedDict, mapped to shift instances they depend on
    #get_shifts = functools.partial(get_shifts_from_sources, cfg)
    cfg.x.event_weights = DotDict({
        "normalization_weight"  : [],
        "pu_weight"             : [],
        #"electron_weight_nom"       : [],
        "muon_weight_nom"           : [],
        "tau_weight_nom"             : [],
    })
    

    # versions per task family, either referring to strings or to callables receving the invoking
    # task instance and parameters to be passed to the task family
    def set_version(cls, inst, params):
        # per default, use the version set on the command line
        version = inst.version 
        return version if version else 'dev1'
            
        
    cfg.x.versions = {
        "cf.CalibrateEvents"    : set_version,
        "cf.SelectEvents"       : set_version,
        "cf.MergeSelectionStats": set_version,
        "cf.MergeSelectionMasks": set_version,
        "cf.ReduceEvents"       : set_version,
        "cf.MergeReductionStats": set_version,
        "cf.MergeReducedEvents" : set_version,
    }
    # channels
    cfg.add_channel(name="etau",   id=1)
    cfg.add_channel(name="mutau",  id=2)
    #cfg.add_channel(name="emu"  ,  id=3)
    cfg.add_channel(name="tautau", id=4)
    
    if cfg.campaign.x("custom").get("creator") == "desy":  
        def get_dataset_lfns(dataset_inst: od.Dataset, shift_inst: od.Shift, dataset_key: str) -> list[str]:
            # destructure dataset_key into parts and create the lfn base directory
            print(f"Creating custom get_dataset_lfns for {config_name}")   
            try:
               basepath = cfg.campaign.x("custom").get("location")
            except:
                print("Did not find any basebath in the campaigns")
                basepath = "" 
            lfn_base = law.wlcg.WLCGDirectoryTarget(
                f"{basepath}{dataset_key}",
                fs="wlcg_fs_eoscms_redirector",
            )
            print(f"lfn basedir:{lfn_base}")
            # loop though files and interpret paths as lfns
            return [
                lfn_base.child(basename, type="f").path
                for basename in lfn_base.listdir(pattern="*.root")
            ]
        # define the lfn retrieval function
        cfg.x.get_dataset_lfns = get_dataset_lfns
        # define a custom sandbox
        cfg.x.get_dataset_lfns_sandbox = dev_sandbox("bash::$CF_BASE/sandboxes/cf.sh")
        # define custom remote fs's to look at
        cfg.x.get_dataset_lfns_remote_fs =  lambda dataset_inst: "wlcg_fs_eoscms_redirector"
        
    # add categories using the "add_category" tool which adds auto-generated ids
    from httcp.config.categories import add_categories
    add_categories(cfg)
        
    from httcp.config.variables import add_variables
    add_variables(cfg)
    
    
