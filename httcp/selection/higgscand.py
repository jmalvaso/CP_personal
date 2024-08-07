# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.util import enforce_hcand_type,trigger_object_matching,hlt_path_fired

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")



@selector(
    produces={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.charge", "hcand.rawIdx", "hcand.decayMode",
        #"hcand.IPx", "hcand.IPy", "hcand.IPz",
    },
    exposed=False,
)
def higgscand(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        hcand_pair: ak.Array,
        domatch: Optional[bool] = False,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:

    #Extraction of the indices from Hcand_pair
    etau_objects_idx = hcand_pair[:,0].rawIdx
    etau_e_idx = etau_objects_idx[:,0:1]
    etau_tau_idx = etau_objects_idx[:,1:2]
    
    mutau_objects = hcand_pair[:,1].rawIdx
    mutau_mu_idx = mutau_objects[:,0:1]
    mutau_tau_idx = mutau_objects[:,1:2]
    
    tautau_objects = hcand_pair[:,2].rawIdx
    tautau_tau1_idx = tautau_objects[:,0:1]
    tautau_tau2_idx = tautau_objects[:,1:2]
    
    #Objects to match
    e_to_match = events.Electron[etau_e_idx]
    tau_etau_to_mach = events.Tau[etau_tau_idx]
    
    mu_to_match = events.Muon[mutau_mu_idx]
    tau_mutau_to_mach = events.Tau[mutau_tau_idx]
    
    tau1_tautau_to_match = events.Tau[tautau_tau1_idx]
    tau2_tautau_to_mach = events.Tau[tautau_tau2_idx]
    
    false_mask = (abs(events.event) < 0)
    single_electron_triggered = false_mask
    cross_electron_triggered  = false_mask
    single_muon_triggered = false_mask
    cross_muon_triggered  = false_mask
    cross_tau_triggered  = false_mask
    
    hlt_path_fired_e   = {}
    hlt_path_fired_mu  = {}   
    hlt_path_fired_tau = {}      

    if  domatch:
        # perform each lepton election step separately per trigger
        for trigger, trigger_fired, leg_masks in trigger_results.x.trigger_data:
            print(f"trigger: {trigger}")
            print(f"trigger_fired: {trigger_fired}")
            print(f"leg_masks:  {leg_masks}")
            
            is_single_el = trigger.has_tag("single_e")
            is_cross_el  = trigger.has_tag("cross_e_tau")
            is_single_mu = trigger.has_tag("single_mu")
            is_cross_mu  = trigger.has_tag("cross_mu_tau")
            is_cross_tau = trigger.has_tag("cross_tau_tau")
            
            if is_single_mu or is_cross_mu:
                mu_matches_leg0 = None
                # muon selection
                muons = mu_to_match
                taus = tau_mutau_to_mach
                # start per-muon mask with trigger object matching
                if is_single_mu:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 1
                    assert abs(trigger.legs[0].pdg_id) == 13
                    # match leg 0
                    mu_matches_leg0 = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                    single_muon_triggered = ak.where(trigger_fired & is_single_mu, True, single_muon_triggered)
                    hlt_path_fired_mu[trigger.hlt_field]= ak.where(mu_matches_leg0, trigger.id,-1)                    
                elif is_cross_mu:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 2
                    assert abs(trigger.legs[0].pdg_id) == 13
                    assert abs(trigger.legs[1].pdg_id) == 15
                    # match leg 0
                    mu_matches_leg0 = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                    tau_matches_leg1 = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                    cross_muon_triggered = ak.where(trigger_fired & is_cross_mu, True, cross_muon_triggered)
                    hlt_path_fired_mu[trigger.hlt_field]= ak.where(mu_matches_leg0 & tau_matches_leg1, trigger.id,-1)                    

            if is_single_el or is_cross_el:
                el_matches_leg0 = None
                # electron selection
                electrons = e_to_match
                taus = tau_etau_to_mach
                # start per-muon mask with trigger object matching
                if is_single_el:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 1
                    assert abs(trigger.legs[0].pdg_id) == 11
                    # match leg 0
                    el_matches_leg0 = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                    single_electron_triggered = ak.where(trigger_fired & is_single_el, True, single_electron_triggered)
                    hlt_path_fired_e[trigger.hlt_field]= ak.where(el_matches_leg0, trigger.id,-1)
                elif is_cross_el:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 2
                    assert abs(trigger.legs[0].pdg_id) == 11
                    assert abs(trigger.legs[1].pdg_id) == 15
                    # match leg 0
                    el_matches_leg0 = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                    tau_matches_leg1 = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                    cross_electron_triggered = ak.where(trigger_fired & is_cross_el, True, cross_electron_triggered)
                    hlt_path_fired_e[trigger.hlt_field]= ak.where(el_matches_leg0 & tau_matches_leg1, trigger.id,-1)

            if is_cross_tau:
                # catch config errors
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 15
                assert abs(trigger.legs[1].pdg_id) == 15
                taus1 = tau1_tautau_to_match 
                tau2 = tau2_tautau_to_mach
                # match both legs
                tau_matches_leg0 = trigger_object_matching(taus1, events.TrigObj[leg_masks[0]])
                tau_matches_leg1 = trigger_object_matching(taus2, events.TrigObj[leg_masks[1]])
                is_cross_tau = tau_matches_leg0 & tau_matches_leg1
                cross_tau_triggered = ak.where(trigger_fired & is_cross_tau, True, cross_tau_triggered)
                hlt_path_fired_tau[trigger.hlt_field]= ak.where(is_cross_tau, trigger.id,-1)
        
        triggerID_e   = hlt_path_fired(hlt_path_fired_e  ) if len(hlt_path_fired_e  ) > 0 else print("No electrons match any trigger or no electron trigger required") 
        triggerID_mu  = hlt_path_fired(hlt_path_fired_mu ) if len(hlt_path_fired_mu ) > 0 else print("No muons match any trigger or no muons trigger required")
        triggerID_tau = hlt_path_fired(hlt_path_fired_tau) if len(hlt_path_fired_tau) > 0 else print("No taus match any trigger or no taus trigger required")   
    
    
    
    # sel_hcand = ak.sum(ak.num(hcand_pair.pt, axis=-1), axis=1) == 2
    
    #channel_id 1,2,4 etau,mutau,tautau
    empty_hcand_pair = hcand_pair[:,:0][:,None]
    etau_channel_mask = (single_electron_triggered | cross_electron_triggered)
    mutau_channel_mask = (single_muon_triggered | cross_muon_triggered)
    tautau_channel_mask = (cross_tau_triggered)
    
    hcand_pair_etau   = ak.where(etau_channel_mask, hcand_pair[:,0][:,None], empty_hcand_pair)
    hcand_pair_mutau  = ak.where(mutau_channel_mask, hcand_pair[:,1][:,None], empty_hcand_pair)
    hcand_pair_tautau = ak.where(tautau_channel_mask, hcand_pair[:,2][:,None], empty_hcand_pair)
    
    
    hcand_array = ak.concatenate([hcand_pair_etau, hcand_pair_mutau, hcand_pair_tautau], axis=1)
    
    # #channel_id 3,5,6 etau_mutau,etau_tautau,mutau_tautau 
    # hcand_pair_concat = ak.where(events.channel_id == 3, 
    #                              ak.concatenate([hcand_pair[:,0][:,None], hcand_pair[:,1][:,None]], axis=1),
    #                              hcand_pair_concat)
    # hcand_pair_concat = ak.where(events.channel_id == 5,
    #                              ak.concatenate([hcand_pair[:,0][:,None], hcand_pair[:,2][:,None]], axis=1),
    #                              hcand_pair_concat)
    # hcand_pair_concat = ak.where(events.channel_id == 6, 
    #                              ak.concatenate([hcand_pair[:,1][:,None], hcand_pair[:,2][:,None]], axis=1),
    #                              hcand_pair_concat)

    hcand_array_type_fix = enforce_hcand_type(hcand_array, 
                                     {"pt"            : "float64",
                                      "eta"           : "float64",
                                      "phi"           : "float64",
                                      "mass"          : "float64",
                                      "charge"        : "int64",
                                      "decayMode"     : "int64",
                                      "rawIdx"        : "int64",
                                      #"IPx"           : "float64",
                                      #"IPy"           : "float64",
                                      #"IPz"           : "float64"
                                       })
        
    
    #Extraction of the indices from Hcand_pair after matching
    etau_objects_idx = hcand_array_type_fix[:,0].rawIdx
    etau_e_idx = etau_objects_idx[:,0:1]
    etau_tau_idx = etau_objects_idx[:,1:2]
    
    mutau_objects = hcand_array_type_fix[:,1].rawIdx
    mutau_mu_idx = mutau_objects[:,0:1]
    mutau_tau_idx = mutau_objects[:,1:2]
    
    tautau_objects = hcand_array_type_fix[:,2].rawIdx
    tautau_tau1_idx = tautau_objects[:,0:1]
    tautau_tau2_idx = tautau_objects[:,1:2]
    
    tau_indices =  ak.concatenate([etau_tau_idx,mutau_tau_idx,tautau_tau1_idx,tautau_tau2_idx],axis=1)
    
    
    ########################################################################
    # Custom function to remove duplicates while preserving order
    def remove_duplicates(array):
        return ak.Array([list(dict.fromkeys(subarray)) for subarray in array])

    # This is not the proper way to make sure that we are not double counting taus,
    #    but it works, even though it is very slow
    unique_tau_indices = remove_duplicates(tau_indices)
    ########################################################################
    
    hcand_electron_indices = ak.values_astype(etau_e_idx,np.int64)
    hcand_muon_indices = ak.values_astype(mutau_mu_idx,np.int64)
    hcand_tau_indices = ak.values_astype(unique_tau_indices,np.int64)
    
    sel_hcand = ak.sum(ak.num(hcand_array_type_fix,axis=2),axis=1) == 2
    
    empty_hcand_array = hcand_array_type_fix[:,:0]
    
    hcand_array_skim = ak.where(sel_hcand,hcand_array_type_fix,empty_hcand_array)
     
    hcand_array_One_Higgs_cand = ak.flatten(hcand_array_skim,axis=2)

    
    events = set_ak_column(events, "hcand", hcand_array_One_Higgs_cand)

    return events,hcand_muon_indices,hcand_electron_indices,hcand_tau_indices,etau_channel_mask,mutau_channel_mask,tautau_channel_mask,hcand_array_One_Higgs_cand, SelectionResult(
        steps={
            "One_higgs_cand_per_event": sel_hcand,
        },
        objects={
            "Muon": {
                "Muon": hcand_muon_indices,
            },
            "Electron": {
                "Electron": hcand_electron_indices,
            },
            "Tau": {
                "Tau": hcand_tau_indices,
            },
        },
    )

#########################################################################################################
def select_tauprods(hcand_idx, tauprods):
    hcand_idx_brdcst, tauprod_tauIdx = ak.broadcast_arrays(ak.firsts(hcand_idx,axis=1), tauprods.tauIdx)
    hcandprod_mask                   = tauprod_tauIdx == hcand_idx_brdcst
    hcandprods                       = tauprods[hcandprod_mask]

    return hcandprods


is_pion         = lambda prods : ((np.abs(prods.pdgId) == 211) | (np.abs(prods.pdgId) == 321))
is_photon       = lambda prods : prods.pdgId == 22
has_one_pion    = lambda prods : (ak.sum(is_pion(prods),   axis = 1) == 1)[:,None]
has_three_pions = lambda prods : (ak.sum(is_pion(prods),   axis = 1) == 3)[:,None]
has_photons     = lambda prods : (ak.sum(is_photon(prods), axis = 1) >  0)[:,None]
has_no_photons  = lambda prods : (ak.sum(is_photon(prods), axis = 1) == 0)[:,None]


@selector(
    uses={
        "TauProd.pdgId",
    },
    produces={
        "TauProd.mass", "TauProd.charge",
    },
    exposed=False,
)
def assign_tauprod_mass_charge(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    pionp  =  211
    pionm  = -211
    kaonp  =  321
    kaonm  = -321
    gamma  =  22

    mass = ak.where(np.abs(events.TauProd.pdgId) == pionp, 
                    0.13957, 
                    ak.where(np.abs(events.TauProd.pdgId) == kaonp,
                             0.493677,
                             ak.where(events.TauProd.pdgId == gamma,
                                      0.0, 0.0)) 
    )
    charge = ak.where(((events.TauProd.pdgId == pionp) | (events.TauProd.pdgId == kaonp)),
                      1.0,
                      ak.where(((events.TauProd.pdgId == pionm) | (events.TauProd.pdgId == kaonm)),
                               -1.0,
                               0.0)
                  )

    events = set_ak_column(events, "TauProd.mass", mass)
    events = set_ak_column(events, "TauProd.charge", charge)    

    return events



@selector(
    uses={
        "channel_id", "TauProd.*", assign_tauprod_mass_charge,
    },
    produces={
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.charge", "hcand.rawIdx", "hcand.decayMode",
        #"hcand.IPx", "hcand.IPy", "hcand.IPz",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass", "hcandprod.charge", "hcandprod.pdgId", "hcandprod.tauIdx",
        assign_tauprod_mass_charge,
    },
    exposed=False,
)
def higgscandprod(
        self: Selector,
        events: ak.Array,
        hcand_array: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    
    One_higgs_cand = ak.num(hcand_array,axis=-1) == 2
    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id

    events   = self[assign_tauprod_mass_charge](events)

    tauprods = events.TauProd[One_higgs_cand]
    hcand = hcand_array[:,0]
    hcand =  hcand[One_higgs_cand]
    
    hcand1 = hcand[:,0:1]
    hcand2 = hcand[:,1:2]
    # hcand1 = ak.firsts(hcand[:,:,0:1], axis=1)
    # hcand2 = ak.firsts(hcand[:,:,1:2], axis=1)
    hcand_concat = ak.concatenate([hcand1, hcand2], axis=1)
    
    
    hcand1_idx = hcand1.rawIdx
    hcand2_idx = hcand2.rawIdx

    hcand1prods = ak.where(events.channel_id == tautau_id,
                           select_tauprods(hcand1_idx, tauprods), 
                           tauprods[:,:0])
    hcand2prods = ak.where(((events.channel_id == etau_id) | (events.channel_id == mutau_id) | (events.channel_id == tautau_id)), 
                           select_tauprods(hcand2_idx, tauprods),
                           tauprods[:,:0])

    dummy = (events.event >= 0)[:,np.newaxis]
    hcand1_mask = ak.where(hcand1.decayMode == 0,
                           (has_one_pion(hcand1prods) & has_no_photons(hcand1prods)),
                           ak.where(((hcand1.decayMode == 1) | (hcand1.decayMode == 2)),
                                    (has_one_pion(hcand1prods) & has_photons(hcand1prods)),
                                    ak.where(hcand1.decayMode == 10,
                                             has_three_pions(hcand1prods),
                                             ak.where(hcand1.decayMode == 11,
                                                      (has_three_pions(hcand1prods) & has_photons(hcand1prods)),
                                                      dummy)
                                         )
                                )
                       )
    hcand2_mask = ak.where(hcand2.decayMode == 0,
                           (has_one_pion(hcand2prods) & has_no_photons(hcand2prods)),
                           ak.where(((hcand2.decayMode == 1) | (hcand2.decayMode == 2)),
                                    (has_one_pion(hcand2prods) & has_photons(hcand2prods)),
                                    ak.where(hcand2.decayMode == 10,
                                             has_three_pions(hcand2prods),
                                             ak.where(hcand2.decayMode == 11,
                                                      (has_three_pions(hcand2prods) & has_photons(hcand2prods)),
                                                      dummy)
                                         )
                                )
                       )
    
    hcand_prod_mask = ak.concatenate([hcand1_mask, hcand2_mask], axis=1)
    
    hcand_prods = ak.concatenate([hcand1prods[:, np.newaxis], hcand2prods[:, np.newaxis]], axis=1)

    # hcand_prods_array = enforce_hcand_type(ak.from_regular(hcand_prods),
    #                                        {"pt"            : "float64",
    #                                         "eta"           : "float64",
    #                                         "phi"           : "float64",
    #                                         "mass"          : "float64",
    #                                         "charge"        : "int64",
    #                                         "pdgId"         : "int64",
    #                                         "tauIdx"        : "int64"}
    #                                    )
    hcand_prods_array = ak.from_regular(hcand_prods)
    events = set_ak_column(events, "hcand",     hcand_concat)
    events = set_ak_column(events, "hcandprod", hcand_prods_array)

    #for i in range(500): 
    #    print(f"ch : {events.channel_id[i]}\t{hcand1.decayMode[i]}\t{hcand2.decayMode[i]}\t{hcand1_mask[i]}\t{hcand2_mask[i]}\t{hcand1prods.pdgId[i]}\t{hcand2prods.pdgId[i]}")


    return events, SelectionResult(
        steps={
            "has_proper_tau_decay_products": ak.sum(hcand_prod_mask, axis=1) == 2,
        },
    )
#########################################################################################################