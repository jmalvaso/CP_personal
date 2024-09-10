# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from httcp.util import enforce_hcand_type,enforce_tauprods_type,trigger_object_matching,hlt_path_fired,has_pt_greater_equal,has_delta_r_less_equal,HLT_path_matching

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
    tau_etau_to_match = events.Tau[etau_tau_idx]
    
    mu_to_match = events.Muon[mutau_mu_idx]
    tau_mutau_to_match = events.Tau[mutau_tau_idx]
    
    tau1_tautau_to_match = events.Tau[tautau_tau1_idx]
    tau2_tautau_to_match = events.Tau[tautau_tau2_idx]
    
    ### Perfoming the trigger object matching is required###
    if domatch:
        hcand_array, etau_channel_mask, mutau_channel_mask, tautau_channel_mask =  HLT_path_matching(events, trigger_results, mu_to_match, e_to_match, tau_mutau_to_match, tau_etau_to_match, tau1_tautau_to_match, tau2_tautau_to_match, hcand_pair)
    else:
        hcand_array = hcand_pair
    
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
    
    ### Ensure that we save all the taus and we remove the duplicates
    tau_indices =  ak.concatenate([etau_tau_idx,mutau_tau_idx,tautau_tau1_idx,tautau_tau2_idx],axis=1)
    max_length = ak.max(ak.num(tau_indices, axis=1))
    tau_indices = ak.pad_none(tau_indices, target=max_length)
    tau_indices = ak.fill_none(tau_indices,-1)[:,:,None]
    tau_indices = ak.flatten(tau_indices, axis=2)
    tau_indices = tau_indices[tau_indices >= 0]

    hcand_electron_indices = ak.values_astype(etau_e_idx,np.int64)
    hcand_muon_indices = ak.values_astype(mutau_mu_idx,np.int64)
    hcand_tau_indices = ak.values_astype(tau_indices,np.int64)
    
    sel_hcand = ak.sum(ak.num(hcand_array_type_fix,axis=2),axis=1) == 2
    
    empty_hcand_array = hcand_array_type_fix[:,:0]
    
    hcand_array_skim = ak.where(sel_hcand,hcand_array_type_fix,empty_hcand_array)
     
    hcand_array_One_Higgs_cand = ak.flatten(hcand_array_skim,axis=2)

    
    events = set_ak_column(events, "hcand", hcand_array_One_Higgs_cand)

    return events,hcand_muon_indices,hcand_electron_indices,hcand_tau_indices,etau_channel_mask,mutau_channel_mask,tautau_channel_mask,hcand_array_One_Higgs_cand, SelectionResult(
        steps={
            "One_higgs_cand_per_event": sel_hcand,
        },
        # objects={
        #     "Muon": {
        #         "Muon": hcand_muon_indices,
        #     },
        #     "Electron": {
        #         "Electron": hcand_electron_indices,
        #     },
        #     "Tau": {
        #         "Tau": hcand_tau_indices,
        #     },
        # },
    )

def select_tauprods(hcand_idx, tauprods):
    

    field_type_dict = {
        'eta': 'float64',
        'pdgId': 'int64',
        'phi': 'float64',
        'pt': 'float64',
        'tauIdx': 'int64',
        'mass': 'float64',
        'charge': 'int64'
    }

    tauprods = enforce_tauprods_type(tauprods, field_type_dict)
    tauprods_idx = tauprods.tauIdx
    
    max_length = ak.max(ak.num(tauprods_idx, axis=1))
    tauprods_idx = ak.pad_none(tauprods_idx, target=max_length)
    hcand_idx = ak.pad_none(hcand_idx, target=max_length)
    hcandprod_mask = hcand_idx == tauprods_idx
    hcandprod_mask = ak.fill_none(hcandprod_mask,-1)

    hcandprods = tauprods[hcandprod_mask >= 0]

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
    etau_id   = self.config_inst.get_channel("etau").id
    mutau_id  = self.config_inst.get_channel("mutau").id
    tautau_id = self.config_inst.get_channel("tautau").id

    events   = self[assign_tauprod_mass_charge](events)

    tauprods = events.TauProd
    
    hcand1 = events.hcand[:,0:1]
    hcand2 = events.hcand[:,1:2]
    
    hcand1_idx = hcand1.rawIdx
    hcand2_idx = hcand2.rawIdx

    hcand1prods = ak.where(events.channel_id == tautau_id,
                           select_tauprods(hcand1_idx, tauprods), 
                           tauprods[:,:0])
    hcand2prods = ak.where(((events.channel_id == etau_id) | (events.channel_id == mutau_id) | (events.channel_id == tautau_id)), 
                           select_tauprods(hcand2_idx, tauprods),
                           tauprods[:,:0])

    dummy = (events.event >= 0)[:,None]
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
    
    hcand_prods = ak.concatenate([hcand1prods[:,None], hcand2prods[:,None]], axis=1)

    hcand_prods_array = enforce_hcand_type(ak.from_regular(hcand_prods),
                                           {"pt"            : "float64",
                                            "eta"           : "float64",
                                            "phi"           : "float64",
                                            "mass"          : "float64",
                                            "charge"        : "int64",
                                            "pdgId"         : "int64",
                                            "tauIdx"        : "int64"}
                                       )
    hcand_prods_array = ak.from_regular(hcand_prods)
    # events = set_ak_column(events, "hcand",     hcand_concat)
    events = set_ak_column(events, "hcandprod", hcand_prods_array)

    #for i in range(500): 
    #    print(f"ch : {events.channel_id[i]}\t{hcand1.decayMode[i]}\t{hcand2.decayMode[i]}\t{hcand1_mask[i]}\t{hcand2_mask[i]}\t{hcand1prods.pdgId[i]}\t{hcand2prods.pdgId[i]}")


    return events, SelectionResult(
        steps={
            "has_proper_tau_decay_products": ak.sum(hcand_prod_mask, axis=1) == 2,
        },
    )
