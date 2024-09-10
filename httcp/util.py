# coding: utf-8

"""
Collection of helpers
"""

from __future__ import annotations


import law
import order as od
from typing import Any
from columnflow.util import maybe_import
from columnflow.columnar_util import ArrayFunction, deferred_column

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


@deferred_column
def IF_NANO_V9(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None


@deferred_column
def IF_NANO_V11(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version >= 10 else None


def transverse_mass(lepton: ak.Array, met: ak.Array) -> ak.Array:
    dphi_lep_met = lepton.delta_phi(met)
    mt = np.sqrt(2 * lepton.pt * met.pt * (1 - np.cos(dphi_lep_met)))
    return mt


def trigger_object_matching(
    vectors1: ak.Array,
    vectors2: ak.Array,
    threshold: float = 0.5,
    axis: int = 2,
) -> ak.Array:
    """
    Helper to check per object in *vectors1* if there is at least one object in *vectors2* that
    leads to a delta R metric below *threshold*. The final reduction is applied over *axis* of the
    resulting metric table containing the full combinatorics. When *return_all_matches* is *True*,
    the matrix with all matching decisions is returned as well.
    """
    # delta_r for all combinations
    dr = vectors1.metric_table(vectors2)
    # check per element in vectors1 if there is at least one matching element in vectors2
    any_match = (dr < threshold)
    #ak.any(, axis=axis)
    return any_match

# def hlt_path_fired(dictionary):
#     from IPython import embed; embed()
#     if len(dictionary) > 0:
#         max_length = 0
#         for key in dictionary.keys():
#             temp_length = ak.max(ak.num(dictionary[key], axis=1))
#             if temp_length > max_length: max_length = temp_length

#         hlt_condition = {}
#         for key in dictionary.keys():
#             hlt_condition[key] = ak.pad_none(dictionary[key], target=max_length)
#             hlt_condition[key] = ak.fill_none(hlt_condition[key],-1)[:,:,None]

#         hlt_condition_values = list(hlt_condition.values())
#         hlt_condition_values_concat = ak.concatenate(hlt_condition_values, axis=-1)
#         HLT_path_fired = ak.max(hlt_condition_values_concat, axis=-1)
#         return HLT_path_fired 

def hlt_path_fired(events,dictionary):
    HLT_path_fired = -1*ak.ones_like(ak.local_index(events.event), dtype=np.int64)
    if len(dictionary) > 0:
        # Start by collecting all the arrays from the dictionary
        # Collect arrays from the dictionary
        array_list = []
        for key in dictionary.keys():
            array_list.append(dictionary[key])

        # Stack all the arrays along the second axis (axis=1)
        hlt_condition_values_concat = np.stack(array_list, axis=1)
        HLT_path_fired = ak.fill_none(ak.max(hlt_condition_values_concat, axis=-1),-1)
    return HLT_path_fired 
    
# Define the fill_missing function
def fill_missing(arr, fill_value):
    return ak.fill_none(arr, fill_value)
   
def get_dataset_lfns(
        dataset_inst: od.Dataset,
        shift_inst: od.Shift,
        dataset_key: str,
) -> list[str]:
    # destructure dataset_key into parts and create the lfn base directory
    lfn_base = law.wlcg.WLCGDirectoryTarget(
        dataset_key,
        fs=f"local",
    )
    # loop though files and interpret paths as lfns
    paths = [lfn_base.child(basename, type="f").path for basename in lfn_base.listdir(pattern="*.root")]

    return paths


def getGenTauDecayMode(prod: ak.Array):
    pids = prod.pdgId

    is_ele  = np.abs(pids) == 11
    is_muon = np.abs(pids) == 13
    is_charged = ((np.abs(pids) == 211) | (np.abs(pids) == 321))
    is_neutral = ((pids == 111) | (pids == 311) | (pids == 130) | (pids == 310))

    edecay = ak.sum(is_ele,  axis=-1) > 0
    mdecay = ak.sum(is_muon, axis=-1) > 0
    hdecay = (ak.sum(is_charged, axis=-1) > 0) | (ak.sum(is_neutral, axis=-1) >= 0)

    Nc = ak.sum(is_charged, axis=-1)
    Np = ak.sum(is_neutral, axis=-1)

    dm = ak.where(edecay, 
                  -1, 
                  ak.where(mdecay, 
                           -2, 
                           ak.where(hdecay, 
                                    (5 * (Nc - 1) + Np),
                                    -9)
                       )
              )

    return dm


def enforce_hcand_type(hcand_pair_concat, field_type_dict):
    temp = {}
    for field, typename in field_type_dict.items():
        temp[field] = ak.enforce_type(ak.values_astype(hcand_pair_concat[field], typename), f"var * var * {typename}")
    hcand_array = ak.zip(temp)
    return hcand_array


# def enforce_hcand_type(hcand_array, field_type_dict):
#     temp = {}
    
#     for field, typename in field_type_dict.items():
#         # Extract the field from the collection
#         field_data = hcand_array[field]

#         # Handle empty lists by filling them with a default value
#         if 'float' in typename:
#             default_value = 0.0
#         elif 'int' in typename:
#             default_value = 0
#         else:
#             raise ValueError(f"Unsupported type: {typename}")
        
#         # Fill empty lists with the default value
#         filled_data = ak.fill_none(field_data, default_value, axis=-1)
        
#         # If the data is deeply nested, consider flattening
#         if isinstance(ak.type(filled_data), ak.types.ListType):
#             try:
#                 enforced_field = ak.enforce_type(ak.values_astype(filled_data, typename), f"var * {typename}")
#             except Exception as e:
#                 print(f"Failed to enforce type for {field}: {e}")
#                 enforced_field = filled_data  # Fallback to the filled data without enforcing
#         else:
#             enforced_field = ak.values_astype(filled_data, typename)
        
#         # Store the enforced field back in the temp dictionary
#         temp[field] = enforced_field
    
#     # Zip the fields back into a NanoCollection-like structure
#     hcand_array = ak.zip(temp)
    
#     return hcand_array


    
def enforce_tauprods_type(tauprods, field_type_dict):
    temp = {}
    
    # Iterate over each field and apply type enforcement
    for field, typename in field_type_dict.items():
        # Extract the field from the collection
        field_data = tauprods[field]
        
        # Convert the field to the specified type
        # Use ak.values_astype to ensure the correct type at each level
        # Enforce type with "var * typename"
        enforced_field = ak.enforce_type(ak.values_astype(field_data, typename), f"var * {typename}")
        
        # Store the enforced field back in the temp dictionary
        temp[field] = enforced_field
    
    # Zip the fields back into a NanoCollection-like structure
    tauprods = ak.zip(temp)
    
    return tauprods


# def has_matching_pt(vector1, vector2, threshold=1.0):
#     """
#     Checks if there is at least one object in vector1 that has a pt
#     difference with any object in vector2 not greater than the threshold.
    
#     Parameters:
#         vector1 (awkward.Array): The pt values of the first vector (nested).
#         vector2 (awkward.Array): The pt values of the second vector (nested).
#         threshold (float): The maximum allowable difference in pt.
    
#     Returns:
#         awkward.Array: A boolean array with True if at least one pair of objects
#                        has a pt difference not greater than the threshold, otherwise False.
#     """
#     # Cartesian product to get all pairs
#     cartesian_product = ak.cartesian({"v1": vector1, "v2": vector2}, axis=-1, nested=True)
    
#     # Extract pt differences
#     pt_diff_values = np.abs(cartesian_product["v1"].pt - cartesian_product["v2"].pt)
    
#     # Check if any of the differences are within the threshold
#     has_match = ak.any(pt_diff_values <= threshold, axis=-1)
    
#     return has_match

def has_pt_greater_equal(vector1, vector2, offset):
    """
    Checks if the pt of any object in vector1 is greater than or equal to 
    the pt of any object in vector2 plus an offset.
    
    Parameters:
        vector1 (awkward.Array): The pt values of the first vector (nested).
        vector2 (awkward.Array): The pt values of the second vector (nested).
        offset (float): The value to add to the pt of objects in vector2.
    
    Returns:
        awkward.Array: A boolean array with True if at least one object in vector1 
                       has pt greater than or equal to pt in vector2 + offset, otherwise False.
    """
    # Cartesian product to get all pairs
    cartesian_product = ak.cartesian({"v1": vector1, "v2": vector2}, axis=-1, nested=True)
    
    # Extract pt values and apply the condition
    pt_comparison = (cartesian_product["v1"].pt >= (cartesian_product["v2"].pt + offset))
    
    # Check if the condition is met for any pair
    has_match = pt_comparison
    
    return has_match

def delta_phi(phi1, phi2):
    delta_phi = phi1 - phi2
    return ak.where(delta_phi > np.pi, delta_phi - 2*np.pi, 
                    ak.where(delta_phi < -np.pi, delta_phi + 2*np.pi, delta_phi))

def delta_eta(eta1, eta2):
    return eta1 - eta2

def delta_r(vector1, vector2):
    dphi = delta_phi(vector1.phi, vector2.phi)
    deta = delta_eta(vector1.eta, vector2.eta)
    return np.sqrt(deta**2 + dphi**2)

def has_delta_r_less_equal(vector1, vector2, threshold):
    """
    Checks if the delta_R between any object in vector1 and any object in vector2
    is less than or equal to the threshold.
    
    Parameters:
        vector1 (awkward.Array): The first vector (nested) containing eta and phi.
        vector2 (awkward.Array): The second vector (nested) containing eta and phi.
        threshold (float): The delta_R threshold value.
    
    Returns:
        awkward.Array: A boolean array with True if at least one pair of objects
                       has delta_R less than or equal to the threshold, otherwise False.
    """
    # Cartesian product to get all pairs
    cartesian_product = ak.cartesian({"v1": vector1, "v2": vector2}, axis=-1, nested=True)
    
    # Calculate delta_R
    delta_r_values = delta_r(cartesian_product["v1"], cartesian_product["v2"])
    
    # Check if the delta_R is less than or equal to the threshold for any pair
    has_match = ak.any(delta_r_values < threshold, axis=-1)
    
    return has_match

def HLT_path_matching(events, triggers, mu_to_match, e_to_match, tau_mutau_to_match, tau_etau_to_match, 
                      tau1_tautau_to_match, tau2_tautau_to_match, hcand_pair):
    
    # Initialize masks and dictionaries
    false_mask = ak.zeros_like(ak.local_index(events.event), dtype=np.bool_)
    single_electron_triggered = false_mask
    cross_electron_triggered  = false_mask
    single_muon_triggered = false_mask
    cross_muon_triggered  = false_mask
    cross_tau_triggered  = false_mask
    
    hlt_path_fired_e     = {}
    hlt_path_fired_etau  = {}
    hlt_path_fired_mu    = {}
    hlt_path_fired_mutau = {}
    hlt_path_fired_tau   = {}

    # Perform each lepton election step separately per trigger
    for trigger, trigger_fired, leg_masks in triggers.x.trigger_data:
        
        is_single_el = trigger.has_tag("single_e")
        is_cross_el  = trigger.has_tag("cross_e_tau")
        is_single_mu = trigger.has_tag("single_mu")
        is_cross_mu  = trigger.has_tag("cross_mu_tau")
        is_cross_tau = trigger.has_tag("cross_tau_tau")
        
        if is_single_mu or is_cross_mu:
            muons = mu_to_match
            taus = tau_mutau_to_match
            
            if is_single_mu:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 13
                
                mu_matches_pt = has_pt_greater_equal(muons, events.TrigObj[leg_masks[0]], 0)
                dr_matching = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                
                single_mu_matches_leg0 = (mu_matches_pt & dr_matching)
                single_mu_matches_leg0 = ak.any(ak.flatten(single_mu_matches_leg0, axis=-1), axis=1)
                single_muon_triggered = ak.where(trigger_fired & is_single_mu, True, single_muon_triggered)
                hlt_path_fired_mu[trigger.hlt_field] = ak.where(single_mu_matches_leg0, trigger.id, -1)
            
            elif is_cross_mu:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 13
                assert abs(trigger.legs[1].pdg_id) == 15
                
                mu_matches_pt = has_pt_greater_equal(muons, events.TrigObj[leg_masks[0]], 0)
                dr_matching_mu = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                cross_mu_matches_leg0 = (mu_matches_pt & dr_matching_mu)
                
                tau_matches_pt = has_pt_greater_equal(taus, events.TrigObj[leg_masks[1]], 0)
                dr_matching_tau = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                cross_mu_tau_matches_leg1 = (tau_matches_pt & dr_matching_tau)
                
                cross_mu_tau_matched = (ak.any(ak.flatten(cross_mu_matches_leg0, axis=-1), axis=1) & 
                                        ak.any(ak.flatten(cross_mu_tau_matches_leg1, axis=-1), axis=1))
                cross_muon_triggered = ak.where(trigger_fired & is_cross_mu, True, cross_muon_triggered)
                hlt_path_fired_mutau[trigger.hlt_field] = ak.where(cross_mu_tau_matched, trigger.id, -1)

        if is_single_el or is_cross_el:
            electrons = e_to_match
            taus = tau_etau_to_match

            if is_single_el:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 11
                
                el_matches_pt = has_pt_greater_equal(electrons, events.TrigObj[leg_masks[0]], 0)
                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                
                single_e_matches_leg0 = (el_matches_pt & dr_matching_e)
                single_e_matches_leg0 = ak.any(ak.flatten(single_e_matches_leg0, axis=-1), axis=1)
                single_electron_triggered = ak.where(trigger_fired & is_single_el, True, single_electron_triggered)
                hlt_path_fired_e[trigger.hlt_field] = ak.where(single_e_matches_leg0, trigger.id, -1)

            elif is_cross_el:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 11
                assert abs(trigger.legs[1].pdg_id) == 15
                
                el_matches_pt = has_pt_greater_equal(electrons, events.TrigObj[leg_masks[0]], 0)
                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                cross_e_matches_leg0 = (el_matches_pt & dr_matching_e)
                
                tau_matches_pt = has_pt_greater_equal(taus, events.TrigObj[leg_masks[1]], 0)
                dr_matching_tau = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                cross_e_tau_matches_leg1 = (tau_matches_pt & dr_matching_tau)
                
                cross_e_tau_matched = (ak.any(ak.flatten(cross_e_matches_leg0, axis=-1), axis=1) & 
                                       ak.any(ak.flatten(cross_e_tau_matches_leg1, axis=-1), axis=1))
                cross_electron_triggered = ak.where(trigger_fired & is_cross_el, True, cross_electron_triggered)
                hlt_path_fired_etau[trigger.hlt_field] = ak.where(cross_e_tau_matched, trigger.id, -1)
        
        if is_cross_tau:
            assert trigger.n_legs == len(leg_masks) == 2
            assert abs(trigger.legs[0].pdg_id) == 15
            assert abs(trigger.legs[1].pdg_id) == 15
            
            taus1 = tau1_tautau_to_match
            taus2 = tau2_tautau_to_match
            
            tau1_matches_pt = has_pt_greater_equal(taus1, events.TrigObj[leg_masks[0]], 0)
            tau2_matches_pt = has_pt_greater_equal(taus2, events.TrigObj[leg_masks[1]], 0)
            
            dr_matching_tau1 = trigger_object_matching(taus1, events.TrigObj[leg_masks[0]])
            dr_matching_tau2 = trigger_object_matching(taus2, events.TrigObj[leg_masks[1]])
            
            cross_tau1_matches_leg0 = (tau1_matches_pt & dr_matching_tau1)
            cross_tau2_matches_leg1 = (tau2_matches_pt & dr_matching_tau2)
            cross_tau_tau_matched = (ak.any(ak.flatten(cross_tau1_matches_leg0, axis=-1), axis=1) & 
                                     ak.any(ak.flatten(cross_tau2_matches_leg1, axis=-1), axis=1))
            cross_tau_triggered = ak.where(trigger_fired & is_cross_tau, True, cross_tau_triggered)
            hlt_path_fired_tau[trigger.hlt_field] = ak.where(cross_tau_tau_matched, trigger.id, -1)
    
    # Process trigger IDs for each type of lepton pair
    triggerID_e     = hlt_path_fired(events, hlt_path_fired_e)
    triggerID_etau  = hlt_path_fired(events, hlt_path_fired_etau)
    triggerID_mu    = hlt_path_fired(events, hlt_path_fired_mu)
    triggerID_mutau = hlt_path_fired(events, hlt_path_fired_mutau)
    triggerID_tau   = hlt_path_fired(events, hlt_path_fired_tau)
    
    # Generate candidate pairs based on matching triggers
    empty_hcand_pair = hcand_pair[:, :0][:, None]
    etau_channel_mask = ((triggerID_e > 0) | (triggerID_etau > 0))
    mutau_channel_mask = ((triggerID_mu > 0) | (triggerID_mutau > 0))
    tautau_channel_mask = (triggerID_tau > 0)
    
    hcand_pair_etau   = ak.where(etau_channel_mask, hcand_pair[:, 0][:, None], empty_hcand_pair)
    hcand_pair_mutau  = ak.where(mutau_channel_mask, hcand_pair[:, 1][:, None], empty_hcand_pair)
    hcand_pair_tautau = ak.where(tautau_channel_mask, hcand_pair[:, 2][:, None], empty_hcand_pair)

    hcand_array = ak.concatenate([hcand_pair_etau, hcand_pair_mutau, hcand_pair_tautau], axis=1)

    return hcand_array, etau_channel_mask, mutau_channel_mask, tautau_channel_mask
