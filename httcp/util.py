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


def TetraVec(arr: ak.Array) -> ak.Array:
    TetraVec = ak.zip({"pt": arr.pt, "eta": arr.eta, "phi": arr.phi, "mass": arr.mass},
                      with_name="PtEtaPhiMLorentzVector",
                      behavior=coffea.nanoevents.methods.vector.behavior)
    return TetraVec


def _invariant_mass(events: ak.Array) -> ak.Array:
    empty_events = ak.zeros_like(events, dtype=np.uint16)[:, 0:0]
    where = ak.num(events, axis=1) == 2
    events_with_2 = ak.where(where, events, empty_events)
    mass = ak.fill_none(ak.firsts((TetraVec(events_with_2[:, :1]) + TetraVec(events_with_2[:, 1:2])).mass), 0)
    return mass


def deltaR(objects1: ak.Array,
           objects2: ak.Array) -> ak.Array:
    dr = objects1.metric_table(objects2)
    #print(dr)
    return ak.fill_none(ak.firsts(ak.firsts(dr, axis=-1)),0)
    

def invariant_mass(events: ak.Array):
    empty_events = ak.zeros_like(events, dtype=np.uint16)[:, 0:0]
    where = ak.num(events, axis=1) == 2
    events_2 = ak.where(where, events, empty_events)
    mass = ak.fill_none(ak.firsts((1 * events_2[:, :1] + 1 * events_2[:, 1:2]).mass), 0)
    return mass


def new_invariant_mass(objs1: ak.Array, objs2: ak.Array):
    lep_lep = ak.concatenate([objs1, objs2], axis=1)
    empty_events = ak.zeros_like(lep_lep, dtype=np.float32)[:, 0:0]
    where = ak.num(lep_lep, axis=1) == 2
    events = ak.where(where, lep_lep, empty_events)
    mass = ak.fill_none(ak.firsts((1 * events[:, :1] + 1 * events[:, 1:2]).mass), 0)
    return mass


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
    #print(f"Vectors1: {vectors1[0,:]}")
    #print(f"Vectors2: {vectors2[0,:]}")
    # delta_r for all combinations
    dr = vectors1.metric_table(vectors2)
    #print(f"DeltaR: {dr[0,:]}")
    # check per element in vectors1 if there is at least one matching element in vectors2
    any_match = ak.any(dr < threshold, axis=axis)
    #print(f"any_match: {any_match[0,:]}")
    return any_match


def get_dataset_lfns(
        dataset_inst: od.Dataset,
        shift_inst: od.Shift,
        dataset_key: str,
) -> list[str]:
    # destructure dataset_key into parts and create the lfn base directory

    #dataset_id, full_campaign, tier = dataset_key.split("/")[1:]
    #main_campaign, sub_campaign = full_campaign.split("-", 1)
    #print(dataset_id, full_campaign, tier)
    lfn_base = law.wlcg.WLCGDirectoryTarget(
        #f"/store/{dataset_inst.data_source}/{main_campaign}/{dataset_id}/{tier}/{sub_campaign}/0",
        #f"/eos/cms/store/group/phys_tau/TauFW/nanoV10/UL2017/{dataset_id}/{full_campaign}/{tier}",
        #"/eos/cms/store/group/phys_tau/TauFW/nanoV10/Run2_2017/DY2JetsToLL_M-50-madgraphMLM/",
        #f"/eos/user/g/gsaha3/Exotic/HCP_Test/UL2017/{dataset_id}/{full_campaign}/{tier}",
        dataset_key,
        fs=f"local",
    )
    #print(f"/eos/cms/store/group/phys_tau/TauFW/nano/UL2017/{dataset_id}/{full_campaign}/{tier}")

    # loop though files and interpret paths as lfns
    paths = [lfn_base.child(basename, type="f").path for basename in lfn_base.listdir(pattern="*.root")]
    #print(paths)
    return paths
