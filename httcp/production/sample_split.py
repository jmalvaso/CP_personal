import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column,remove_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses={f"Tau.{var}" for var in [
                "decayMode", "genPartFlav"
                ] 
    } | {"process_id"},
    produces={
        "process_id"
    },
    mc_only=True,
)
def split_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    

    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }
    hcand_ele_mu_DM = events.hcand.decayMode[:,:1]
    # hcand_Tau_idx = events.hcand.rawIdx[:,1:2]
    # mask_hcand_tau_idx = ak.flatten(hcand_Tau_idx)
    match = events.Tau.genPartFlav
    mu2tau_fakes_mask = ((match == tau_part_flav["prompt_mu"]) | (match == tau_part_flav["tau->mu"])) 
    e2tau_fakes_mask = ((match == tau_part_flav["prompt_e"]) | (match == tau_part_flav["tau->e"]))
    genuine_tau_mask = (match == tau_part_flav["tau_had"])
    process_id = np.array(events.process_id,dtype=np.int64)
    # np.array(events.process_id,dtype=np.int64)*ak.ones_like(events.process_id)
    process_id = ak.where((mu2tau_fakes_mask & (hcand_ele_mu_DM == -2)), 51001, 51004)
    process_id = ak.where((e2tau_fakes_mask & (hcand_ele_mu_DM == -1)), 51003, process_id)
    events = remove_ak_column(events, "process_id")
    events = set_ak_column(events, "process_id", process_id, value_type=np.int32)
    return events


