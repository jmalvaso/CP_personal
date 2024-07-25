#!/usr/bin/env python3
import matplotlib.pyplot as plt
import mplhep
import pickle
import numpy as np


def create_cutflow_histogram(cuts, data, xlabel="Selections", ylabel="Selection efficiency", title="", log=False, rel=False, save_path=None):
    """
    Create a cutflow histogram using Matplotlib with CMS style and save it to a PDF file.

    Parameters:
    - cuts: List of strings representing the names of the cuts.
    - data: List of integers representing the corresponding event counts for each cut.
    - xlabel: Label for the x-axis (default is "Cuts").
    - ylabel: Label for the y-axis (default is "Events").
    - title: Title of the plot (default is "Cutflow Histogram").
    - save_path: Path to save the PDF file. If None, the plot will be displayed but not saved.

    Returns:
    - fig: Matplotlib figure object.
    - ax: Matplotlib axis object.
    """

    # Set CMS style
    plt.style.use(mplhep.style.CMS)
   
    # Create Matplotlib figure and axis
    fig, ax = plt.subplots()
    if log: plt.yscale('log')
    # Plotting the cutflow histogram
    color = ['black','red']
    
    for i, (name, n_evt) in enumerate(data.items()):
        
        n_evt = np.array(n_evt)
        for cut_name, the_n_evt in zip(cuts,n_evt):
            print(f"{cut_name}: {the_n_evt}")
        if rel:
            x = cuts[1:]
            y = n_evt[1:]/n_evt[:-1]
        else:
            x = cuts
            y = n_evt/n_evt[0]
        print(f'Event numbers:')
        ax.scatter(x, y , color=color[i], marker='o', alpha=0.8, label=f"Data {name}")
        n_evt_prev = n_evt

    if log: ax.set_ylim((1*10**-2,2*10**0))
    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticklabels(cuts[1:] if rel else cuts, rotation=45, ha='right')
    
    # Add legend
    ax.legend()
    ax.grid(True, which='both', axis='y')
    ax.grid(True, which='major', axis='x')
    label_options = {
    "wip": "Work in progress",
    "pre": "Preliminary",
    "pw": "Private work",
    "sim": "Simulation",
    "simwip": "Simulation work in progress",
    "simpre": "Simulation preliminary",
    "simpw": "Simulation private work",
    "od": "OpenData",
    "odwip": "OpenData work in progress",
    "odpw": "OpenData private work",
    "public": "",
    }
    cms_label_kwargs = {
        "ax": ax,
        "llabel": label_options.get("pw"),
        "fontsize": 22,
        "data": True,
        'rlabel': "Data to Data"
    }
    mplhep.cms.label(**cms_label_kwargs)
   
    plt.tight_layout()
    
     # Save to PDF if save_path is provided
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    else:
        # Show the plot if save_path is not provided
        plt.show()


    return fig, ax

def get_hist_values(pickle_file):
    file_obj = open(pickle_file, 'rb')
    data = pickle.load(file_obj)
    hist = data.profile(axis=0)
    cuts = []
    values = []
    print(hist)
    for cut_name in hist.axes[1]:
        cuts.append(f'{cut_name}')
        values.append(hist[0,f'{cut_name}',0,0].count)
    return cuts, values

 



path22 = "/afs/cern.ch/user/j/jmalvaso/CP_CF_old/CPinHToTauTau/data/cf_store/analysis_httcp/cf.CreateCutflowHistograms/run3_2022_preEE_nano_tau_v12/data_mu_c/nominal/calib__main/sel__main__steps_trigger_met_filter_b_veto_PreTrigObjMatch_PostTrigObjMatch_One_higgs_cand_per_event_extra_lepton_veto/test/cutflow_hist__event.pickle"
cuts, values22 = get_hist_values(path22)

create_cutflow_histogram(cuts, 
                         data={"2022_postEE": values22},
                         save_path="mutau_cutflow_histogram_l.pdf",
                         log=False,
                         rel=True)