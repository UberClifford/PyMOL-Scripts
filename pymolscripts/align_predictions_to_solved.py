### IMPORTS AND STATIC STUFF ###
from pymol import cmd
from pathlib import Path
import os
import pdb
import argparse
import pickle
import colorsys
import numpy as np

### COMMANDLINE ARGUMENTS ###

parser = argparse.ArgumentParser("Create PyMOL session file of predicted structures to solved structure")
parser.add_argument("-ip", required=True, action="store", dest="predicted_dir", type=Path, help="Path to directory with predicted structures.")
parser.add_argument("-is", required=True, action="store", dest="solved_file", type=Path, help="Path to solved structure.")
parser.add_argument("-o", required=True, action="store", dest="output_session", type=Path, help="File path to save PyMOL session file (.pse).")
parser.add_argument("-p_chain", required=True, action="store", dest="pred_chains", type=str, help="Chains to use for alignment for predicted structures, (e.g. 'A', 'AB'). For predicted structure these chains must have the same name.")
parser.add_argument("-s_chain", required=True, action="store", dest="solved_chains", type=str, help="Chains to use for alignment on solved structure. e.g. (e.g. 'A', 'AB')")
parser.add_argument("-rm_algn_p_chains", action="store_true", dest="rm_algn_p_chains", help="Remove predicted structure chain used for alignment, after aligning it to solved structure chain (default:False)")
parser.add_argument("-color_predstruc_conf", dest="color_predstruc_conf", default=None, help="Color predicted structures by confidence. Default is None (not do it)")
parser.add_argument("-no_norm", action="store_false", dest="norm_conf", help="If set, will not normalize confidence scores for predicted structures. Default is true min-max normalize")
parser.add_argument("-set_normmin", action="store", type=float, default=None, dest="min_score", help="If set, wil use this as the min value for min-max normalization. Default is to use min of confidences for predicted structures.")
parser.add_argument("-set_normmax", action="store", type=float, default=None, dest="max_score", help="If set, wil use this as the max value for min-max normalization. Default is to use maximum of confidences for predicted structures.")


args = parser.parse_args()

predicted_dir = args.predicted_dir
solved_file = args.solved_file
output_session= args.output_session
pred_chains  = args.pred_chains
solved_chains = args.solved_chains
rm_algn_p_chains = args.rm_algn_p_chains
color_predstruc_conf = args.color_predstruc_conf
norm_conf = args.norm_conf
min_score = args.min_score
max_score = args.max_score


### FUNCTIONS ###
def get_rgb(score, min_score=None, max_score=None):
    
    if min_score != None and max_score != None: norm = (score - min_score) / (max_score - min_score + 1e-10)  # Avoid division by zero
    else: norm = score
    
    return colorsys.hsv_to_rgb(0.7 * (1 - norm), 1, 1)  # Blue → Red

### MAIN ### 


## LOAD SOLVED STRUCTURE ##
solved_chains = [c for c in solved_chains]
solved_obj = "solved"
cmd.load(solved_file, solved_obj)

# Create selection string for chains in solved structure
solved_sel = f"{solved_obj} and chain {'+'.join(solved_chains)}"

## LOAD AND ALIGN PREDICTED STRUCTURES ##

pred_chains = [c for c in pred_chains]
structure_files = list(predicted_dir.glob("**/*.pdb"))
N = len(structure_files)
obj_names = [f"{i}_{structure_files[i].stem}" for i in range(N)]

# if get confidence scores coloring
if color_predstruc_conf:
    # if chai-1 predictions, confidence scores are stored in the same directory as .npz files
    
    if color_predstruc_conf == "chai":
        scores = [] 
        for i in range(N):
            structure_file = structure_files[i]
            structure_filename = structure_file.stem
            score_file = structure_file.parent / f"{structure_filename.replace("pred", "scores")}.npz"
            np_data = np.load(score_file)
            score = 0.2*np_data["ptm"].item() + 0.8*np_data["iptm"].item()
            scores.append(score)
        
        # create score dict
        scores = np.asarray(scores)
        structurefile_conf_lookup = {structure_files[i]:scores[i] for i in range(N)}

    # if something else specified, assume it is path to user python dict in pickle format, where {key:value} = {structure_file:score}
    else:
        with open(color_predstruc_conf, "rb") as infile:
            structurefile_conf_lookup = pickle.load(infile)

# align structures 
for i in range(N):
    structure_file = structure_files[i]
    obj_name = obj_names[i]
    # Load predicted structure
    cmd.load(structure_file, obj_name)

    # Create selection string for predicted structure's chains
    pred_sel = f"{obj_name} and chain {'+'.join(pred_chains)}"
    
    # Align predicted chains onto solved chains
    cmd.align(pred_sel, solved_sel)
    
    # remove predicted aligned chains
    if rm_algn_p_chains: cmd.remove(pred_sel)

    print(f"Aligned structures: {i+1}/{N}")

# color structures by confidence
if color_predstruc_conf is not None:
    scores = np.asarray(list(structurefile_conf_lookup.values()))
    
    # use min and max confidences for predicted structures 
    if min_score is None: min_score = np.min(scores)
    if max_score is None: max_score = np.max(scores)
    
    print(f"Minimum confidence value {min_score} Maximum confidence value {max_score}")

    for i in range(N):
        structure_file = structure_files[i]
        obj_name = obj_names[i]
        
        # get color for score
        score = structurefile_conf_lookup[structure_file]
        if norm_conf: r, g, b = get_rgb(score, min_score=min_score, max_score=max_score)
        
        else: r, g, b = get_rgb(score)
        # color structure
        color_name = f"score_color_{obj_name}"
        cmd.set_color(color_name, [r, g, b])
        cmd.color(color_name, obj_name)

        print(f"Colored predicted structures by confidence: {i+1}/{N}")

# === SAVE SESSION ===
cmd.save(output_session)
print(f"Session saved as: {output_session}")
