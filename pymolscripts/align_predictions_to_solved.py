from pymol import cmd
import os

# === USER SETTINGS ===
predicted_dir = "predicted_pdbs"  # directory with predicted PDBs
solved_pdb = "solved_structure.pdb"  # path to solved structure
antigen_chains = ["A"]  # list of antigen chain IDs to align on
output_session = "aligned_complexes.pse"

# === LOAD SOLVED STRUCTURE ===
solved_obj = "solved"
cmd.load(solved_pdb, solved_obj)

# Create selection string for antigen chains in solved structure
solved_sel = f"{solved_obj} and chain {'+'.join(antigen_chains)}"

# === LOAD AND ALIGN PREDICTED STRUCTURES ===
for fname in os.listdir(predicted_dir):
    if not fname.endswith(".pdb"):
        continue

    pdb_path = os.path.join(predicted_dir, fname)
    obj_name = os.path.splitext(fname)[0]
    
    # Load predicted structure
    cmd.load(pdb_path, obj_name)

    # Create selection string for predicted structure's antigen chains
    pred_sel = f"{obj_name} and chain {'+'.join(antigen_chains)}"
    
    # Align predicted antigen chains onto solved antigen chains
    cmd.align(pred_sel, solved_sel)

# === SAVE SESSION ===
cmd.save(output_session)
print(f"Session saved as: {output_session}")
