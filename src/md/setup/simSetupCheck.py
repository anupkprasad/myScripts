import MDAnalysis as mda
import numpy as np

def analyze_ligand_properties(path, residues=("resname ATP")):
    """
    Load a GROMACS system and compute ligand charge, total charge, bonds, and angles.

    Parameters
    ----------
    path : str
        Directory containing md.tpr and md.gro files.
    residues : tuple of str
        Residue names to analyze (default = common nucleotides and metal ions).

    Returns
    -------
    dict
        Dictionary containing charges, net charge, bonds, angles, and angle values.
    """
    # Load system
    u = mda.Universe(f"{path}/min.tpr", f"{path}/min.gro")

    # Select desired residues
    selection_query = residues

    print(f"Selection query: {selection_query}")
    selected = u.select_atoms(selection_query)

    # Extract data
    charges = selected.charges
    net_charge = np.sum(charges)
    bonds = selected.bonds
    angles = selected.angles

    # Print summary
    print(f"✅ Selected residues: {residues}")
    print(f"Number of atoms: {len(selected)}")
    print(f"Net charge: {net_charge:.3f} e")
    print(f"First 10 charges: {charges[:10]}")
    print(f"Number of bonds: {len(bonds)}")
    print(f"Number of angles: {len(angles)}")

    # Print example bonds
    if len(bonds) > 0:
        print("\n🔗 Example bonds (first 5):")
        for i, b in enumerate(bonds):
            atom_names = [u.atoms[idx].name for idx in b.indices]
            bond_length = b.length()
            print(f"  {i+1}. {atom_names[0]:>4s} — {atom_names[1]:>4s} : {bond_length:.3f} Å")

    # Print example angles with values
    if len(angles) > 0:
        print("\n📐 Example angles (first 10):")
        for i, a in enumerate(angles):
            atom_names = [u.atoms[idx].name for idx in a.indices]
            # Calculate angle value in degrees
            angle_value = a.angle()  # Returns value in degrees by default
            print(f"  {i+1}. {atom_names[0]:>4s} - {atom_names[1]:>4s} - {atom_names[2]:>4s} : {angle_value:.2f}°")
        

    return {
        "universe": u,
        "selection": selected,
        "charges": charges,
        "net_charge": net_charge,
        "bonds": bonds,
        "angles": angles,
        "angle_values": [a.angle() for a in angles] if len(angles) > 0 else []
    }

# Example usage
if __name__ == "__main__":
    path = "/mnt/mydrive/pseudokinase/pseudoNcontrol_AF3/0_MD_simulation/test/atp_p25092_kd1__amsa_atemp/equilibration"
    result = analyze_ligand_properties(path)
    
    # Access angle values if needed
    if result["angle_values"]:
        print(f"\n✅ Total angles analyzed: {len(result['angle_values'])}")


