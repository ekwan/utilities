import glob
import cctk
from cctk import GaussianFile, ConformationalEnsemble

filenames = sorted(list(glob.iglob("path*.gjf")))
ensemble = ConformationalEnsemble()
for filename in filenames:
    output_file = GaussianFile.read_file(filename)
    molecule = output_file.ensemble.molecules[-1]
    ensemble.add_molecule(molecule)
template = ensemble.molecules[0]
atoms = [1,2,3,4,5,6]
#ensemble = ensemble.align(to_geometry=0, comparison_atoms="heavy", compute_RMSD=False)
ensemble = ensemble.align(to_geometry=0, comparison_atoms=atoms, compute_RMSD=False)
template = ensemble.molecules[0]
for i,molecule in enumerate(ensemble.molecules):
    filename = f"aligned_path_{i+1:03d}.gjf"
    header = "#p"
    footer = "\n"
    GaussianFile.write_molecule_to_file(filename, molecule, header, footer)
