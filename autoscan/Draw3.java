import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw3
{
    // Make a series of gjf files where the specified bond distance has been set and frozen.
    public static void main(String[] args)
        {
            // read input file
	    
            GJFfile input_gjf = new GJFfile("nico_Br_F_template.gjf");

            Molecule molecule = input_gjf.molecule;

            // define the atom numbers to be frozen
            // the fragment attached to "toAtomNumber" will be moved
            int fromAtomNumber1 = 16;
            int toAtomNumber1   = 18;
            int fromAtomNumber2 = 16;
            int toAtomNumber2   = 17;

            molecule = molecule.addBond(fromAtomNumber1, toAtomNumber1);
            molecule = molecule.addBond(fromAtomNumber2, toAtomNumber2);

            // define the keywords for the input files that will be generated
            String keywords = "%mem=3GB\n%nprocshared=4\n#p b3lyp 6-31+g* scrf=(pcm,solvent=thf) empiricaldispersion=gd3bj pop=none opt=modredundant freq=noraman";
            String tail = String.format("B %d %d F\nB %d %d F\n\n", fromAtomNumber1, toAtomNumber1, fromAtomNumber2, toAtomNumber2);
            
            for (double distance1 = 1.65; distance1 <= 3.20; distance1 += 0.05)
                {
                    for ( double distance2 = 1.30; distance2 <= 1.85; distance2 += 0.05)
                        {
                            Molecule newMolecule = molecule.setDistance(fromAtomNumber1, toAtomNumber1, distance1);
                            newMolecule = newMolecule.setDistance(fromAtomNumber2, toAtomNumber2, distance2);

                            // check the distance has been adjusted correctly
                            Atom atom1 = newMolecule.contents.get(fromAtomNumber1-1);
                            Atom atom2 = newMolecule.contents.get(toAtomNumber1-1);
                            Atom atom3 = newMolecule.contents.get(fromAtomNumber2-1);
                            Atom atom4 = newMolecule.contents.get(toAtomNumber2-1);

                            String filename = String.format("gjf/nico_Br_F-%03.0f_%03.0f-b3lyp_d3bj-631+gd-thf_pcm.gjf", distance1 * 100, distance2 * 100);
                            
                            // molecule, name, keywords, tail, writeConnectivity, charge, multiplicity
                            GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail, false, -1, 1);
                            output_gjf.write(filename);
                            
                            System.out.printf("%.2f\t%.2f\t%s\n", newMolecule.getDistance(atom1, atom2), newMolecule.getDistance(atom3, atom4), filename);
                        }
                }
        }
}
