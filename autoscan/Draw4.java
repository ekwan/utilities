import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw4
{
    // Make a series of gjf files where the specified bond distance has been set and frozen.
    public static void main(String[] args)
        {
            // read input file
            GJFfile input_gjf = new GJFfile("template_galactose_b3lyp.gjf");
            Molecule molecule = input_gjf.molecule;

            // define the atom numbers to be frozen
            // the fragment attached to "toAtomNumber" will be moved
            int fromAtomNumber = 119;
            int toAtomNumber   = 125;

            molecule = molecule.addBond(fromAtomNumber, toAtomNumber);

            // define the keywords for the input files that will be generated
            String keywords = "%mem=8GB\n%nprocshared=12\n#p geom=connect b3lyp/6-31g(d) empiricaldispersion=gd3bj scrf=(pcm,solvent=toluene) pop=none freq=noraman opt=modredundant";
            String tail = String.format("B %d %d F\n\n", fromAtomNumber, toAtomNumber);
            
            for (double distance = 1.80; distance <= 2.20; distance += 0.05)
                {
                    Molecule newMolecule = molecule.setDistance(fromAtomNumber, toAtomNumber, distance);

                    // check the distance has been adjusted correctly
                    Atom atom1 = newMolecule.contents.get(fromAtomNumber-1);
                    Atom atom2 = newMolecule.contents.get(toAtomNumber-1);

                    String filename = String.format("output4/galactose_scan_b3lyp_d3bj_%3.0f.gjf", distance * 100);
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                    output_gjf.write(filename);
                    
                    System.out.printf("%.2f\t%s\n", newMolecule.getDistance(atom1, atom2), filename);
                }
        }
}
