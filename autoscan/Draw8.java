import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw8
{
    // Make a series of gjf files where the specified bond distance has been set and frozen.
    public static void main(String[] args)
        {
            // read input file
	    
            GJFfile input_gjf = new GJFfile("clay-stripped-product-b3lyp_d3bj-631gd-thf_pcm.gjf");

            Molecule molecule = input_gjf.molecule;

            // define the atom numbers to be frozen
            // the fragment attached to "toAtomNumber" will be moved
            // 1-indexed
            int fromAtomNumber1 = 4;
            int toAtomNumber1   = 1;

            molecule = molecule.addBond(fromAtomNumber1, toAtomNumber1);

            // define the keywords for the input files that will be generated
            String keywords = "%mem=3GB\n%nprocshared=4\n#p b3lyp 6-31g* scrf=(pcm,solvent=thf) empiricaldispersion=gd3bj pop=none opt=modredundant freq=noraman";
            String tail = String.format("B %d %d F\n\n", fromAtomNumber1, toAtomNumber1);
            
            for (double distance1 = 1.30; distance1 <= 4.0; distance1 += 0.10)
                {
                    Molecule newMolecule = molecule.setDistance(fromAtomNumber1, toAtomNumber1, distance1);

                    // check the distance has been adjusted correctly
                    Atom atom1 = newMolecule.contents.get(fromAtomNumber1-1);
                    Atom atom2 = newMolecule.contents.get(toAtomNumber1-1);

                    String filename = String.format("gjf/clay-stripped-product_scan-%03.0f-b3lyp_d3bj-631gd-thf_pcm.gjf", distance1 * 100);
                    
                    // molecule, name, keywords, tail, writeConnectivity, charge, multiplicity
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail, false, 0, 1);
                    output_gjf.write(filename);
                    
                    System.out.printf("%.2f\t%s\n", newMolecule.getDistance(atom1, atom2), filename);
                }
        }
}
