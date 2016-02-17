import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class SimpleDraw
{
    // Make a series of gjf files where the specified bond distance has been set and frozen.
    public static void main(String[] args)
        {
            // read input file
            GJFfile input_gjf = new GJFfile("diax_template.gjf");
            Molecule molecule = input_gjf.molecule;

            // define the atom numbers to be frozen
            // the fragment attached to "toAtomNumber" will be moved
            int fromAtomNumber = 138;
            int toAtomNumber   = 119;

            // these atoms should be connected in the template molecule
            // and there should be no loop from toAtomNumber back to fromAtomNumber
            if ( ! molecule.directlyConnected(fromAtomNumber, toAtomNumber) )
                throw new IllegalArgumentException("atoms must be directly connected");
            molecule.getHalfGraph(fromAtomNumber, toAtomNumber);

            // define the keywords for the input files that will be generated
            String keywords = "%chk=checkpoint.chk\n%mem=32GB\n%nprocshared=16\n#p geom=connect m062x/6-31g(d) scrf=(pcm,solvent=toluene) pop=none freq=noraman opt=modredundant";
            String tail = String.format("B %d %d F\n\n", fromAtomNumber, toAtomNumber);
            
            for (double distance = 2.22; distance <= 2.40; distance += 0.02)
                {
                    Molecule newMolecule = molecule.setDistance(fromAtomNumber, toAtomNumber, distance);

                    // check the distance has been adjusted correctly
                    Atom atom1 = newMolecule.contents.get(fromAtomNumber-1);
                    Atom atom2 = newMolecule.contents.get(toAtomNumber-1);

                    String filename = String.format("output/indocat_mannosyl_diax_BnOH_%.0f.gjf", distance * 100);
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                    output_gjf.write(filename);
                    
                    System.out.printf("%.2f\t%s\n", newMolecule.getDistance(atom1, atom2), filename);
                }
        }
}
