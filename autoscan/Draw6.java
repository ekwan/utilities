import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw6
{
    // Make a series of gjf files where the specified bond distance has been set and frozen.
    public static void main(String[] args)
        {
            // read input file
            GJFfile input_gjf = new GJFfile("pi_complex.gjf");
            Molecule molecule = input_gjf.molecule;

            // define the bond angle to be frozen
            int atomNumber1 = 14;
            int atomNumber2 = 16;
            int atomNumber3 = 12;
            //molecule = molecule.addBond(fromAtomNumber, toAtomNumber);

            // define the keywords for the input files that will be generated
            String keywords = "%mem=3GB\n%nprocshared=4\n#p geom=connect b3lyp/genecp empiricaldispersion=gd3bj scrf=(pcm,solvent=water) pop=none freq=noraman opt=modredundant";
            String tail = String.format("A %d %d %d F\n\n@/n/jacobsen_lab/ekwan/basis/basis1.bas\n\n", atomNumber1, atomNumber2, atomNumber3);
            
            for (double angle = 80.0; angle <= 130.0; angle += 2.0)
                {
                    Molecule newMolecule = molecule.setAngle(atomNumber1, atomNumber2, atomNumber3, angle);

                    String filename = String.format("gjf/pi_complex_autoscan_%03.0f.gjf", angle);
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail, false, 1, 1);  // false = no connectivity
                    output_gjf.write(filename);
                    
                    System.out.printf("%.2f\t%s\n", newMolecule.getAngle(atomNumber1, atomNumber2, atomNumber3), filename);
                }
        }
}
