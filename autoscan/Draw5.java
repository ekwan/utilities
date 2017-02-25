import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw5
{
    // Make a series of gjf files where the specified bond distance has been set and frozen.
    public static void main(String[] args)
        {
            // read input file
            GJFfile input_gjf = new GJFfile("endo-b3lyp_d3bj-tz-pcm.gjf");
            Molecule molecule = input_gjf.molecule;

            // define the atom numbers to be frozen
            // the fragment attached to "toAtomNumber" will be moved
            int fromAtomNumber1 = 2;
            int toAtomNumber1   = 6;
            int fromAtomNumber2 = 4;
            int toAtomNumber2   = 5;

            molecule = molecule.addBond(fromAtomNumber1, toAtomNumber1);
            molecule = molecule.removeBond(fromAtomNumber2, toAtomNumber2);
            molecule = molecule.removeBond(6,18);

            // define the keywords for the input files that will be generated
            String keywords = "%mem=3GB\n%nprocshared=4\n#p geom=connect b3lyp/cc-pvtz empiricaldispersion=gd3bj scrf=(pcm,solvent=o-xylene) pop=none opt=modredundant";
            String tail = String.format("B %d %d F\nB %d %d F\n\n", fromAtomNumber1, toAtomNumber1, fromAtomNumber2, toAtomNumber2);
            
            for (double distance1 = 2.16; distance1 <= 2.28; distance1 += 0.01)
                {
                    for (double distance2 = 2.22; distance2 <= 2.36; distance2 += 0.01)
                        {
                            Molecule newMolecule = molecule.setDistance(fromAtomNumber1, toAtomNumber1, distance1);
                            
                            double bestTheta = 0.0;
                            double bestDelta = 100.0;
                            for (double theta = 70; theta < 130.0; theta += 5.0)
                                {
                                    Molecule testMolecule = newMolecule.setAngle(fromAtomNumber1, toAtomNumber1, toAtomNumber2, theta);
                                    Atom atom1 = testMolecule.contents.get(fromAtomNumber1-1);
                                    Atom atom2 = testMolecule.contents.get(toAtomNumber1-1);
                                    Atom atom3 = testMolecule.contents.get(fromAtomNumber2-1);
                                    Atom atom4 = testMolecule.contents.get(toAtomNumber2-1);
                                    double distance1actual = testMolecule.getDistance(atom1, atom2);
                                    double distance2actual = testMolecule.getDistance(atom3, atom4);
                                    if ( Math.abs(distance2-distance2actual) < bestDelta )
                                        {
                                            bestDelta = Math.abs(distance2-distance2actual);
                                            bestTheta = theta;
                                        }
                                    //System.out.printf("%6.2f %6.2f %6.2f\n", theta, distance1actual, distance2actual);
                               }
                             
                             for (double theta = bestTheta-5.0; theta <= bestTheta+5.0; theta += 0.001)
                                {
                                    Molecule testMolecule = newMolecule.setAngle(fromAtomNumber1, toAtomNumber1, toAtomNumber2, theta);
                                    Atom atom1 = testMolecule.contents.get(fromAtomNumber1-1);
                                    Atom atom2 = testMolecule.contents.get(toAtomNumber1-1);
                                    Atom atom3 = testMolecule.contents.get(fromAtomNumber2-1);
                                    Atom atom4 = testMolecule.contents.get(toAtomNumber2-1);
                                    double distance1actual = testMolecule.getDistance(atom1, atom2);
                                    double distance2actual = testMolecule.getDistance(atom3, atom4);
                                    if ( Math.abs(distance2-distance2actual) < bestDelta )
                                        {
                                            bestDelta = Math.abs(distance2-distance2actual);
                                            bestTheta = theta;
                                        }
                                    //System.out.printf("%6.2f %6.2f %6.2f\n", theta, distance1actual, distance2actual);
                               }
                            if ( bestDelta > 0.001 )
                                {
                                    System.out.printf("did not converge for %.2f, %.2f\n", distance1, distance2);
                                    continue;
                                }
                            
                            newMolecule = newMolecule.setAngle(fromAtomNumber1, toAtomNumber1, toAtomNumber2, bestTheta);
                            newMolecule = newMolecule.addBond(6,18);
                            newMolecule = newMolecule.removeBond(fromAtomNumber1, toAtomNumber1);

                            Atom atom1 = newMolecule.contents.get(fromAtomNumber1-1);
                            Atom atom2 = newMolecule.contents.get(toAtomNumber1-1);
                            Atom atom3 = newMolecule.contents.get(fromAtomNumber2-1);
                            Atom atom4 = newMolecule.contents.get(toAtomNumber2-1);
                            
                            String filename = String.format("output5/da_scan_%3.0f_%3.0f.gjf", distance1 * 100, distance2 * 100);
                            GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                            output_gjf.write(filename);
                            System.out.printf("%.4f\t%.4f\t%s\n", newMolecule.getDistance(atom1, atom2), newMolecule.getDistance(atom3, atom4), filename);
                        }
                }
        }
}
