import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class DrawTorsion
{
    public static void main(String[] args)
        {
            // read input file
            GJFfile input_gjf = new GJFfile("torsion.gjf");
            Molecule molecule = input_gjf.molecule;

            //molecule = molecule.addBond(fromAtomNumber1, toAtomNumber1);
            //molecule = molecule.addBond(fromAtomNumber2, toAtomNumber2);

            // define the keywords for the input files that will be generated
            String keywords = "%mem=4GB\n%nprocshared=4\n#p geom=connect b3lyp genecp empiricaldispersion=gd3bj scrf=(pcm,solvent=n-hexane) pop=none opt=(modredundant,maxcyc=50)";
            String tail = "D 11 12 13 14 F\nB 1 10 F\n\n@/n/jacobsen_lab/ekwan/basis/simple_basis.bas\n\n";
            
            for (double torsion = 136.0; torsion >= 74.0; torsion -= 2.0)
                {
                    for (double distance = 2.00; distance <= 5.10; distance += 0.10)
                        {
                            if ( distance < 4.0 && torsion >= 128.0 )
                                continue;
                            if ( distance > 4.0 && torsion <= 126.0 )
                                continue;
                            IndexTorsion indexTorsion = IndexTorsion.createIndexTorsion(11,12,13,14, molecule);
                            AtomTorsion atomTorsion = indexTorsion.getAtomTorsion(molecule);
                            Molecule newMolecule = molecule.setDihedral(atomTorsion, torsion);

                            Atom Tiatom = newMolecule.contents.get(1-1);
                            Atom Clatom = newMolecule.contents.get(10-1);
                            Atom C6atom = newMolecule.contents.get(14-1);

                            double a = newMolecule.getDistance(Tiatom, Clatom);
                            double b = newMolecule.getDistance(Tiatom, C6atom);
                            double c = distance;
                            //System.out.println(a + "   " + b + "   " + c);
                            double arg = Math.pow(a,2)+Math.pow(b,2)-Math.pow(c,2);
                            arg = arg / (2.0*a*b);
                            double theta = Math.toDegrees(Math.acos(arg));
                            if (Double.isNaN(theta))
                                continue;
                            //System.out.println(theta);

                            newMolecule = newMolecule.setAngle(14,1,10, theta);
                            
                            String filename = String.format("output_torsion/simple_torsion_scan_%03.0f_%03.0f.gjf", torsion, distance * 100);
                            GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                            output_gjf.write(filename);
                            
                            double actualTorsion = indexTorsion.getDihedralAngle(newMolecule);
                            Clatom = newMolecule.contents.get(10-1);
                            C6atom = newMolecule.contents.get(14-1);
                            double actualDistance = newMolecule.getDistance(Clatom,C6atom);
                            System.out.printf("%.0f,%.0f\t\t%.2f,%.2f\t%s\n", torsion, actualTorsion, distance, actualDistance, filename);
                        }
                }
        }
}
