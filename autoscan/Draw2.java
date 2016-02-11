//Draw4 by Fritz Seidl 2016-2-10, adapted from Draw.java and SimpleDraw.java by Eugene
//Usage: compile2-forPC Draw4.java.
//Check all the strings prior to running. Place template in working directory. Make sure your output folder exists and is empty.
//Suggested 4 processors, normal queue, 360 mins (6 hr)
import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw2
{
    public static void main(String[] args)
        {
			//Make template molecule from a .gjf in the working directory.
            GJFfile input_gjf = new GJFfile("template3.gjf");
            Molecule molecule = input_gjf.molecule;
			
			//Assign atom numbers.
			int C5atomNumber = 6; //distal C
            int C6atomNumber = 5; //proximal C
            int ClatomNumber = 18;
			int BratomNumber = 17;
			Atom C5atom = molecule.contents.get(C5atomNumber-1);
			Atom C6atom = molecule.contents.get(C6atomNumber-1);
			double distCC = molecule.getDistance(C5atom, C6atom);
			molecule = molecule.addBond(C5atomNumber, ClatomNumber);
						
			//Make header and tail. 
			String keywords = "%mem=4GB\n%nprocshared=4\n#p geom=connect b3lyp/6-31g(d) scrf=(pcm,solvent=toluene) pop=none opt=modredundant freq=noraman";
            String tail = String.format("B %d %d F\nB %d %d F\n", C5atomNumber, ClatomNumber, C6atomNumber, ClatomNumber);

			//bond 1 is the long, distal C5-Cl bond, which is not formed during the reaction (leg C in triangle.)
			//bond 2 is the shorter, proximal C6-Cl bond, which is being formed.
			//These two for-loops are centered around the bond lengths in full_chlorination1_TS.out:
			//(bond 1 was 3.9A, and bond 2 was 3.05A.)
            for (double bond1 = 3.4; bond1 <= 4.4; bond1 += 0.10)
                {
                    for (double bond2 = 2.5; bond2 <= 3.5; bond2 += 0.10)
                    {
						//initialize thisMolecule
						Molecule thisMolecule = molecule;
					
						thisMolecule = thisMolecule.setDistance(C5atomNumber,ClatomNumber,bond1);
						
						// rotate the atom until the bond2 distance is met
						double desiredAngle = Math.pow(bond1,2) + Math.pow(distCC,2) - Math.pow(bond2, 2);
						desiredAngle = desiredAngle / (2.0 * bond1 * distCC);
						desiredAngle = Math.toDegrees(Math.acos(desiredAngle));
						//System.out.printf("bond1 %.2f   CCbond %.2f   bond2 %.2f   delta %.4f\n", bond1, distCC, bond2, bond1+distCC-bond2);
						// double arg1 = Math.pow(bond1,2) + Math.pow(distCC,2) - Math.pow(bond2, 2);
						// double arg2 = 2.0 * bond1 * distCC;
						// double arg3 = arg1/arg2;
						// System.out.printf("%.4f  %.4f  %.4f\n", arg1, arg2, arg3);
						// System.out.printf("%.4f  %.4f  %.4f  %.4f\n", bond1, distCC, bond2, bond1-distCC-bond2);

						if ( Double.isNaN(desiredAngle) )
							{
								continue;
							}
						//System.out.println(desiredAngle);
						thisMolecule = thisMolecule.setAngle(C6atomNumber, C5atomNumber, ClatomNumber, desiredAngle);						
						
						// check the desired distances have been met
						C5atom = thisMolecule.contents.get(C5atomNumber-1);
						C6atom = thisMolecule.contents.get(C6atomNumber-1);
						Atom Clatom = thisMolecule.contents.get(ClatomNumber-1);
						double dist1 = thisMolecule.getDistance(C5atom, Clatom);
						double dist2 = thisMolecule.getDistance(C6atom, Clatom);

						//Write this file. Note: first number is proximal C-Cl, second number is distal C-Cl. 
						String filename = String.format("output2/b3lyp-631gd-scan_%3.0f_%3.0f.gjf", bond2 * 100, bond1 * 100);
						System.out.printf("%s %.2f %.2f\n", filename, dist2, dist1);
						GaussianInputFile output_gjf = new GaussianInputFile(thisMolecule, "title", keywords, tail);
						output_gjf.write(filename);
                    }
                }
		}
}