import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw
{
    public static void main(String[] args)
        {
            // define a grid
            // bond1 is the C5-Cl distance
            // bond2 is the C6-Cl distance
            List<Coordinate> desiredCoordinates = new ArrayList<>();
            for (double bond1 = 1.90; bond1 <= 3.05; bond1 += 0.05)
                {
                    for (double bond2 = 1.90; bond2 <= 3.80; bond2 += 0.05)
                    {
                        Coordinate c = new Coordinate(bond1, bond2);
                        desiredCoordinates.add(c);
                    }
                }

            // these are the atom numbers to be moved
            // the sceond atom number in each pair is the atom that will be moved
            int C5atomNumber = 17;
            int C6atomNumber = 18;
            int ClatomNumber = 14;

            // figure out which points have already been done
            String directory = "input/";
            List<Coordinate> doneCoordinates = new ArrayList<>();
            Map<Coordinate,Molecule> doneMap = new HashMap<>();
            for (File f : new File(directory).listFiles())
                {
                    String filename = f.getName();
                    if ( f.isFile() && filename.endsWith(".out") )
                        {
                            // ensure this file terminated normally
                            GaussianOutputFile out = null;
                            try { out = new GaussianOutputFile(directory + filename); }
                            catch (Exception e) { continue; }

                            int normalTerminationCount = 0;
                            for (List<String> fields : out.fileContents)
                                {
                                    if ( fields.size() > 5 && fields.get(0).equals("Normal") && fields.get(1).equals("termination") )
                                        normalTerminationCount++;
                                }
                            if ( normalTerminationCount == 0 )
                                {
                                    System.out.printf("%s failed and will be ignored.\n", filename);
                                    continue;
                                }

                            // determine which point this is
                            Molecule molecule = out.molecule;
                            if ( molecule.contents.size() != 50 )
                                {
                                    System.out.printf("File %s has the wrong number of atoms and will be ignored.\n");
                                    continue;
                                } 
                            Atom C5atom = molecule.contents.get(C5atomNumber-1);
                            Atom C6atom = molecule.contents.get(C6atomNumber-1);
                            Atom Clatom = molecule.contents.get(ClatomNumber-1);
                            double C5distance = molecule.getDistance(C5atom, Clatom);
                            double C6distance = molecule.getDistance(C6atom, Clatom);
                            double energy = molecule.energy;
                            Coordinate c = new Coordinate(C5distance, C6distance);
                            System.out.printf("Read file %s.  C5-Cl: %.2f  C6-Cl: %.2f  Energy: %.8f\n", filename, C5distance, C6distance, energy);
                            doneCoordinates.add(c);
                            doneMap.put(c,molecule);
                        }
                }

            // these points are still left
            List<Coordinate> coordinatesRequired = new ArrayList<>(desiredCoordinates);
            coordinatesRequired.removeAll(doneCoordinates);
            System.out.printf("%d coordinates still needed.\n", coordinatesRequired.size());

            // read template
            GaussianOutputFile r = new GaussianOutputFile("template.out");
            Molecule templateMolecule = r.molecule;
            
            // remove the bonds to avoid connectivity loops
            templateMolecule = templateMolecule.removeBond(1,14);
            templateMolecule = templateMolecule.removeBond(C5atomNumber, ClatomNumber);
            templateMolecule = templateMolecule.removeBond(C6atomNumber, ClatomNumber);

            // all files will be made with these keywords
            String keywords = "%mem=20GB\n%nprocshared=8\n#p geom=connect m062x/6-31g(d) scrf=(pcm,solvent=toluene) pop=none opt=modredundant";
            String tail = String.format("B %d %d F\nB %d %d F\n", C5atomNumber, ClatomNumber, C6atomNumber, ClatomNumber);
 
            // create up to the specified number of input files
            int desiredFiles = 5;
            int filesMade = 0;
            double distanceThreshold = 0.30;  // must use a template that is at least this close
            for (Coordinate desiredPoint : coordinatesRequired)
                {
                    // for each required point, find out which points is closest
                    Double minimumDistance = null;
                    Molecule closestMolecule = null;
                    Coordinate closestCoordinate = null;
                    for (Coordinate c : doneMap.keySet())
                        {
                            double thisDistance = desiredPoint.distanceFrom(c);
                            Molecule thisMolecule = doneMap.get(c);
                            if ( minimumDistance == null || thisDistance < minimumDistance )
                                {
                                    minimumDistance = thisDistance;
                                    closestMolecule = thisMolecule;
                                    closestCoordinate = c;
                                }
                        }

                    if ( minimumDistance > distanceThreshold )
                        continue;
                    System.out.printf("Closest reference point for target (%s) is (%s): distance = %.3f.\n",
                                      desiredPoint.toString(), closestCoordinate.toString(), minimumDistance);

                    // fix connectivity
                    // makes a molecule that has the same connectivity as the template
                    // but the same geometry as the closest molecule
                    Map<Atom,Atom> atomMap = new HashMap<>();
                    for (int i=0; i < closestMolecule.contents.size(); i++)
                        {
                            Atom fromAtom = templateMolecule.contents.get(i);
                            Atom toAtom   = closestMolecule.contents.get(i);
                            atomMap.put(fromAtom, toAtom);
                        }
                    Molecule thisTemplateMolecule = templateMolecule.moveAtoms(atomMap);

                    // these are the distances to set
                    double d1 = desiredPoint.x;
                    double d2 = desiredPoint.y;
                    Atom C5atom = thisTemplateMolecule.contents.get(C5atomNumber-1);
                    Atom C6atom = thisTemplateMolecule.contents.get(C6atomNumber-1);
                    double d3 = thisTemplateMolecule.getDistance(C5atom, C6atom);
                    //System.out.println(d1);
                    //System.out.println(d2);
                    //System.out.println(d3);
                    // set the first bond distance
                    Molecule newMolecule = thisTemplateMolecule.addBond(C5atomNumber, ClatomNumber);
                    newMolecule = newMolecule.setDistance(C5atomNumber, ClatomNumber, d1);
                    
                    // rotate the atom until the second bond distance is met
                    double desiredAngle = Math.pow(d1,2) + Math.pow(d3,2) - Math.pow(d2, 2);
                    desiredAngle = desiredAngle / (2.0 * d1 * d3);
                    desiredAngle = Math.toDegrees(Math.acos(desiredAngle));
                    //System.out.println(desiredAngle);
                    newMolecule = newMolecule.setAngle(C6atomNumber, C5atomNumber, ClatomNumber, desiredAngle);

                    // check we did it correctly
                    Atom atom1 = newMolecule.contents.get(C5atomNumber-1);
                    Atom atom2 = newMolecule.contents.get(ClatomNumber-1);
                    double checkDistance1 = newMolecule.getDistance(atom1, atom2);

                    atom1 = newMolecule.contents.get(C6atomNumber-1);
                    atom2 = newMolecule.contents.get(ClatomNumber-1);
                    double checkDistance2 = newMolecule.getDistance(atom1, atom2);

                    // make the file and write it to disk
                    String filename = String.format("output/scan1_%3.0f_%3.0f.gjf", d1 * 100, d2 * 100);
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                    output_gjf.write(filename);
                    
                    // update counter
                    filesMade++;
                    System.out.printf("Made %d files.  %s    %.2f  %.2f\n", filesMade, filename, checkDistance1, checkDistance2);
                    if ( filesMade >= desiredFiles )
                        {
                            System.out.println("Reached maximum number of requested files.");
                            break;
                        }
                }
/*
            // parse input file
            //*** atoms must be connected in template ***
            GJFfile input_gjf = new GJFfile("template6.gjf");
            Molecule molecule = input_gjf.molecule;

            // move atoms on grid
            int fromAtomNumber1 = 119;
            int toAtomNumber1 = 138;

            int fromAtomNumber2 = 138;
            int toAtomNumber2 = 67;

            //String keywords = "%chk=checkpoint.chk\n%mem=32GB\n%nprocshared=16\n#p geom=connect m062x/6-31g(d) scrf=(pcm,solvent=toluene) pop=none freq=noraman opt=modredundant";
            String keywords = "%chk=checkpoint.chk\n%mem=12GB\n%nprocshared=8\n#p geom=connect m062x/6-31g(d) scrf=(pcm,solvent=toluene) pop=none freq=noraman opt=modredundant";
            String tail = String.format("B %d %d F\nB %d %d F\n", fromAtomNumber1, toAtomNumber1, fromAtomNumber2, toAtomNumber2);
            
            for (double forming = 1.80; forming <= 3.00; forming += 0.1)
                {
                    Molecule newMolecule = molecule.setDistance(fromAtomNumber1,toAtomNumber1, forming);
                    for (double breaking = 1.80; breaking <= 3.40; breaking += 0.2)
                        {
                            if ( breaking - forming >= 0.50 )
                                continue;

                            newMolecule = newMolecule.setDistance(fromAtomNumber2, toAtomNumber2, breaking);

                            Atom atom1 = newMolecule.contents.get(fromAtomNumber1-1);
                            Atom atom2 = newMolecule.contents.get(toAtomNumber1-1);
                            System.out.printf("%.2f\t", newMolecule.getDistance(atom1, atom2)*100.0);

                            atom1 = newMolecule.contents.get(fromAtomNumber2-1);
                            atom2 = newMolecule.contents.get(toAtomNumber2-1);
                            System.out.printf("%.2f\n", newMolecule.getDistance(atom1, atom2)*100.0);

                            String filename = String.format("output/indocat_mannosyl_4E_grid_%.0f_%.0f.gjf", forming * 100, breaking * 100);
                            GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                            output_gjf.write(filename);
                        }
                }
*/
        }
}
