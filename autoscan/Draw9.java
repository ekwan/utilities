import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw9
{
    public static void main(String[] args)
        {
            // define a grid
            // bond1 is the C5-Cl distance
            // bond2 is the C6-Cl distance
            List<Coordinate> desiredCoordinates = new ArrayList<>();
            for (double bond1 = 1.30; bond1 <= 3.20; bond1 += 0.10)
                {
                    for (double bond2 = 1.30; bond2 <= 3.40; bond2 += 0.10)
                    {
                        if ( bond1 + bond2 < 3.30 ) // || bond1 + bond2 > 4.60 )
                            continue;
                        Coordinate c = new Coordinate(bond1, bond2);
                        desiredCoordinates.add(c);
                    }
                }

            // these are the atom numbers to be moved
            // the sceond atom number in each pair is the atom that will be moved

            // figure out which points have already been done
            String directory = "/Users/ekwan/research/clay/grid/good";
            List<Coordinate> doneCoordinates = new ArrayList<>();
            Map<Coordinate,Molecule> doneMap = new HashMap<>();
            for (File f : new File(directory).listFiles())
                {
                    String filename = f.getName();
                    if ( f.isFile() && filename.endsWith(".out") )
                        {
                            // ensure this file terminated normally
                            GaussianOutputFile out = null;
                            try { out = new GaussianOutputFile(directory + "/" + filename); }
                            catch (Exception e) { e.printStackTrace(); continue; }

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
                            if ( molecule.contents.size() != 47 )
                                {
                                    System.out.printf("File %s has the wrong number of atoms and will be ignored.\n");
                                    continue;
                                } 
                            Atom atomC = molecule.contents.get(4-1);
                            Atom atomNuc = molecule.contents.get(1-1);
                            Atom atomLG = molecule.contents.get(13-1);
                            double formingDistance = molecule.getDistance(atomC, atomNuc);
                            double breakingDistance = molecule.getDistance(atomC, atomLG);
                            double energy = molecule.energy;
                            Coordinate c = new Coordinate(formingDistance, breakingDistance);
                            System.out.printf("Read file %s.  forming: %.2f  breaking: %.2f  Energy: %.8f\n", filename, formingDistance, breakingDistance, energy);
                            doneCoordinates.add(c);
                            doneMap.put(c,molecule);
                        }
                }

            // these points are still left
            List<Coordinate> coordinatesRequired = new ArrayList<>(desiredCoordinates);
            coordinatesRequired.removeAll(doneCoordinates);
            System.out.printf("%d coordinates still needed.\n", coordinatesRequired.size());


            // read template
            GJFfile input_gjf = new GJFfile("clay-stripped.gjf");
            Molecule templateMolecule = input_gjf.molecule;
            templateMolecule = templateMolecule.addBond(4, 1);
            templateMolecule = templateMolecule.addBond(4, 13);
            //GaussianOutputFile r = new GaussianOutputFile("template.out");
            //Molecule templateMolecule = r.molecule;
            
            // remove the bonds to avoid connectivity loops
            //templateMolecule = templateMolecule.removeBond(1,14);
            //templateMolecule = templateMolecule.removeBond(C5atomNumber, ClatomNumber);
            //templateMolecule = templateMolecule.removeBond(C6atomNumber, ClatomNumber);

            // all files will be made with these keywords
            String keywords = "%mem=8GB\n%nprocshared=8\n#p empiricaldispersion=gd3bj b3lyp/6-31g(d) scrf=(pcm,solvent=tetrahydrofuran) pop=none opt=(calcfc,modredundant) freq=noraman";
            String tail = "B 4 1 F\nB 4 13 F\n";
 
            // create up to the specified number of input files
            int desiredFiles = 500;
            int filesMade = 0;
            double distanceThreshold = 0.25;  // must use a template that is at least this close
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

                    if ( minimumDistance > distanceThreshold ) {
                        System.out.printf("Skipping target %s because nothing is close enough.\n", desiredPoint.toString());
                        continue;
                    }

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
                    
                    Molecule newMolecule = thisTemplateMolecule.setDistance(4,1,d1);
                    newMolecule = newMolecule.setDistance(4,13,d2);
                    double checkDistance1 = newMolecule.getDistance(newMolecule.contents.get(4-1),newMolecule.contents.get(1-1));
                    double checkDistance2 = newMolecule.getDistance(newMolecule.contents.get(4-1),newMolecule.contents.get(13-1));

                    // make the file and write it to disk
                    String filename = String.format("gjf/clay-stripped-surface-%3.0f_%3.0f-b3lyp_d3bj-631gd-thf_pcm.gjf", d1 * 100, d2 * 100);
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail, false, -1, 1);
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
