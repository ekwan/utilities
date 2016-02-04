import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

public class Draw
{
    public static class Coordinate implements Immutable
    {
        public final double x;
        public final double y;

        public Coordinate(double x, double y)
        {
            this.x = round(x);
            this.y = round(y);
        }

        public static double round(double x)
        {
            return (double)Math.round(x*100.0)/100.0;
        }

        public double distanceFrom(Coordinate c)
        {
            return Math.sqrt( Math.pow(x-c.x, 2) + Math.pow(y-c.y,2) );
        }

        @Override
        public String toString()
        {
            return String.format("%.2f,%.2f", x, y);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(x,y);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            else if ( obj == this )
                return true;
            else if ( !(obj instanceof Coordinate) )
                return false;

            Coordinate c = (Coordinate)obj;
            return (x == c.x && y == c.y);
        }
    }

    public static void main(String[] args)
        {
            // define a grid
            // bond1 is the C5-Cl distance
            // bond2 is the C6-Cl distance
            List<Coordinate> desiredCoordinates = new ArrayList<>();
            for (double bond1 = 1.60; bond1 <= 3.05; bond1 += 0.05)
                {
                    for (double bond2 = 1.80; bond2 <= 3.80; bond2 += 0.05)
                    {
                        Coordinate c = new Coordinate(bond1, bond2);
                        desiredCoordinates.add(c);
                    }
                }

            // remove the points that have already been found
            // done.txt should contain entries like "1.4 3.5" corresponding to bond1 and bond2 points
            // lines starting with # are ignored
            OutputFileFormat alreadyDoneFile = new OutputFileFormat("done.txt") {};
            List<Coordinate> doneCoordinates = new ArrayList<>();
            for (List<String> line : alreadyDoneFile.fileContents)
                {
                    if (line.size() != 3 || line.get(0).startsWith("#"))
                        {
                            System.out.printf("Ignored line %s.\n", line.toString());
                            continue;
                        }
                    double bond1 = Double.parseDouble(line.get(0));
                    double bond2 = Double.parseDouble(line.get(1));
                    Coordinate c = new Coordinate(bond1, bond2);
                    doneCoordinates.add(c);
                    if (! desiredCoordinates.contains(c))
                        System.out.printf("Warning, coordinate %s is not in target range.\n", c.toString());
                    else
                        System.out.println(c);
                }

            // these points are still left
            List<Coordinate> coordinatesRequired = new ArrayList<>(desiredCoordinates);
            coordinatesRequired.removeAll(doneCoordinates);
            System.out.printf("%d coordinates still needed.\n", coordinatesRequired.size());

            // read template and ensure the atoms to be moved are connected
            GJFfile input_gjf = new GJFfile("template.gjf");
            Molecule molecule = input_gjf.molecule;

            int fromAtomNumber1 = 119;
            int toAtomNumber1 = 138;

            int fromAtomNumber2 = 138;
            int toAtomNumber2 = 67;

            if (! molecule.directlyConnected(fromAtomNumber1, toAtomNumber2))
                throw new IllegalArgumentException("bond1 atoms are not connected");
            if (! molecule.directlyConnected(fromAtomNumber2, toAtomNumber2))
                throw new IllegalArgumentException("bond2 atoms are not connected");

            // all files will be made with these keywords
            String keywords = "%mem=20GB\n%nprocshared=8\n#p geom=connect m062x/6-31g(d) scrf=(pcm,solvent=toluene) pop=none opt=modredundant";
            String tail = String.format("B %d %d F\nB %d %d F\n", fromAtomNumber1, toAtomNumber1, fromAtomNumber2, toAtomNumber2);
 
            // create up to the specified number of input files
            int filesMade = 0;
            for (Coordinate desiredPoint : coordinatesRequired)
                {
                    // for each required point, find out which points is closest
                    
                    // move the atoms
                    double bondDistance1 = desiredPoint.x;
                    double bondDistance2 = desiredPoint.y;
                    Molecule newMolecule = closestMolecule.setDistance(fromAtomNumber2, toAtomNumber2, bondDistance1);
                    newMolecule = newMolecule.setDistance(fromAtomNumber2, toAtomNumber2, bondDistance2); 

                    Atom atom1 = newMolecule.contents.get(fromAtomNumber1-1);
                    Atom atom2 = newMolecule.contents.get(toAtomNumber1-1);
                    double checkDistance1 = newMolecule.getDistance(atom1, atom2);

                    atom1 = newMolecule.contents.get(fromAtomNumber2-1);
                    atom2 = newMolecule.contents.get(toAtomNumber2-1);
                    double checkDistance2 = newMolecule.getDistance(atom1, atom2);

                    // make the file and write it to disk
                    String filename = String.format("output/scan1_%3d_%3d.gjf", bondDistance1 * 100, bondDistance2 * 100);
                    GaussianInputFile output_gjf = new GaussianInputFile(newMolecule, "title", keywords, tail);
                    output_gjf.write(filename);
                    
                    // update counter
                    filesMade++
                    System.out.printf("Made %d files.  %s    %.2f  %.2f\n", filesMade, filename, checkDistance1, checkDistance2);
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
