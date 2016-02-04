import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;

/**
 * This class reads Gaussian output files.  It reads the last energy and geometry that was
 * printed out.
 */
public class GaussianOutputFile extends OutputFileFormat
{
    /** The molecule. */
    public final Molecule molecule;

    /**
     * Reads a g09 output file from disk.  Reads the last geometry and potential energy, as well
     * as forces and normal modes, if they are available.  Verbose output should be called with
     * #p so that the masses can be read.
     * @param filename the g09 output file to read from
     */
    public GaussianOutputFile(String filename)
    {
        super(filename);

        // check for Link1
        int link1count = 0;
        for (int i=0; i < fileContents.size(); i++)
            {
                List<String> fields = fileContents.get(i);
                if ( fields.size() > 4 && fields.get(0).equals("Entering") && fields.get(1).equals("Link") && fields.get(2).equals("1") )
                    {
                        link1count++;
                        if (link1count > 1)
                            {
                                System.out.printf("Warning, %s has more than one Link1 directive.\n", filename);
                                break;
                            }
                    }
            }

        // read atom symbols
        List<String> atomSymbols = new ArrayList<>();
        for (int i=0; i < fileContents.size(); i++)
            {
                List<String> fields = fileContents.get(i);
                if ( fields.size() == 6 && fields.get(0).equals("Charge") && fields.get(3).equals("Multiplicity") )
                    {
                        while (true)
                            {
                                i++;
                                fields = fileContents.get(i);
                                if ( fields.size() < 4 || fields.size() > 5 )
                                    break;
                                atomSymbols.add(fields.get(0));
                            }
                        break;
                    }
            }
        //System.out.println(atomSymbols.size() + " / " + atomSymbols);

        // read last geometry
        List<Vector3D> positions = new ArrayList<>(atomSymbols.size());
        for (int i=0; i < fileContents.size(); i++)
            {
                List<String> fields = fileContents.get(i);
                if ( fields.size() == 2 && fields.get(0).equals("Standard") && fields.get(1).equals("orientation:") )
                    {
                        i += 5;
                        positions.clear();
                        for (int j=i; j < i+atomSymbols.size(); j++)
                            {
                                fields = fileContents.get(j);
                                int numberOfFields = fields.size();
                                if ( numberOfFields < 5 || numberOfFields > 6 )
                                    throw new IllegalArgumentException("unexpected number of geometry fields:\n" + fields.toString());
                                double x = Double.valueOf(fields.get(numberOfFields-3));
                                double y = Double.valueOf(fields.get(numberOfFields-2));
                                double z = Double.valueOf(fields.get(numberOfFields-1));
                                positions.add(new Vector3D(x,y,z));
                            }
                        i += atomSymbols.size();
                    }
            }
        if ( positions.size() == 0 || positions.size() != atomSymbols.size() )
            throw new IllegalArgumentException("error reading positions!");
        //for (Vector3D v : positions)
        //    System.out.println(v);

        // construct atoms
        List<Atom> contents = new ArrayList<>();
        for (int i=0; i < atomSymbols.size(); i++)
            {
                String symbol = atomSymbols.get(i);
                double mass = 0.0;
                Vector3D position = positions.get(i);
                Atom atom = new Atom(symbol, position, 0);
                contents.add(atom);
                //System.out.println(atom);
            }

        // read connectivity
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        for (Atom a : contents)
            connectivity.addVertex(a);
        for (int i=0; i < fileContents.size(); i++)
            {
                List<String> fields = fileContents.get(i);
                if ( fields.size() == 4 && fields.get(1).equals("Initial") && fields.get(2).equals("Parameters") )
                    {
                        i+=5;
                        while (true)
                            {
                                fields = fileContents.get(i);
                                if ( fields.size() >= 4 )
                                    {
                                        String firstLetter = fields.get(1).substring(0,1);
                                        if ( firstLetter.equals("R") )
                                            {
                                                String[] fields2 = fields.get(2).split(",");
                                                int atomNumber1 = Integer.parseInt(fields2[0].replaceAll("[^\\d.]", ""));
                                                int atomNumber2 = Integer.parseInt(fields2[1].replaceAll("[^\\d.]", ""));
                                                //System.out.print(fields);
                                                //System.out.printf("   %d, %d\n", atomNumber1, atomNumber2);
                                                Atom fromAtom = contents.get(atomNumber1-1);
                                                Atom toAtom = contents.get(atomNumber2-1);
                                                double bondOrder = 1.0;                                                
                                                DefaultWeightedEdge thisEdge = connectivity.addEdge(fromAtom, toAtom);
                                                connectivity.setEdgeWeight(thisEdge, bondOrder);
                                            }
                                    }
                                else if ( fields.size() == 1 )
                                    break;
                                i++;
                            }
                    }
            }

        // read any virtual bonds
        /*
        for (List<String> fields : fileContents)
            {
                if ( fields.size() > 5 && fields.get(1).equals("virtual") && fields.get(2).equals("bond") )
                    {
                        int atomNumber1 = Integer.parseInt(fields.get(5).replaceAll("[^\\d.]", ""));
                        int atomNumber2 = Integer.parseInt(fields.get(7).replaceAll("[^\\d.]", ""));
                        //System.out.print(fields);
                        //System.out.printf("   %d, %d\n", atomNumber1, atomNumber2);
                        Atom fromAtom = contents.get(atomNumber1-1);
                        Atom toAtom = contents.get(atomNumber2-1);
                        double bondOrder = 1.0;
                        DefaultWeightedEdge thisEdge = connectivity.addEdge(fromAtom, toAtom);
                        System.out.println(connectivity.getAllEdges(fromAtom,toAtom).size());
                        if (thisEdge == null)
                            continue;
                        connectivity.setEdgeWeight(thisEdge, bondOrder);
                    }
            }
        */

        // read last potential energy
        double potentialEnergy = 0.0;
        for (List<String> fields : fileContents)
            {
                if ( fields.size() > 5 && fields.get(0).equals("SCF") && fields.get(1).equals("Done:") )
                    potentialEnergy = Double.parseDouble(fields.get(4));
            }
        //System.out.println(potentialEnergy);

        // create object
        this.molecule = new Molecule("title", contents, connectivity, potentialEnergy);
    }

    @Override 
    public String toString()
    {
        return molecule.toString();
    }

    /** For testing. */
    public static void main(String[] args)
    {
        GaussianOutputFile r = new GaussianOutputFile("template.out");
        String keywords = "# geom=connect";
        String tail = " ";
        Molecule molecule = r.molecule;
        molecule = molecule.removeBond(1,14);
        molecule = molecule.removeBond(17,14);
        molecule = molecule.removeBond(18,14);
        GaussianInputFile output_gjf = new GaussianInputFile(molecule, "title", keywords, tail);
        output_gjf.write("template.gjf");
        //System.out.println(r);
        //System.out.println(r.molecule.energy);
    }
}
