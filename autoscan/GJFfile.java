import java.io.*;
import java.util.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;

/**
 * Represents a Gaussian input file for an octahedral complex.  Stores geometry,
 * connectivity, ligand atom labels, requested memory, requested processors,
 * method, and filename of basis set.  These will be specified as @LIGAND_NAME@ATOM_NUMBER, @MEM@X(GB),
 * @NPROCSHARED@X, @method@method, and @BASIS@basis.
 */
public class GJFfile extends OutputFileFormat implements Immutable
{
    /** The Molecule representation of this file. */
    public final Molecule molecule;

    /** Atoms at which titanium connects to each ligand. i=1,2,...*/
    public int O1Number;
    public int O2Number;
    public int N3Number;
    public int Cl1Number;
    public int Su2Number;
    public int Ol3Number;

    /** Memory, processors, method, basis. */
    public int mem;
    public int nprocshared;
    public String method;
    public String basis;
    
    /**
     * Reads the geometry and connectivity.
     * @param filename the location of the gjf file
     */
    public GJFfile(String filename)
    {
	    super(filename);

        // read geometry
	    String name = "";
        List<Atom> contents = new ArrayList<>();
	    SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        int blanks = 0;
        boolean lastBlank = false;
        boolean inGeometryBlock = false;
        for (List<String> line : fileContents)
            {
                // keep track of how many blanks we have seen
                if ( line.size() == 1 && line.get(0).length() == 0 )
                    {
                        if ( lastBlank == false )
                            {
                                blanks++;
                                lastBlank = true;
                            }
                        continue;
                    }
                else
                    lastBlank = false;
                
                // read the metadata 
                if ( blanks == 1 )
                    {   
                        for (String s : line)
                        {  
                            String[] fields = s.split("@");
                            if ( fields.length != 3 )
                                continue;
                            String identifier = fields[1].toLowerCase();
                            String value = fields[2];
                            //System.out.println(s);
                            //System.out.println(identifier + " : " + value);

                            if ( identifier.equals("o1") )
                                O1Number = Integer.parseInt(value);
                            else if ( identifier.equals("o2") )
                                O2Number = Integer.parseInt(value);
                            else if ( identifier.equals("n3") )
                                N3Number = Integer.parseInt(value);
                            else if ( identifier.equals("cl1") )
                                Cl1Number = Integer.parseInt(value);
                            else if ( identifier.equals("su2") )
                                Su2Number = Integer.parseInt(value);
                            else if ( identifier.equals("ol3") )
                                Ol3Number = Integer.parseInt(value);
                            else if ( identifier.equals("mem") )
                                mem = Integer.parseInt(value);
                            else if ( identifier.equals("nprocshared") )
                                nprocshared = Integer.parseInt(value);
                            else if ( identifier.equals("method") )
                                method = value;
                            else if ( identifier.equals("basis") )
                                basis = value;
                            else
                                System.out.println("unrecognized entry: " + s);
                        }

                        continue;
                    }
                else if ( blanks != 2 )
                    continue;

                // deal with the charge and multiplicity card (by ignoring it)
                if ( line.size() == 2 && inGeometryBlock == false )
                    {
                        inGeometryBlock = true;
                        continue;
                    }

                if ( line.size() != 4 && inGeometryBlock == false )
                    throw new IllegalArgumentException("unexpected text in geometry block in " + filename + ":\n" + line.toString());

                // create atom
                // tinker atom types will be nonsense, of course
                Atom newAtom = new Atom(line.get(0), new Vector3D(Double.parseDouble(line.get(1)), Double.parseDouble(line.get(2)), Double.parseDouble(line.get(3))), 1);
                contents.add(newAtom);
                connectivity.addVertex(newAtom);
            }
        
        // read connectivity
        blanks = 0;
        lastBlank = false;
        for (List<String> line : fileContents)
            {
                // read the fourth block of text
                if ( line.size() == 1 && line.get(0).length() == 0 )
                    {
                        if ( lastBlank == false )
                            {
                                blanks++;
                                lastBlank = true;
                            }
                        continue;
                    }
                else
                    lastBlank = false;
               
               // only read connectivity lines
                if ( blanks != 3 )
                    continue;

                Atom fromAtom = contents.get(Integer.parseInt(line.get(0))-1);
                for (int i=1; i < line.size(); i+=2)
                    {
                        int toAtomIndex = Integer.parseInt(line.get(i))-1;
                        Atom toAtom = contents.get(toAtomIndex);
                        double bondOrder = Double.parseDouble(line.get(i+1));
                        DefaultWeightedEdge thisEdge = connectivity.addEdge(fromAtom, toAtom);
                        connectivity.setEdgeWeight(thisEdge, bondOrder);
                    }

            }

        // create the molecule
	    molecule = new Molecule(name, contents, connectivity, 0.0);
    }

    /** Writes out the GJF file to disk.  Metadata not reproduced.  Charge and multiplicity assumed to be 1. */
    public void write(String filename)
    {
        InputFileFormat.writeStringToDisk(this.toString(), filename);
    }

    public void write(Molecule newMolecule, String filename)
    {
        InputFileFormat.writeStringToDisk(this.toString(newMolecule), filename);
    }

    public String toString()
    {
        return this.toString(this.molecule);
    }

    public String toString(Molecule molecule)
    {
        // write header
        String returnString = String.format("%%mem=%dGB\n%%nprocshared=%d\n#p geom=connect %s/genecp empiricaldispersion=gd3bj opt\n", 
                                            mem, nprocshared, method);
        returnString += "\ntitle\n\n0 1\n";
        
        // write geometry
        returnString = returnString + molecule.getGJFstring() + "\n";
        
        // initialize some variables
        Atom fromAtom = null;
        Atom toAtom = null;
        Integer fromAtomNumber = 0;
        Integer toAtomNumber = 0;
        Double bondOrder = 0.0;
        DefaultWeightedEdge thisEdge = null;

        // read connectivity data into parallel arrays
        ArrayList<Integer> fromAtoms = new ArrayList<Integer>();
        ArrayList<Integer> toAtoms = new ArrayList<Integer>();
        ArrayList<Double> bondOrders = new ArrayList<Double>();
        ArrayList<Boolean> visited = new ArrayList<Boolean>();
        for (DefaultWeightedEdge e : molecule.connectivity.edgeSet())
            {
                fromAtom = molecule.connectivity.getEdgeSource(e);
                fromAtomNumber = molecule.contents.indexOf(fromAtom) + 1;
                toAtom = molecule.connectivity.getEdgeTarget(e);
                toAtomNumber = molecule.contents.indexOf(toAtom) + 1;
                bondOrder = molecule.connectivity.getEdgeWeight(e);

                fromAtoms.add(fromAtomNumber);
                toAtoms.add(toAtomNumber);
                bondOrders.add(bondOrder);
                visited.add(false);
            }

        // write connectivity data
        for (int i=0; i < molecule.contents.size(); i++)
            {
                returnString = returnString + (i+1) + " ";
                for (int j=0; j < fromAtoms.size(); j++)
                    {
                        if (fromAtoms.get(j) == i+1 && visited.get(j) == false)
                            {
                                returnString = returnString + toAtoms.get(j) + " " + String.format("%.1f ", bondOrders.get(j));
                                visited.set(j, true);
                            }
                        if (toAtoms.get(j) == i+1 && visited.get(j) == false)
                            {
                                returnString = returnString + fromAtoms.get(j) + " " + String.format("%.1f ", bondOrders.get(j));
                                visited.set(j, true);
                            }
                    }
                returnString = returnString + "\n";
            }

        // write footer
        returnString += String.format("\n@%s\n\n", basis);
        return returnString;
    }

    public void printDebugString()
    {
        System.out.printf("O1:  %d\n", O1Number);
        System.out.printf("O2:  %d\n", O2Number);
        System.out.printf("N3:  %d\n", N3Number);
        System.out.printf("Cl1: %d\n", Cl1Number);
        System.out.printf("Su2: %d\n", Su2Number);
        System.out.printf("Ol3: %d\n", Ol3Number);
        System.out.printf("mem: %d\n", mem);
        System.out.printf("prc: %d\n", nprocshared);
        System.out.printf("mtd: %s\n", method);
        System.out.printf("bas: %s\n", basis);
    }

    /** for testing */
    public static void main(String args[])
    {
        GJFfile g = new GJFfile("templates/merSR33.gjf");
        //g.printDebugString();
        System.out.println(g);
        g.write("test.gjf");
    }
}
