import java.io.*;
import java.util.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;

/** Represent a Gaussian input file for debugging purposes. */
public class GaussianInputFile extends InputFileFormat
{
    public GaussianInputFile(Molecule molecule)
    {
        super(getGaussianString(molecule,molecule.name,"",""));
    }

    public GaussianInputFile(Molecule molecule, String name)
    {
        super(getGaussianString(molecule,name,"",""));
    }

    public GaussianInputFile(Molecule molecule, String name, String keywords)
    {
        super(getGaussianString(molecule,name,keywords,""));
    }

    public GaussianInputFile(Molecule molecule, String name, String keywords, String tail)
    {
        super(getGaussianString(molecule,name,keywords,tail));
    }

    public GaussianInputFile(Molecule molecule, Set<Atom> includeAtoms)
    {
        super(getGaussianString(molecule,includeAtoms));
    }

    private static String getGaussianString(Molecule molecule, String name, String keywords, String tail)
    {
        //String returnString = "%chk=checkpoint.chk\n%mem=32GB\n%nprocshared=16\n#p geom=connect m062x/6-31g(d) scrf=(pcm,solvent=toluene) pop=none freq=noraman opt=modredundant " + keywords + "\n";
        String returnString = "";
        if ( keywords.length() == 0 )
            returnString = "#\n";
        else
            returnString = keywords + "\n";
        returnString += "\n" + name + "\n\n0 1\n";
        returnString = returnString + molecule.toOmnisolString() + "\n";
        
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
        returnString = returnString + "\n" + tail + "\n\n";
        return returnString;
    }

    /**
     * Writes a gjf file that only includes part of the molecule.
     * @param molecule the molecule whose geometry we want to write partially
     * @param includeAtoms the atoms to include the result
     * @return a String containing a gjf that only includes includeAtoms
     */
    private static String getGaussianString(Molecule molecule, Set<Atom> includeAtoms)
    {
        String returnString = "%mem=1GB\n%nprocshared=12\n#p geom=connect\n";
        returnString += "\nname\n\n0 1\n";

        List<Atom> newContents = new LinkedList<>();
        Set<Atom> excludedAtoms = new HashSet<>();
        for (Atom a : molecule.contents)
            {
                if ( includeAtoms.contains(a) )
                    newContents.add(a);
                else
                    excludedAtoms.add(a);
            }
        Collections.sort(newContents);

        Molecule tempMolecule = molecule.moveAtoms(new HashMap<Atom,Atom>());
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = tempMolecule.connectivity;
        for (Atom a : excludedAtoms)
            connectivity.removeVertex(a);

        Molecule tempMolecule2 = new Molecule(molecule.name, newContents, connectivity, 0.0);

        return getGaussianString(tempMolecule2, molecule.name, "", "");
    }
}
