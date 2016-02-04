import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import Jama.*;

/**
 * Represents a molecule.  This class is effectively immutable and serializable.
 * The connectivity graph is not exposed because it is potentially mutable.
 */
public class Molecule implements Immutable, Serializable
{
    public static final long serialVersionUID = 1L;

    public final String name;
    public final List<Atom> contents;
    protected final SimpleWeightedGraph<Atom, DefaultWeightedEdge> connectivity;
    public final double energy;

    /**
     * Factory method to create a molecule given a map of old atoms to new atoms.  Should be used
     * to move atoms.
     * @param atomMap a map from old atoms to new atoms (does not have to include all atoms)
     */
    public Molecule moveAtoms(Map<Atom,Atom> atomMap)
    {
        // copy the list of vertices
        List<Atom> newContents = new LinkedList<Atom>();
        for (Atom a : contents)
            {
                if ( atomMap.containsKey(a) )
                    newContents.add(atomMap.get(a));
                else
                    newContents.add(a);
            }

        // populate a new connectivity graph
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
        for (Atom newAtom : newContents)
            newConnectivity.addVertex(newAtom);
        for (DefaultWeightedEdge e : connectivity.edgeSet())
            {
                // get old edge data
                Double bondOrder = connectivity.getEdgeWeight(e);
                Atom fromAtom    = connectivity.getEdgeSource(e);
                Atom toAtom      = connectivity.getEdgeTarget(e);

                // replace any changes
                if ( atomMap.containsKey(fromAtom) )
                    fromAtom = atomMap.get(fromAtom);
                if ( atomMap.containsKey(toAtom)   )
                    toAtom   = atomMap.get(toAtom);

                // create new edge
                if (! newContents.contains(fromAtom) || ! newContents.contains(toAtom))
                    System.out.println("FromAtom: " + getAtomString(fromAtom) +  " ToAtom: " + getAtomString(toAtom));
                DefaultWeightedEdge newEdge = newConnectivity.addEdge(fromAtom,toAtom);
                newConnectivity.setEdgeWeight(newEdge, bondOrder);
            }

        // return result
        return new Molecule(name, newContents, newConnectivity, energy);
    }

    /**
     * Factory method to create a new molecule by transforming this one.
     * The rotation is applied before the translation.  
     * @param rot a three-dimensional rotation that we apply to the Atoms in this to get the Atoms in the output
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output
     * @return this rotated by rot, plus shift
     */
    public Molecule transform(Rotation rot, Vector3D shift)
    {
        // copy the list of vertices
        List<Atom> newContents = new LinkedList<Atom>();
        for (Atom a : contents) newContents.add(a.transform(rot,shift));

        // populate a new connectivity graph
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
        for (Atom newAtom : newContents) newConnectivity.addVertex(newAtom);
        for (DefaultWeightedEdge e : connectivity.edgeSet())
            {
                // get old edge data
                Double bondOrder = connectivity.getEdgeWeight(e);
                Atom fromAtom = connectivity.getEdgeSource(e);
                Atom toAtom = connectivity.getEdgeTarget(e);

                // create new edge
                DefaultWeightedEdge newEdge = newConnectivity.addEdge(newContents.get(contents.indexOf(fromAtom)),newContents.get(contents.indexOf(toAtom)));
                newConnectivity.setEdgeWeight(newEdge, bondOrder);
            }

        // return result
        return new Molecule(name, newContents, newConnectivity, energy);
    }

    /**
     * Factory method to create a new molecule by shifting this one.
     * Uses Rotation.IDENTITY as a parameter in the transform method
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output
     * @return this plus shift
     */
    public Molecule shift(Vector3D shift)
    {
        return transform(Rotation.IDENTITY, shift);
    }

    /**
     * Factory method to create a new molecule by shifting this one to have its center at the origin.
     * @return this shifted by minus the barycenter of this 
     */
    public Molecule normalize()
    {
        Vector3D barycenter = this.getCentroid();
        return shift(barycenter.negate());
    }

    /**
     * Returns the centroid of this molecule.
     * @return the centroid vector
     */
    public Vector3D getCentroid()
    {
        Vector3D barycenter = Vector3D.ZERO;
        for (Atom a : contents) barycenter = barycenter.add(a.position);
        barycenter = barycenter.scalarMultiply(1.0/contents.size());
        return barycenter;
    }

    /**
     * Create a new Molecule that is a transformed copy of this one, with the Molecule
     * other added to it.  There will be no bonds between the transformed copy and other.
     * The rotation is applied before the translation.  
     * @param rot a three-dimensional rotation that we apply to the Atoms in this to get the Atoms in the output
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output. 
     * @param other a Molecule that we throw in at the end
     * @return two Molecules combined into a single Molecule; null if they are too close
     */
    @SuppressWarnings("unchecked")
    public Molecule transformAndCombine(Rotation rot, Vector3D shift, Molecule other)
    {
        Element O = Element.getElement("O");
        Element N = Element.getElement("N");
        Element H = Element.getElement("H");
        
        Molecule output = transform(rot,shift);
        //if ( output.tooClose(other) || (!output.isInteresting(other)) )
        if ( !output.isInteresting(other) )
            return null;
        
        // custom tooClose
        for (Atom a : output.contents)
            {
                for (Atom b : other.contents)
                    {
                        double minimumDistance = 2.5;
                        // a   b  a   b
                        // H...N  H...O
                        // N...H  O...H
                        if ( ( a.element.equals(H) && b.element.equals(N) && other.getAdjacentAtoms(b).size() == 2 ) ||
                             ( a.element.equals(N) && b.element.equals(H) && output.getAdjacentAtoms(a).size() == 2 ) ||
                             ( a.element.equals(H) && b.element.equals(O) ) ||
                             ( a.element.equals(O) && b.element.equals(H) ) )
                            minimumDistance = 2.0;
                        double distance = Vector3D.distance(a.position, b.position);
                        //if ( distance < 3.0 )
                        //    System.out.printf("%.2f\n", distance);
                        if ( distance < minimumDistance )
                            return null;
                    }
            }

        SimpleWeightedGraph<Atom, DefaultWeightedEdge> c = (SimpleWeightedGraph<Atom, DefaultWeightedEdge>) output.connectivity.clone();
        Graphs.addGraph(c,other.connectivity);
        List<Atom> newContents = new LinkedList<>();
        for (Atom a : output.contents) newContents.add(a);
        for (Atom a : other.contents) newContents.add(a);

        return new Molecule(output.name + "/" + other.name, newContents, c, 0.0);
    }

    /**
     *  called by TemplateFile to create Molecule.
     */
    @SuppressWarnings("unchecked")
    public Molecule(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity, double energy)
    {
	    this.name = name;
	    this.contents = ImmutableList.copyOf(contents);

        // will throw an unchecked cast exception
        this.connectivity = (SimpleWeightedGraph<Atom,DefaultWeightedEdge>)connectivity.clone();

        this.energy = energy;
    }

    /**
     * Determines if this Atom is contained in this Molecule.
     * @param atom the atom that is supposed to be in this molecule
     * @return true if the atom is in the molecule
     */
    public boolean containsAtom(Atom atom)
    {
        if ( contents.contains(atom) && connectivity.vertexSet().contains(atom) )
            return true;
        return false;
    }

    /**
     * Returns the number of the atom.  Answers are given in 1,2,...n where n is
     * the number of atoms in the molecule.
     * @param atom an Atom in this Molecule
     * @return the requested atom number (-1 if it is not in the molecule)
     */
    public int getAtomNumber(Atom atom)
    {
        if ( containsAtom(atom) )
            return contents.indexOf(atom) + 1;
        else
            return -1;
    }

    /**
     * Returns the element and atom number of an atom.  e.g., C5.
     * Throws an IllegalArgumentException if the Atom is not in the molecule.
     * @param atom the atom whose string is desired
     * @return the description of the atom
     */
    public String getAtomString(Atom atom)
    {
        if ( atom == null)
            throw new NullPointerException("atom cannot be null");

        int atomNumber = getAtomNumber(atom);
        
        if ( atomNumber == -1 )
            throw new IllegalArgumentException("cannot print out atom because it is not in the molecule!\n" + atom.toString() + " TAT: " + atom.tinkerAtomType);
        return atom.element.symbol + atomNumber;
    }

    public String getClosestAtomString(Atom atom)
    {
        double closestDistance = 100.0;
        Atom closestAtom = null;
        for (Atom a : contents)
            {
                if ( a.tinkerAtomType != atom.tinkerAtomType )
                    continue;
                double distance = Vector3D.distance(atom.position, a.position);
                if ( distance < closestDistance )
                    {
                        closestDistance = distance;
                        closestAtom = a;
                    }
            }
        if ( closestDistance > 0.1 || closestAtom == null )
            return atom.toString();
        else
            return getAtomString(closestAtom);
    }

    /**
     * Returns a human-readable string for a bond like C5-C6 (1.0) where
     * C5 and C6 are bonded together with an order of 1.0.
     * Throws an IllegalArgumentException if the edge is not in the molecule.
     * @param e the bond
     * @return the string
     */
    public String getBondString(DefaultWeightedEdge e)
    {
        if ( ! connectivity.containsEdge(e) )
            throw new IllegalArgumentException("cannot print out bond string because the corresponding edge is not in graph!");
        Atom atom1 = connectivity.getEdgeSource(e);
        Atom atom2 = connectivity.getEdgeTarget(e);
        double bondOrder = connectivity.getEdgeWeight(e);
        return String.format("%s-%s (.1f)", getAtomString(atom1), getAtomString(atom2), bondOrder);
    }

    /**
     * Returns the number of atoms in this Molecule.
     * @return the size of this molecule
     */
    public int getSize()
    {
        return contents.size();
    }

    /**
     * Given a single atom, returns a set of all connected atoms, including the
     * specified atom.  Uses a breadth-first search algorithm.
     * @param startingAtom the atom to start exploring the connectivity from
     * @return the atoms in the subgraph of the startingAtom
     */
    public Set<Atom> exploreGraph(Atom startingAtom)
    {
        if ( ! contents.contains(startingAtom) )
            throw new IllegalArgumentException("cannot search the connectivity graph because the specified atom is not in this molecule");
        Set<Atom> returnSet = new HashSet<>();
        LinkedList<Atom> searchQueue = new LinkedList<>();
        searchQueue.add(startingAtom);

        // breadth-first search
        while (searchQueue.size() > 0)
            {
                Atom currentNode = searchQueue.remove();
                for (Atom a : getAdjacentAtoms(currentNode))
                    {
                        if ( ! returnSet.contains(a) )
                            {
                                returnSet.add(a);
                                searchQueue.add(a);
                            }
                    }
            }
        return returnSet;
    }

    /**
     * Given a bond between includeAtom and excludeAtom, returns a set of atoms
     * containing all the atoms on the includeAtom side of the bond, including
     * includeAtom.  Returns an empty set if includeAtom and excludeAtom are not
     * directly bonded.  Will throw an exception if includeAtom and excludeAtom
     * form a ring.
     * @param excludeAtom this atom will not be included in the result
     * @param includeAtom this atom will be included in the result
     * @return the atoms on the includeAtom side of the graph
     */
    public Set<Atom> getHalfGraph(Atom excludeAtom, Atom includeAtom)
    {
        //System.out.println("excluded from result: " + getAtomNumber(excludeAtom));
        //System.out.println("included in result: " + getAtomNumber(includeAtom));
        Set<Atom> returnSet = new HashSet<Atom>();

        // if these atoms are not directly bonded, then return an empty set
        if ( directlyConnected(includeAtom, excludeAtom) == false )
            return returnSet;

        // preform a breadth-first search of one branch of the graph only
        LinkedList<Atom> searchQueue = new LinkedList<Atom>();
        Set<Atom> searched = new HashSet<>();
        for (Atom a : getAdjacentAtoms(includeAtom))
            {
                searchQueue.add(a);
                returnSet.add(a);
            }
        searchQueue.remove(excludeAtom);
        returnSet.remove(excludeAtom);
        returnSet.add(includeAtom);
        searched.add(includeAtom);
        //System.out.println("******************");
        //int i=0;
        Atom lastNode = includeAtom;
        while (searchQueue.size() > 0)
            {
                /*i++;
                System.out.printf("Iteration %d\n", i);
                System.out.println("search queue: ");
                for (Atom a : searchQueue)
                    System.out.print(getAtomString(a) + ", ");
                System.out.println("return set: ");
                for (Atom a : returnSet)
                    System.out.print(getAtomString(a) + ", ");
                System.out.printf("Current node is: %s\n", getAtomString(currentNode));
                System.out.println("\n=========\n");*/
                
                Atom currentNode = searchQueue.remove();
                Set<Atom> adjacent = getAdjacentAtoms(currentNode);
                adjacent.remove(lastNode);
                for (Atom a : adjacent)
                    {   
                        // if the excluded atom is found, this is a ring!
                        if ( a == excludeAtom)
                        {
                            GaussianInputFile gjf = new GaussianInputFile(this);
                            gjf.write("error.gjf");
                            throw new IllegalArgumentException("includeAtom " + getAtomString(includeAtom) +
                                                " and excludeAtom " + getAtomString(excludeAtom) + " cannot form a ring!");
                        }

                        // if this is an atom we haven't already searched, mark it
                        if ( ! searched.contains(a) )
                            {
                                searchQueue.add(a);
                                searched.add(a);
                                returnSet.add(a);
                            }
                    }
            }
/*    System.out.println("includeAtom: " + getAtomString(includeAtom));
    System.out.println("excludeAtom: " + getAtomString(excludeAtom));
 System.out.println("return set: ");
                for (Atom a : returnSet)
                    System.out.print(getAtomString(a) + ", ");
System.out.println();
 GaussianInputFile gjf = new GaussianInputFile(this);
                            gjf.write("error.gjf");
*/
        return(returnSet);
    }

    /** convenience method.  1,2,...n */
    public Set<Atom> getHalfGraph(int index1, int index2)
    {
        Atom atom1 = contents.get(index1-1);
        Atom atom2 = contents.get(index2-1);
        if ( atom1 == null || atom2 == null )
            throw new NullPointerException("atoms not found in this molecule for half graph call");
        return getHalfGraph(atom1, atom2);
    }
    
    /** convenience method that returns atom numbers instead of atoms */
    public Set<Integer> getHalfGraphNumbers(Atom atom1, Atom atom2)
    {
        Set<Atom> set = getHalfGraph(atom1, atom2);
        Set<Integer> returnSet = new HashSet<>();
        for (Atom a : set)
            returnSet.add(getAtomNumber(a));
        return returnSet;
    }

    /**
     * Convenience method that returns atom numbers instead of atoms.
     * @param atomNumber1 the atom number of the first atom (1, 2, ..., n)
     * @param atomNumber2 the atom number of the second atom (1, 2, ..., n)
     * @return the set of atom numbers on the atomNumber2 side
     */
    public Set<Integer> getHalfGraphNumbers(int atomNumber1, int atomNumber2)
    {
        return getHalfGraphNumbers(getAtom(atomNumber1), getAtom(atomNumber2));
    }

    /** returns the atom given an atom number */
    public Atom getAtom(int atomNumber)
    {
        return contents.get(atomNumber-1);
    }

    /**
     * Determines whether atom1 and atom2 share an edge (i.e., are bonded).
     * No exception is thrown if these atoms aren't in the graph.
     * @param atom1 test whether this atom is connected to the other atom 
     * @param atom2 test whether this atom is connected to the other atom
     * @return true if the atoms are bonded
     */
    public boolean directlyConnected(Atom atom1, Atom atom2)
    {
        DefaultWeightedEdge e = connectivity.getEdge(atom1,atom2);
        if ( e == null )
            return false;
        return true;
        /*Set<Atom> set1 = getAdjacentAtoms(atom1);
        Set<Atom> set2 = getAdjacentAtoms(atom2);
        if ( set1.contains(atom2) && set2.contains(atom1) )
            return true;
        return false;*/
    }

    /**
     * Alias method.  Indices: 1, 2, ..., n.  No checks.
     */
    public boolean directlyConnected(int i, int j)
    {
        return directlyConnected(contents.get(i-1), contents.get(j-1));
    }

    /**
     * Returns the distance between two atoms.
     * @param atom1 one of the two atoms
     * @param atom2 one of the two atoms
     * @return the distance between the atoms in angstroms
     */
    public double getDistance(Atom atom1, Atom atom2)
    {
        return Vector3D.distance(atom1.position, atom2.position);
    }

    /**
     * Returns the angle between three atoms.
     * @param atom1 one of the three atoms
     * @param atom2 one of the three atoms
     * @param atom3 one of the three atoms
     * @return the angle between the atoms in degrees
     */
    public static double getAngle(Atom atom1, Atom atom2, Atom atom3)
    {
        Vector3D v1 = atom1.position.subtract(atom2.position);
        Vector3D v3 = atom3.position.subtract(atom2.position);
        return Math.toDegrees(Vector3D.angle(v1, v3));
    }

    /** returns the angle between three vectors */
    public static double getAngle(Vector3D v1, Vector3D v2, Vector3D v3)
    {
        Vector3D v1prime = v1.subtract(v2);
        Vector3D v3prime = v3.subtract(v2);
        return Math.toDegrees(Vector3D.angle(v1prime, v3prime));
    }

    /** convenience method 1,2,...,N*/
    public double getAngle(int atom1number, int atom2number, int atom3number)
    {
        return getAngle(contents.get(atom1number-1),contents.get(atom2number-1),contents.get(atom3number-1));
    }

    /**
     * Moves the group associated with atom2 to the specified distance.
     * Motion occurs along the atom1-atom2 bond vector.  Note that this returns a new
     * molecule.  No checks are made.
     * @param atom1 this atom will be held fixed
     * @param atom2 this atom and anything connected to it will be moved
     * @param requestedDistance the requested distance in Angstroms
     * @return a new Molecule containing the same connectivity but new positions
     */
    public Molecule setDistance(Atom atom1, Atom atom2, double requestedDistance)
    {
       Map<Atom,Atom> atomMap = setDistanceMap(atom1, atom2, requestedDistance);
       return moveAtoms(atomMap); 
    }

    /**
     * Creates an atom map for setDistance() methods.
     */
    protected Map<Atom,Atom> setDistanceMap(Atom atom1, Atom atom2, double requestedDistance)
    {
        // determine which atoms have to be moved
        Set<Atom> toBeMoved = getHalfGraph(atom1, atom2);

        // determine how much to move the atoms
        Vector3D oldPosition1 = atom1.position;
        Vector3D oldPosition2 = atom2.position;
        double currentDistance = getDistance(atom1, atom2);

        Vector3D translateVector = oldPosition2.subtract(oldPosition1);
        double scaling = (requestedDistance - currentDistance)/currentDistance;

        Vector3D requiredTranslation = translateVector.scalarMultiply(scaling);

        Map<Atom,Atom> atomMap = new HashMap<>();
        for (Atom oldAtom : toBeMoved)
            {
                Vector3D oldPosition = oldAtom.position;
                Vector3D newPosition = oldPosition.add(requiredTranslation);
                Atom newAtom         = oldAtom.moveAtom(newPosition);
                atomMap.put(oldAtom, newAtom);
            }
        return atomMap;
    }

    /**
     * Alias method.  Atom indices are 1, 2, ..., n.  No checks.
     */
    public Molecule setDistance(int i, int j, double requestedDistance)
    {
        return setDistance(contents.get(i-1), contents.get(j-1), requestedDistance);
    }

    /**
     * Rotates the atom1-atom2-atom3 angle, moving only atom3 and anything in its
     * attached subgraph.  No checks.
     * @param atom1 will not be moved
     * @param atom2 will not be moved
     * @param atom3 will be moved
     * @param theta rotation in degrees
     */
    public Molecule rotateAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        Map<Atom,Atom> atomMap2 = rotateAngleMap(atom1, atom2, atom3, theta);
        return moveAtoms(atomMap2);
    }

    /**
    * Creates an atom map for rotateAngle() methods.
    */
    protected Map<Atom,Atom> rotateAngleMap(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        // figure out which atoms to move
        Set<Atom> toBeMoved = getHalfGraph(atom2, atom3);

        // create atom map
        LinkedHashMap<Atom,Atom> atomMap = new LinkedHashMap<Atom,Atom>();

        // move everything to put atom2 at the origin
        Vector3D v1 = null;
        Vector3D v3 = null;
        for (Atom a : contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom2.position);
                if ( toBeMoved.contains(a) )
                    atomMap.put(a,a.moveAtom(newPosition));
                if ( a == atom1 )
                    v1 = newPosition;
                else if ( a == atom3 )
                    v3 = newPosition;
            }
        // form the rotation axis and matrix
        Vector3D rotationAxis = Vector3D.crossProduct(v1, v3);
        Rotation rotation = new Rotation(rotationAxis, Math.toRadians(theta));

        // apply rotation and undo translation
        LinkedHashMap<Atom,Atom> atomMap2 = new LinkedHashMap<Atom,Atom>();
        for (Atom a : atomMap.keySet())
            {
                Vector3D oldPosition = atomMap.get(a).position;
                Vector3D newPosition = rotation.applyTo(oldPosition);
                newPosition = newPosition.add(atom2.position);
                atomMap2.put(a,a.moveAtom(newPosition));
            }

        return atomMap2;
    }

    /**
     * Method alias.  Indices are 1,2,...,n.  No checks.
     */
    public Molecule rotateAngle(int i, int j, int k, double theta)
    {
        return rotateAngle(contents.get(i-1), contents.get(j-1), contents.get(k-1), theta);
    }

    /**
     * Set the atom1-atom2-atom3 angle to theta degrees, moving atom3 and its subgraph only.
     * New molecule returned.  No checks.
     * @param atom1 not moved
     * @param atom2 not moved
     * @param atom3 moved
     * @param theta desired angle in degrees
     */
    public Molecule setAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        double currentAngle = getAngle(atom1, atom2, atom3);
        double requiredRotation = theta - currentAngle;
        return rotateAngle(atom1, atom2, atom3, requiredRotation);
    }

    /**
     * Method alias.  Indices are 1,2,...,n.  No checks.
     */
    public Molecule setAngle(int i, int j, int k, double theta)
    {
        return setAngle(contents.get(i-1), contents.get(j-1), contents.get(k-1), theta);
    }

    /**
     * Returns a new Molecule with a rotated dihedral.
     * Note that the old AtomTorsion will no longer point to the new Molecule.
     * @param theta the desired dihedral angle in degrees
     * @return the new molecule
     */
    public Molecule setDihedral(AtomTorsion atomTorsion, double theta)
    {
        Map<Atom,Atom> atomMap2 = setDihedralMap(atomTorsion, theta);
        return moveAtoms(atomMap2);
    }

    /**
    * Creates an atom map for setDihedral() methods.
    */
    protected Map<Atom,Atom> setDihedralMap(AtomTorsion atomTorsion, double theta)
    {
         // check that this AtomTorsion is the correct one for this Molecule
        if ( atomTorsion.molecule != this )
            throw new IllegalArgumentException("this isn't the right Molecule for this AtomTorsion");

        // get fields
        Atom atom1 = atomTorsion.atom1;
        Atom atom2 = atomTorsion.atom2;
        Atom atom3 = atomTorsion.atom3;
        Atom atom4 = atomTorsion.atom4;
        List<Atom> atomsToRotate = atomTorsion.atomsToRotate;

        // determine how much rotation is needed
        double currentDihedralAngle = atomTorsion.getDihedralAngle();
        double requiredRotation = currentDihedralAngle - theta;

        LinkedHashMap<Atom,Atom> atomMap = new LinkedHashMap<>();
        
        // move atom 3 to the origin
        // define the rotation axis as the vector from atom3 (now at origin) to atom2
        Vector3D rotationAxis = null;
        for (Atom a : contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom3.position);
                if ( atomsToRotate.contains(a) )
                    atomMap.put(a,a.moveAtom(newPosition));
                if ( a == atom2 )
                    rotationAxis = newPosition;
            }

        // rotate the atoms and make a new atom map
        LinkedHashMap<Atom,Atom> atomMap2 = new LinkedHashMap<>();
        Rotation rotation = new Rotation(rotationAxis, Math.toRadians(requiredRotation));
        for (Atom a : atomMap.keySet())
            {
                // update rotation
                Vector3D oldPosition = atomMap.get(a).position;
                Vector3D newPosition = rotation.applyTo(oldPosition);
                
                // undo translation
                newPosition = newPosition.add(atom3.position);

                // update map
                atomMap2.put(a, a.moveAtom(newPosition));
            }

        return atomMap2;
    }

    /**
     * Alias method.  
     */
    public Molecule setDihedral(ProtoTorsion protoTorsion, double theta)
    {
        AtomTorsion atomTorsion = protoTorsion.getAtomTorsion(this);
        return setDihedral(atomTorsion, theta);
    }

    /**
     * Returns a new Molecule with a rotated dihedral.  Alias method.
     * Note that the old IndexTorsion will still be valid for the new Molecule.
     */
    public Molecule setDihedral(IndexTorsion indexTorsion, double theta)
    {
        return setDihedral(indexTorsion.getAtomTorsion(this), theta);
    }

    /**
     * Creates a new Molecule where atom2 and its subgraph have been moved to
     * make atom1 sp2-hybridized (bond angles set at 120 degrees).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @param forceAngle true if we want to force the atom1alpha-atom1-atom1beta angle to 120 (safe for non-prolines)
     * @return a new Molecule with adjusted hybridization
     */
    public Molecule set_sp2(Atom atom1, Atom atom2, boolean forceAngle)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);
        Map<Atom,Atom> newAtomMap = set_sp2_map(atom1, atom2, forceAngle);
        Molecule rotatedMolecule = moveAtoms(newAtomMap);

        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        Molecule returnMolecule = rotatedMolecule.setDistance(atom1number, atom2number, currentLength); 
        return returnMolecule;
    }

    /**
    * Returns an atom map for setting atom 1 to sp2 geometry.  Does NOT include
    * bond length adjustment.
    */
    protected Map<Atom,Atom> set_sp2_map(Atom atom1, Atom atom2, boolean forceAngle)
    {
        // get the neighbors of atom1
        // call them atom1alpha and atom1beta
        List<Atom> atom1neighbors = new LinkedList<>(getAdjacentAtoms(atom1));
        if ( atom1neighbors.size() != 3 )
            throw new IllegalArgumentException("expected 3 neighbors for atom 1, found " + atom1neighbors.size());
        else if ( ! atom1neighbors.contains(atom2) )
            throw new IllegalArgumentException("atoms are not adjacent");
        atom1neighbors.remove(atom2);
        Atom atom1alpha = atom1neighbors.get(0);
        Atom atom1beta = atom1neighbors.get(1);
        int atom1alphaNumber = getAtomNumber(atom1alpha);
        int atom1betaNumber = getAtomNumber(atom1beta);
        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        // force the existing bond angle to 120 degrees
        Molecule newMolecule = this;
        if (forceAngle)
            newMolecule = setAngle(atom1alpha, atom1, atom1beta, 120.0);
        //System.out.println(getAtomString(atom1alpha));
        //System.out.println(getAtomString(atom1));
        //System.out.println(getAtomString(atom1beta));

        // translate atom1 to the origin
        List<Vector3D> newPositions = new LinkedList<>();
        for (Atom a : newMolecule.contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom1.position);
                newPositions.add(newPosition);
            }

        // get unit vectors for atom1alpha and atom1beta
        Vector3D atom1alphaPosition = newPositions.get(atom1alphaNumber-1).normalize();
        Vector3D atom1betaPosition = newPositions.get(atom1betaNumber-1).normalize();

        // get the cross product of atom1alpha and atom1beta
        Vector3D atom1crossPosition = Vector3D.crossProduct(atom1alphaPosition, atom1betaPosition).normalize();

        // find the linear transformation matrix that rotates atom1alphaPosition to a=(-sqrt(3)/2, 0.5, 0.0)
        // and atom1betaPosition to b=(-0.5, sqrt(3)/2, 0.0). If these are two points of an equilateral triangle,
        // the third is at c=(1,0,0).

        // solve the simultaneous matrix equations:
        // T atom1alphaPosition = a
        // T atom1betaPosition = b
        // T atom1cross = c
        //
        // call atom1alphaPosition A, atom1betaPosition B, and atom1cross C.
        // this is equivalent to solving:
        //
        // T [ ABC ] = [ abc ]
        //
        // where this means concatenated column vectors. Therefore,
        // T = [ abc ] [ ABC ]^-1

        double[][] preMatrix_abc = { { -0.5, Math.sqrt(3.0)/2.0, 0.0 }, {-0.5, -1.0 * Math.sqrt(3.0)/2.0, 0.0}, {0.0, 0.0, 1.0} };
        Matrix matrix_abc = new Matrix(preMatrix_abc);
        matrix_abc = matrix_abc.transpose();

        double[][] preMatrix_ABC = { { atom1alphaPosition.getX(), atom1alphaPosition.getY(), atom1alphaPosition.getZ() },
                                     { atom1betaPosition.getX(), atom1betaPosition.getY(), atom1betaPosition.getZ() },
                                     { atom1crossPosition.getX(), atom1crossPosition.getY(), atom1crossPosition.getZ() } };
        Matrix matrix_ABC = new Matrix(preMatrix_ABC);
        matrix_ABC = matrix_ABC.transpose();

        Matrix matrix_ABC_inverse = matrix_ABC.inverse();
        Matrix T = matrix_abc.times(matrix_ABC_inverse);
        Matrix Tinverse = T.inverse();

        // apply the inverse of T to (1,0,0) to get the third vertex of the triangle
        double[][] preMatrix_c = { { 1.0, 0.0, 0.0 } };
        Matrix matrix_c = new Matrix(preMatrix_c);
        matrix_c = matrix_c.transpose();

        Matrix thirdVertex = Tinverse.times(matrix_c);
        Vector3D thirdVertexPosition = new Vector3D( thirdVertex.get(0, 0), thirdVertex.get(1, 0), thirdVertex.get(2, 0) );

        //double angle1 = getAngle(newPositions.get(atom1alphaNumber-1), newPositions.get(atom1number-1), newPositions.get(atom1betaNumber-1));
        //double angle2 = getAngle(newPositions.get(atom1betaNumber-1), newPositions.get(atom1number-1), newPositions.get(atom2number-1));
        //double angle3 = getAngle(newPositions.get(atom2number-1), newPositions.get(atom1number-1), newPositions.get(atom1alphaNumber-1));
        //System.out.println(angle1);
        //System.out.println(angle2);
        //System.out.println(angle3);

        // calculate the necessary rotation
        Vector3D atom2position = newPositions.get(atom2number-1);
        Vector3D rotationAxis = Vector3D.crossProduct( thirdVertexPosition, atom2position );
        double requiredTheta = Vector3D.angle( thirdVertexPosition, atom2position );
        Rotation fixRotation = new Rotation( rotationAxis, -1.0 * requiredTheta );

        // determine which atoms should be moved
        Set<Integer> atomNumbersToMove = getHalfGraphNumbers(atom1, atom2);
        List<Vector3D> newPositions2 = new LinkedList<>();

        for (int i=0; i < newPositions.size(); i++)
            {
                Integer currentAtomNumber = Integer.valueOf(i+1);

                if ( atomNumbersToMove.contains(currentAtomNumber) )
                    {
                        // rotate this atom
                        Vector3D oldPosition = newPositions.get(i);
                        Vector3D newPosition = fixRotation.applyTo(oldPosition);
                        newPositions2.add(newPosition);
                    }
                else
                    {
                        // do not rotate this atom
                        newPositions2.add( newPositions.get(i) );
                    }
            }

        // undo translation
        List<Vector3D> newPositions3 = new LinkedList<>();
        for (Vector3D v : newPositions2)
            newPositions3.add(v.add(atom1.position));

        //angle1 = getAngle(newPositions3.get(atom1alphaNumber-1), newPositions3.get(atom1number-1), newPositions3.get(atom1betaNumber-1));
        //angle2 = getAngle(newPositions3.get(atom1betaNumber-1), newPositions3.get(atom1number-1), newPositions3.get(atom2number-1));
        //angle3 = getAngle(newPositions3.get(atom2number-1), newPositions3.get(atom1number-1), newPositions3.get(atom1alphaNumber-1));
        //System.out.println(angle1);
        //System.out.println(angle2);
        //System.out.println(angle3);
        //System.out.println(atom1alphaNumber + " " + atom1number + " " + atom1betaNumber);
        //System.out.println(atom1betaNumber + " " + atom1number + " " + atom2number);
        //System.out.println(atom2number + " " + atom1number + " " + atom1alphaNumber);

        // create new atom map
        Map<Atom,Atom> newAtomMap = new LinkedHashMap<>();
        for (int i=0; i < contents.size(); i++)
            {
                Integer currentAtomNumber = Integer.valueOf(i+1);
                if ( !atomNumbersToMove.contains(currentAtomNumber) )
                    continue;
                Atom oldAtom = contents.get(i);
                Atom newAtom = oldAtom.moveAtom( newPositions3.get(i) );
                if ( !oldAtom.equals(newAtom) )
                    newAtomMap.put(oldAtom, newAtom);
            }
        return newAtomMap;
    }
    
    /**
    * Alias method.  
    */
    public Molecule set_sp2(Atom atom1, Atom atom2)
    {
        return set_sp2(atom1, atom2, true);
    }

    /**
     * Creates a new Molecule where atom2 and its subgraph have been moved to
     * make atom1 sp3-hybridized (tetrahedral geometry).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @return a new Molecule with adjusted hybridization
     */
    public Molecule set_sp3(Atom atom1, Atom atom2)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);
        Map<Atom,Atom> newAtomMap = set_sp3_map(atom1, atom2);
        Molecule rotatedMolecule = moveAtoms(newAtomMap);

        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        Molecule returnMolecule = rotatedMolecule.setDistance(atom1number, atom2number, currentLength); 
        return returnMolecule;
    }

    /**
    * Returns an atom map for setting atom 1 to sp3 geometry.  Does NOT include
    * bond length adjustment.
    */
    protected Map<Atom,Atom> set_sp3_map(Atom atom1, Atom atom2)
    {
        // get the neighbors of atom1
        // call them atom1alpha, atom1beta, atom1gamma
        List<Atom> atom1neighbors = new LinkedList<>(getAdjacentAtoms(atom1));
        if ( atom1neighbors.size() != 4 )
            throw new IllegalArgumentException("expected 4 neighbors for atom 1, found " + atom1neighbors.size());
        else if ( ! atom1neighbors.contains(atom2) )
            throw new IllegalArgumentException("atoms are not adjacent");
        atom1neighbors.remove(atom2);
              Atom atom1alpha = atom1neighbors.get(0);
        Atom atom1beta = atom1neighbors.get(1);
	    Atom atom1gamma = atom1neighbors.get(2);
        int atom1alphaNumber = getAtomNumber(atom1alpha);
        int atom1betaNumber = getAtomNumber(atom1beta);
	    int atom1gammaNumber = getAtomNumber(atom1gamma);
        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        Molecule newMolecule = this;

        // translate atom1 to the origin
        List<Vector3D> newPositions = new LinkedList<>();
        for (Atom a : newMolecule.contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom1.position);
                newPositions.add(newPosition);
            }

        // get unit vectors
        Vector3D atom1alphaPosition = newPositions.get(atom1alphaNumber-1).normalize();
        Vector3D atom1betaPosition = newPositions.get(atom1betaNumber-1).normalize();
	Vector3D atom1gammaPosition = newPositions.get(atom1gammaNumber-1).normalize();
	Vector3D atom1atom2Position = newPositions.get(atom2number-1).normalize();

        // get the negative sum of atom1alpha, atom1beta, atom1gamma
        Vector3D targetPosition = atom1alphaPosition.add(atom1betaPosition.add(atom1gammaPosition));
        targetPosition = targetPosition.negate();
        Rotation fixRotation = new Rotation(atom1atom2Position, targetPosition);

        // determine which atoms should be moved
        Set<Integer> atomNumbersToMove = getHalfGraphNumbers(atom1, atom2);
        List<Vector3D> newPositions2 = new LinkedList<>();

        for (int i=0; i < newPositions.size(); i++)
            {
                Integer currentAtomNumber = Integer.valueOf(i+1);

                if ( atomNumbersToMove.contains(currentAtomNumber) )
                    {
                        // rotate this atom
                        Vector3D oldPosition = newPositions.get(i);
                        Vector3D newPosition = fixRotation.applyTo(oldPosition);
                        newPositions2.add(newPosition);
                    }
                else
                    {
                        // do not rotate this atom
                        newPositions2.add( newPositions.get(i) );
                    }
            }

        // undo translation
        List<Vector3D> newPositions3 = new LinkedList<>();
        for (Vector3D v : newPositions2)
            newPositions3.add(v.add(atom1.position));

        // create new atom map
        Map<Atom,Atom> newAtomMap = new LinkedHashMap<>();
        for (int i=0; i < contents.size(); i++)
            {
                Integer currentAtomNumber = Integer.valueOf(i+1);
                if ( !atomNumbersToMove.contains(currentAtomNumber) )
                    continue;
                Atom oldAtom = contents.get(i);
                Atom newAtom = oldAtom.moveAtom( newPositions3.get(i) );
                if ( !oldAtom.equals(newAtom) )
                    newAtomMap.put(oldAtom, newAtom);
            }
        return newAtomMap;
    }

    /**
     * Returns the bonded neighbors of includeAtom.  Does not include includeAtom itself.
     * @param includeAtom the atom whose neighbors are to be searched
     * @return the Atoms adjacent to includeAtom
     */
    public Set<Atom> getAdjacentAtoms(Atom includeAtom)
    {
        Set<Atom> returnSet = new HashSet<Atom>();
        if ( ! connectivity.containsVertex(includeAtom) )
            throw new IllegalArgumentException("includeAtom must be within this connectivity graph!");
        for (DefaultWeightedEdge e : connectivity.edgesOf(includeAtom))
            {
                returnSet.add(connectivity.getEdgeSource(e));
                returnSet.add(connectivity.getEdgeTarget(e));
            }
        returnSet.remove(includeAtom);
        return(returnSet);
    }

    /**
     * Returns determine if Atom is adjacent to another Atom of element e.
     * @param includeAtom the atom whose neighbors are to be searched
     * @return true if includeAtom is adjacent to another Atom of element e, false otherwise.
     */
    public boolean adjacentTo(Atom includeAtom, Element e)
    {
        for (Atom a : getAdjacentAtoms(includeAtom)) if (a.element.equals(e)) return true;
        return false;
    }

    /**
     * Returns the bond between atom1 and atom2.  Throws an IllegalArgumentException if there
     * is no such edge.
     * @param atom1 one of the atoms
     * @param atom2 the other atom
     * @return the bond
     */
    public DefaultWeightedEdge getBond(Atom atom1, Atom atom2)
    {
        DefaultWeightedEdge returnEdge = connectivity.getEdge(atom1, atom2);
        if (returnEdge == null)
            throw new IllegalArgumentException("cannot return a bond; atoms not bonded!");
        return returnEdge;
    }

    /**
     * Checks if atom1 and atom2 are more than two bonds apart.
     * @param atom1 the first atom
     * @param atom2 the second atom
     * @return true if atom1 and atom2 are separated by three or more bonds
     */
    public boolean areSeparated(Atom atom1, Atom atom2)
    {
        Set<Atom> atom1neighbors = getAdjacentAtoms(atom1);
        
        // check if direct neighbors
        if (atom1neighbors.contains(atom2))
            return false;
        Set<Atom> atom2neighbors = getAdjacentAtoms(atom2);
       
        // check if geminal
        for (Atom a : atom1neighbors)
            if (atom2neighbors.contains(a))
                return false;
        return true;
    }

    /**
     * Checks if the atoms in this molecule clash.
     * @return true if there is a self-clash in this molecule
     */
    public boolean tooClose()
    {
        for (int i=0; i < contents.size(); i++)
            {
                Atom atom1 = contents.get(i);
                for (int j=i+1; j < contents.size(); j++)
                    {
                        Atom atom2 = contents.get(j);
                        double distance = Vector3D.distance(atom1.position, atom2.position);
                        if ( distance > 1.20 )
                            continue;
                        // ignore atoms that are 1,2 or 1,3
                        if ( !areSeparated(atom1,atom2) )
                            continue;
                        
                        //for debugging
                        //System.out.printf("Atoms too close: %s %s %s %s \n", getAtomString(atom1), atom1.toString(), getAtomString(atom2), atom2.toString());
                        return true;
                    }
            }
        return false;
    }

    /**
     * Computes a rough Lennard-Jones steric energy for this molecule.  Answer is normalized
     * by the number of atoms; i.e. kcal/mol divided by the number of atoms.  Atoms that are
     * separated by one or two bonds are ignored.  Atoms separated by more than Settings.CUTOFF_DISTANCE
     * are also ignored.
     * @return the steric energy
     */
    public double getOPLSenergy()
    {
        double energy = 0.0;
        for (int i=0; i < contents.size(); i++)
            {
                Atom atom1 = contents.get(i);
                Vector3D atom1position = atom1.position;
                for (int j=i+1; j < contents.size(); j++)
                    {
                        Atom atom2 = contents.get(j);
                        // ignore if atoms are too close in the connectivity graph
                        if ( ! areSeparated(atom1,atom2) )
                            continue;

                        // ignore if atoms are too far apart
                        Vector3D atom2position = atom2.position;
                        double distance = Vector3D.distance(atom1position, atom2position);
                        if ( distance > Settings.CUTOFF_DISTANCE )
                            continue;

                        // prevent overflow
                        if ( distance < 0.5 )
                            distance = 0.5;

                        // get parameters
                        double epsilon1 = atom1.element.epsilon;
                        double sigma1   = atom1.element.sigma;
                        double epsilon2 = atom2.element.epsilon;
                        double sigma2   = atom2.element.sigma;

                        // apply combination rules
                        double epsilon = epsilon1;
                        double sigma   = sigma1;
                        if ( epsilon1 != epsilon2 )
                            epsilon = Math.sqrt(epsilon1 * epsilon2);
                        if ( sigma1 != sigma2 )
                            sigma = Math.sqrt(sigma1 * sigma2);

                        // compute energy
                        double temp = Math.pow(sigma/distance, 6);
                        energy += 4.0 * epsilon * temp * (temp - 1.0);
                    }
            }
        return energy / contents.size();
    }

    /**
     * Checks if atoms are too close in a molecule, given another molecule
     * whose atoms we know are not too close.  Intended for assessing the result
     * of dihedral changes.  Does not assume the molecules have the same composition.
     * @param oldMolecule the molecule this molecule was modified from
     * @return true if there is at least one atom that is too close to another atom
     */
    public boolean checkCloseContacts(Molecule oldMolecule)
    {
        // the atoms that have not changed in the new peptide
        ArrayList<Atom> oldAtoms = new ArrayList<>();

        // the atoms that are new in newPeptide
        ArrayList<Atom> newAtoms = new ArrayList<>();

        // populate lists
        for (Atom a : contents)
            {
                if ( oldMolecule.contents.contains(a) )
                    oldAtoms.add(a);
                else
                    newAtoms.add(a);
            }

        // compare distances between old atoms and new atoms
        for (Atom oldAtom : oldAtoms)
            {
                Vector3D oldPosition = oldAtom.position;
                for (Atom newAtom : newAtoms)
                    {
                        Vector3D newPosition = newAtom.position;
                        double distance = Vector3D.distance(oldPosition,newPosition);
                        if ( distance < Settings.MINIMUM_DISTANCE &&
                             connectivity.getEdge(oldAtom, newAtom) == null )
                            return true;
                    }
            }
        return false;
    }

    /**
     * Checks if the atoms are too close in a molecule.  The minimum distance
     * is controlled by Settings.MINIMUM_DISTANCE.  Computations are performed
     * on the triangular distance matrix.
     * @return true if there is at least one atom that is too close to another atom
     */
    public boolean checkCloseContacts()
    {
        // compute the upper triangle of distance contacts
        for (int i=0; i < contents.size(); i++)
            {
                Atom atom1 = contents.get(i);
                Vector3D position1 = atom1.position;
                for (int j=i+1; j < contents.size(); j++)
                    {
                        Atom atom2 = contents.get(j);
                        Vector3D position2 = atom2.position;

                        // ignores distances between directly connected atoms
                        if ( Vector3D.distance(position1, position2) < Settings.MINIMUM_DISTANCE &&
                             connectivity.getEdge(atom1,atom2) == null )
                            return true;
                    }
            }
        return false;
    }

    /**
     * Check if this Molecule is too close to other.
     * @param other the other molecule to which we compare this one
     * @return true if an Atom of this is too close to an Atom of other, and false otherwise
     */
    public boolean tooClose(Molecule other)
    {
        for (Atom a : this.contents) {
            for (Atom b : other.contents) {
                if (Vector3D.distance(a.position,b.position) < 2.0) return true;
            }
        }
        return false;
    }

    /**
     * Check if this Molecule is too far from other.
     * @param other the other molecule to which we compare this one
     * @return true if every Atom of this is too far from every Atom of other, and false otherwise
     */
    public boolean tooFar(Molecule other)
    {
        for (Atom a : this.contents) {
            for (Atom b : other.contents) {
                if (Vector3D.distance(a.position,b.position) < 2.5) return false;
            }
        }
        return true;
    }

    /**
     * Check if this Molecule does not meet other in an interesting way.
     * Interesting is defined as an intermolecular hydrogen bond.
     * Implementation assumes there is only one adjacent atom to each hydrogen.
     * @param other the other molecule to which we compare this one
     * @return true if this is an interesting arrangement 
     */
    public boolean isInteresting(Molecule other)
    {
        Element O = Element.getElement("O");
        Element N = Element.getElement("N");
        Element H = Element.getElement("H");

        Set<Element> polarElements = ImmutableSet.of(O,N);

        double desiredMaxDistance = 3.0;
        double desiredMinAngle = 120.0;

        // this   other
        // X: ... H-X
        // a      b c
        for (Atom a : this.contents)
            {
                if ( ! polarElements.contains(a.element) )
                    continue;

                // ignore nitrogens with three neighbors
                if ( a.element.equals(N) && getAdjacentAtoms(a).size() > 2 )
                    continue;

                for (Atom b : other.contents)
                    {
                        // check if this atom is a polar hydrogen
                        if ( ! b.element.equals(H) )
                            continue;

                        boolean nextToPolar = false;
                        Atom polarAtom = null;
                        for (Atom c : other.getAdjacentAtoms(b))
                            {
                                if ( polarElements.contains(c.element) )
                                    {
                                        nextToPolar = true;
                                        polarAtom = c;
                                        break;
                                    }
                            }
                        if ( ! nextToPolar )
                            continue;

                        double thisDistance = Vector3D.distance(a.position, b.position);
                        double thisAngle = getAngle(a,b,polarAtom);
                        if ( thisDistance < desiredMaxDistance && thisAngle > desiredMinAngle )
                            return true; 
                    }
            }

        // this   other
        // X-H ... X
        // a b     c
        for (Atom b : this.contents)
            {
                // check if this atom is a polar hydrogen
                if ( ! b.element.equals(H) )
                    continue;

                boolean nextToPolar = false;
                Atom polarAtom = null;
                for (Atom a : getAdjacentAtoms(b) )
                    {
                        if ( polarElements.contains(a.element) )
                            {
                                nextToPolar = true;
                                polarAtom = a;
                            }
                    }

                if ( ! nextToPolar )
                    continue;

                for (Atom c : other.contents)
                    {
                        if ( ! polarElements.contains(c.element) )
                            continue;

                        // ignore nitrogens with three neighbors
                        if ( c.element.equals(N) && other.getAdjacentAtoms(c).size() > 2 )
                            continue;

                        double thisDistance = Vector3D.distance(b.position, c.position);
                        double thisAngle = getAngle(polarAtom,b,c);
                        if ( thisDistance < desiredMaxDistance && thisAngle > desiredMinAngle )
                            return true;
                    }
            }

        return false;
    }

    /**
     * Calculates the "distance" between two molecules by comparing their positions.
     * @param molecule1 the first molecule
     * @param molecule2 the second molecule
     * @return the RMS average distance between the atoms in molecule1 and molecule2 in angstroms
     */
    public static double calculateRMSD(Molecule molecule1, Molecule molecule2)
    {
        if ( molecule1.contents.size() != molecule2.contents.size() )
            throw new IllegalArgumentException("molecules are not the same size");
        
        double RMSD = 0.0;
        for (int i=0; i < molecule1.contents.size(); i++)
            {
                Vector3D fromVector = molecule1.contents.get(i).position;
                Vector3D toVector = molecule2.contents.get(i).position;
                RMSD += Math.pow(Vector3D.distance(fromVector,toVector), 2);
            }
        RMSD = RMSD * ( 1.0 / molecule1.contents.size() );
        RMSD = Math.sqrt(RMSD);
        return RMSD;
    }

    /**
     * Calculates the "distance" between two molecules by comparing their positions.
     * Molecules are superimposed first.
     * @param molecule1 the first molecule
     * @param molecule2 the second molecule
     * @param atomNumbers the atom numbers (1...N) of the atoms to be compared
     * @return the RMS average distance between the superimposed molecules, counting only the atoms in the lsit
     */
    public static double calculateRMSD(Molecule molecule1, Molecule molecule2, List<Integer> atomNumbers)
    {
        Molecule m2 = superimpose(molecule1,molecule2,atomNumbers);

        double RMSD = 0.0;
        for (Integer i : atomNumbers)
            {
                int atomIndex = i-1;
                Atom fromAtom = molecule1.contents.get(atomIndex);
                Atom toAtom = m2.contents.get(atomIndex);
                Vector3D fromVector = fromAtom.position;
                Vector3D toVector = toAtom.position;
                RMSD += Math.pow(Vector3D.distance(fromVector,toVector), 2);
            }
        RMSD = RMSD * ( 1.0 / molecule1.contents.size() );
        RMSD = Math.sqrt(RMSD);
        return RMSD;
    }

    /**
     * Creates an atom map that matches atoms to new ones by index.  For instance
     * this is useful if you want to do operations only on the molecule part of
     * a class, then rectify the contents of the class to match the product.
     * @param newMolecule a molecule with the same number of atoms,
                indexed in the same order
     * @return an atom map
     */
     public Map<Atom,Atom> matchMap(Molecule newMolecule)
     {
        Map<Atom,Atom> returnMap = new HashMap<>();
        if ( this.contents.size() != newMolecule.contents.size())
            throw new IllegalArgumentException("Sizes of molecules do not match!");
        for ( int i = 0; i < this.contents.size(); i++ )
            returnMap.put(this.getAtom(i+1), newMolecule.getAtom(i+1));

        return returnMap;
     }

    /**
     * Superimposes two molecules based on a list of atom numbers.
     * Superposition is performed using the
     * <a href="http://en.wikipedia.org/wiki/Kabsch_algorithm">Kabsch algorithm</a>.
     * @param molecule1 one of the molecules
     * @param molecule2 another molecule (superposition transformation will be applied to this vector)
     * @param atomNumbers the atom numbers (1...N) to use for the superposition
     * @return a new Molecule, which is molecule2 rotated to be maximally superimposed with molecule1
     */
    public static Molecule superimpose(Molecule molecule1, Molecule molecule2, List<Integer> atomNumbers)
    {
        // check validity of molecule arguments
        if ( molecule1.contents.size() != molecule2.contents.size() )
            throw new IllegalArgumentException("molecule size mismatch");

        // check validity of atom numbers list
        // does not check the validity of the numbers themselves
        if ( atomNumbers == null || atomNumbers.size() < 3 || atomNumbers.size() > molecule1.contents.size() )
            throw new IllegalArgumentException("invalid atom numbers list");

        // collect centroids
        Vector3D centroid1 = molecule1.getCentroid();
        Vector3D centroid2 = molecule2.getCentroid();

        // move molecules to origin
        Molecule m1 = molecule1.normalize();
        Molecule m2 = molecule2.normalize();

        // calculate superposition on centered molecules
        double[][] pre_matrixP = new double[atomNumbers.size()][3];
        double[][] pre_matrixQ = new double[atomNumbers.size()][3];
        for (int i=0; i < atomNumbers.size(); i++)
            {
                int atomIndex = atomNumbers.get(i) - 1;
                Atom fromAtom = m1.contents.get(atomIndex);
                if ( fromAtom == null )
                    throw new NullPointerException("error in atom index (" + atomIndex + ")");
                pre_matrixP[i][0] = fromAtom.position.getX();
                pre_matrixP[i][1] = fromAtom.position.getY();
                pre_matrixP[i][2] = fromAtom.position.getZ();

                Atom toAtom = m2.contents.get(atomIndex);
                pre_matrixQ[i][0] = toAtom.position.getX();
                pre_matrixQ[i][1] = toAtom.position.getY();
                pre_matrixQ[i][2] = toAtom.position.getZ();
            }

        Matrix matrixP = new Matrix(pre_matrixP);
        Matrix matrixQ = new Matrix(pre_matrixQ);
        Matrix matrixPtranspose = matrixP.transpose();

        // calculate covariance matrix and its singular value decomposition
        Matrix covarianceMatrix = matrixPtranspose.times(matrixQ);
        SingularValueDecomposition SVD = new SingularValueDecomposition(covarianceMatrix);
        Matrix matrixV = SVD.getU();
        Matrix matrixW = SVD.getV();

        // check sign of rotation matrix
        Double determinant = covarianceMatrix.det();
        double sign = 0.0;
        if (determinant > 0)
            sign = 1.0;
        else if (determinant == 0)
            sign = 1.0;
        else if (determinant < 0)
            sign = -1.0;

        double[][] pre_signedMatrix = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,sign}};
        Matrix signedMatrix = new Matrix(pre_signedMatrix);

        // compute optimal rotation matrix
        Matrix rotationMatrix = matrixV.times(signedMatrix);
        rotationMatrix = rotationMatrix.times(matrixW.transpose());

        // translate between packages by converting Vector3D objects to Matrix objects
        double[][] currentVectorPreMatrix = new double[3][1];
        double[][] newVectorPreMatrix = new double[3][1];
        
        // apply rotation to all atoms
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (int i=0; i < m2.contents.size(); i++)
            {
                Atom fromAtom = m2.contents.get(i);
                Vector3D currentVector = fromAtom.position;
                currentVectorPreMatrix[0][0] = currentVector.getX();
                currentVectorPreMatrix[1][0] = currentVector.getY();
                currentVectorPreMatrix[2][0] = currentVector.getZ();
                Matrix currentVectorMatrix = new Matrix(currentVectorPreMatrix);

                Matrix newVectorMatrix = rotationMatrix.times(currentVectorMatrix);
                newVectorPreMatrix = newVectorMatrix.getArrayCopy();
                Vector3D newVector = new Vector3D(newVectorPreMatrix[0][0], newVectorPreMatrix[1][0], newVectorPreMatrix[2][0]);
                
                Atom newAtom = fromAtom.moveAtom(newVector);
                atomMap.put(fromAtom, newAtom);
            }

        // rotate molecule and undo translation
        Molecule rotatedMolecule = m2.moveAtoms(atomMap);
        return rotatedMolecule.shift(centroid2);
    }

    /**
     * Returns a deep copy of this molecule with a new name.
     * @param name the new name
     * @return the new molecule
     */
    public Molecule setName(String name)
    {
        Molecule copiedMolecule = moveAtoms(new HashMap<Atom,Atom>());
        return new Molecule(name, copiedMolecule.contents, copiedMolecule.connectivity, copiedMolecule.energy);
    }

    /**
     * Returns the hash code of this molecule.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(name, contents, connectivity, energy);
    }

    /**
     * Checks for object equality by comparing all fields.
     * @param obj the object to be compared to
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Molecule) )
            return false;

        Molecule anotherMolecule = (Molecule)obj;
        if ( this.name.equals(anotherMolecule.name) &&
             this.contents.equals(anotherMolecule.contents) &&
             this.connectivity.equals(anotherMolecule.connectivity) &&
             this.energy == anotherMolecule.energy )
            return true;
        return false;
    }

    /**
     * Returns a textual description of this molecule.
     * @return the string
     */
    @Override
    public String toString()
    {
        String returnString = name + "\n\n";
        for (Atom a : contents)
            returnString = returnString + (contents.indexOf(a) + 1) + a.toString() + "\n";
        return returnString;
    }

    /**
     * Returns a GaussView readable string.
     */
    public String toGaussianString()
    {
        String geometry = "";
        for (Atom a : contents)
            geometry = geometry + a.toString() + "\n";
        return String.format("#\n\n%s\n\n0 1\n%s\n\n", name, geometry, connectivity);
    }

    /**
     * Returns the geometry string for use in Omnisol.
     */
    public String toOmnisolString()
    {
        String geometry = "";
        for (Atom a : contents)
            geometry = geometry + a.toString() + "\n";
        return geometry;
    }

    /**
       Returns the string representation of a Molecule in XYZ files (to be passed to Tinker)
    */
    public String toXYZString() 
    {
	//creates String to write
        String outputString = "";
        //write number of atoms and molecule name
        outputString = contents.size() + " " + name + "\n";

        //write atom list and connections
        int currentAtomNumber = 1;
        for (Atom currentAtom : contents)
            {
                Set<DefaultWeightedEdge> bonds = connectivity.edgesOf(currentAtom);
                outputString = outputString + String.format("%3d %2s %12.8f %12.8f %12.8f %6d", currentAtomNumber,
                                                            currentAtom.element.symbol,  currentAtom.position.getX(),
                                                            currentAtom.position.getY(), currentAtom.position.getZ(),
                                                            currentAtom.tinkerAtomType);

                for (DefaultWeightedEdge e : bonds)
                    {
                        int edgeTarget = contents.indexOf(connectivity.getEdgeTarget(e)) + 1;
                        int edgeSource = contents.indexOf(connectivity.getEdgeSource(e)) + 1;
                        int edgeTargetToWrite = edgeTarget;
                        //change to edgeSource if edgeTarget is simply the current atom (edges have no particular ordering)
                        if (edgeTarget == currentAtomNumber)
                            edgeTargetToWrite = edgeSource;
                        outputString = outputString+String.format("%6d",edgeTargetToWrite);
                    }
                outputString = outputString + "\n";
                currentAtomNumber++;
            }
        return outputString;
    }

    /**
     * Creates a MOL2 String describing this molecule.
     * @return a String containing the geometry and connectivity in MOL2 format
     */
    public String toMOL2()
    {
        // header information
        String returnString = "@<TRIPOS>MOLECULE\n";
        returnString = returnString + name + "\n";

        int numberOfAtoms = contents.size();
        int numberOfBonds = connectivity.edgeSet().size();

        returnString = returnString + numberOfAtoms + " " + numberOfBonds + "\n";
        returnString = returnString + "SMALL\nNO_CHARGES\n\n\n";

        // print geometry data
        returnString = returnString + "@<TRIPOS>ATOM\n";
        
        Atom currentAtom = null;
        String atomSymbolString = "";
        String atomSymbolString2 = "";
        String atomString = "";
        int formalCharge = 0;
        
        for (int i=0; i < contents.size(); i++)
            {
                currentAtom = contents.get(i);
                atomSymbolString = String.format("%s%s", currentAtom.element.symbol, i);
                atomString = String.format("%5s %-6s %9.4f %9.4f %9.4f %5s\n", i+1, atomSymbolString,
                                           currentAtom.position.getX(), currentAtom.position.getY(),
                                           currentAtom.position.getZ(), currentAtom.element.symbol);
                returnString = returnString + atomString;
            }
                
        // write connectivity data
        Atom fromAtom = null;
        Atom toAtom = null;
        Integer fromAtomNumber = 0;
        Integer toAtomNumber = 0;
        Double bondOrder = 0.0;
        int count = 0;
        String bondString = "";
        
        returnString = returnString + "@<TRIPOS>BOND\n";
        for ( DefaultWeightedEdge e : connectivity.edgeSet() )
            {
                count++;
                fromAtom = connectivity.getEdgeSource(e);
                fromAtomNumber = getAtomNumber(fromAtom);
                toAtom = connectivity.getEdgeTarget(e);
                toAtomNumber = getAtomNumber(toAtom);
                bondOrder = connectivity.getEdgeWeight(e);
                if ( bondOrder == 1.5 ) 
                    bondString = String.format("%6s%5s%5s %s\n", count, fromAtomNumber, toAtomNumber, "Ar");
                else
                    bondString = String.format("%6s%5s%5s %s\n", count, fromAtomNumber, toAtomNumber, Math.round(bondOrder));
                returnString = returnString + bondString;
            }
        return returnString;
    }

    public String getGJFstring()
    {
        String returnString = "";
        for (Atom a : contents)
            returnString += a.toString() + "\n";
        return returnString;
    }

    public Molecule addBond(Atom a1, Atom a2)
    {
        if ( (!contents.contains(a1)) || (!contents.contains(a2)) )
            throw new IllegalArgumentException("atom not in graph");
        Molecule newMolecule = this.moveAtoms(new HashMap<Atom,Atom>());
        DefaultWeightedEdge thisEdge = newMolecule.connectivity.addEdge(a1, a2);
        if (thisEdge != null)
            newMolecule.connectivity.setEdgeWeight(thisEdge, 1.0);
        return newMolecule;
    }

    public Molecule addBond(int atomNumber1, int atomNumber2)
    {
        Atom a1 = contents.get(atomNumber1-1);
        Atom a2 = contents.get(atomNumber2-1);
        return addBond(a1,a2);
    }

    public Molecule removeBond(Atom a1, Atom a2)
    {
        if ( (!contents.contains(a1)) || (!contents.contains(a2)) )
            throw new IllegalArgumentException("atom not in graph");
        Molecule newMolecule = this.moveAtoms(new HashMap<Atom,Atom>());
        newMolecule.connectivity.removeEdge(a1, a2);
        return newMolecule;
    }

    public Molecule removeBond(int atomNumber1, int atomNumber2)
    {
        Atom a1 = contents.get(atomNumber1-1);
        Atom a2 = contents.get(atomNumber2-1);
        return removeBond(a1,a2);
    }

    /** for testing */
    public static void main(String[] args)
    {
        /*TinkerMultipleXYZOutputFile file = new TinkerMultipleXYZOutputFile("g09_scans/combined.xyz");
        List<Molecule> molecules = new LinkedList<Molecule>(file.scanPoints.keySet());
        Molecule m1 = molecules.get(0);
        Molecule m2 = molecules.get(1);
        List<Integer> atomNumbers = new LinkedList<>();
        for (int i=0; i < m1.contents.size(); i++)
            atomNumbers.add(i+1);
        System.out.println(calculateRMSD(m1,m1,atomNumbers));
        System.out.println(calculateRMSD(m1,m2,atomNumbers));
        */
        GJFfile file = new GJFfile("input/alkane.gjf");
        Molecule molecule = file.molecule;
        System.out.println(molecule);
        Set<Atom> atoms = molecule.getHalfGraph(7,5);
        Set<Integer> numbers = new TreeSet<>();
        for (Atom a : atoms)
            {
                int number = molecule.contents.indexOf(a) +1;
                numbers.add(number);
            }
        for (Integer number : numbers)
            System.out.printf("%d, ", number);
        System.out.println();
    }
}
