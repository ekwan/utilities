import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * Represents a torsion composed of atom numbers in a molecule.<p>
 * This allows the IndexTorsion to remain usable if the molecular geometry changes. <p>
 * Torsion is: index1-index2-index3-index4<p>
 * Torsion angle calculated between atoms 2 and 3.<p>
 * If moved, the convention is to hold atom1 fixed and allow atom4 to move.<p>
 * This class contains a pre-computed list of the atoms to be rotated and
 * the rotation axis.
 */
public class IndexTorsion extends AbstractTorsion
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** the atom numbers (1,2,...,n) that comprise this torsion */
    public final int index1, index2, index3, index4;

    /** the Atoms including atom4 that should be rotated */
    public final List<Integer> atomNumbersToRotate;

    /** constructor is designed to be called from AtomTorsion */
    public IndexTorsion(int index1, int index2, int index3, int index4, List<Integer> atomNumbersToRotate)
    {
        this.index1 = index1;
        this.index2 = index2;
        this.index3 = index3;
        this.index4 = index4;
        this.atomNumbersToRotate = atomNumbersToRotate;
        //System.out.println(index1 + ", " + index2 + ", " + index3 + ", " + index4);
        //System.out.println(atomNumbersToRotate);
        if ( index1 < 1 || index2 < 1 || index3 < 1 || index4 < 1 )
            throw new IllegalArgumentException("atom numbers must be greater than zero");
        if ( ImmutableSet.of(index1,index2,index3,index4).size() != 4 )
            throw new IllegalArgumentException("duplicate atom numbers");
        if ( atomNumbersToRotate == null || atomNumbersToRotate.size() == 0 )
            throw new NullPointerException("must have some atoms to rotate");
        if ( !atomNumbersToRotate.contains(index4) )
            throw new IllegalArgumentException("must rotate atom4");
        for (Integer i : atomNumbersToRotate)
            if ( i < 1 )
                throw new IllegalArgumentException("rotate atom numbers must be greater than zero");
    }

    @Override
    public String toString()
    {
        return String.format("%d-%d-%d-%d: %s", index1, index2, index3, index4, atomNumbersToRotate.toString());
    }

    /**
     * Static factory method.
     */
    public static IndexTorsion createIndexTorsion(int index1, int index2, int index3, int index4, Molecule molecule)
    {
        Set<Atom> atomsToRotate = molecule.getHalfGraph(index2, index3);
        List<Integer> atomNumbersToRotate = new LinkedList<>();
        for (Atom a : atomsToRotate)
            atomNumbersToRotate.add(molecule.contents.indexOf(a)+1);
        return new IndexTorsion(index1,index2,index3,index4,ImmutableList.copyOf(atomNumbersToRotate));
    }

    /**
     * Makes an IndexTorsion from a ProtoTorsion.
     * @param torsion the ProtoTorsion to construct this IndexTorsion from
     * @param m the molecule the atoms in the ProtoTorsion belong to
     * @return the corresponding IndexTorsion
     */
    public static IndexTorsion createIndexTorsion(ProtoTorsion torsion, Molecule m)
    {
        int index1 = m.contents.indexOf(torsion.atom1);
        int index2 = m.contents.indexOf(torsion.atom2);
        int index3 = m.contents.indexOf(torsion.atom3);
        int index4 = m.contents.indexOf(torsion.atom4);

        if ( index1 == -1 || index2 == -1 || index3 == -1 || index4 == -1)
            throw new IllegalArgumentException("this atom is not in the molecule");

        return createIndexTorsion(index1,index2,index3,index4,m);
    }

    /** returns the corresponding torsion angle */
    public double getDihedralAngle()
    {
        throw new IllegalArgumentException("Mandor needs to know which Molecule this torsion corresponds to to calculate an angle.");
    }

    /**
     * Returns the corresponding torsion angle.  Could throw a NullPointerException.
     */
    public double getDihedralAngle(Molecule molecule)
    {
        Vector3D v1 = molecule.contents.get(index1-1).position; 
        Vector3D v2 = molecule.contents.get(index2-1).position; 
        Vector3D v3 = molecule.contents.get(index3-1).position; 
        Vector3D v4 = molecule.contents.get(index4-1).position;
        return getDihedralAngle(v1,v2,v3,v4);
    }

    public List<Atom> getAtomsToMove(Molecule molecule)
    {
        Atom atom2 = molecule.contents.get(index2-1);
        Atom atom3 = molecule.contents.get(index3-1);
        return ImmutableList.copyOf(molecule.getHalfGraph(atom2, atom3));
    }

    /** returns an IndexTorsion version of this AtomTorsion */
    public AtomTorsion getAtomTorsion(Molecule molecule)
    {
        Atom atom1 = molecule.contents.get(index1-1);
        Atom atom2 = molecule.contents.get(index2-1);
        Atom atom3 = molecule.contents.get(index3-1);
        Atom atom4 = molecule.contents.get(index4-1);
        return new AtomTorsion(atom1, atom2, atom3, atom4, molecule);
    }

    /**
     * Returns the hash code for this IndexTorsion
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(index1, index2, index3, index4, atomNumbersToRotate);
    }

    /**
     * Returns true if the AtomTorsions are equal.
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof IndexTorsion) )
            return false;

        IndexTorsion t = (IndexTorsion)obj;
        if ( index1 == t.index1 &&
             index2 == t.index2 &&
             index3 == t.index3 &&
             index4 == t.index4 &&
             Objects.equals(atomNumbersToRotate, t.atomNumbersToRotate) )
            return true;
        return false;
    }
}
