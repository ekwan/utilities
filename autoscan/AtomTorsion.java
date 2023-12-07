import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * Represents a torsion composed of concrete atoms.<p>
 * Torsion is: atom1-atom2-atom3-atom4<p>
 * Torsion angle calculated between atoms 2 and 3.<p>
 * If moved, the convention is to hold atom1 fixed and allow atom4 to move.<p>
 * This class contains a pre-computed list of the atoms to be rotated and
 * the rotation axis.
 */
public class AtomTorsion extends AbstractTorsion
{
    /** for serialization */
    public static final long serialVersionUID = 1L; 

    /** the Molecule corresponding to this torsion */
    public final Molecule molecule;

    /** the Atoms that comprise this torsion */
    public final Atom atom1, atom2, atom3, atom4;

    /** the Atoms including atom4 that should be rotated */
    public final List<Atom> atomsToRotate;

    /** constructor */
    public AtomTorsion(Atom atom1, Atom atom2, Atom atom3, Atom atom4, Molecule molecule)
    {
        this.atom1 = atom1;
        this.atom2 = atom2;
        this.atom3 = atom3;
        this.atom4 = atom4;
        this.molecule = molecule;
        atomsToRotate = ImmutableList.copyOf(molecule.getHalfGraph(atom2, atom3));
    }

    /** returns the corresponding torsion angle */
    public double getDihedralAngle()
    {
        return AbstractTorsion.getDihedralAngle(atom1.position, atom2.position, atom3.position, atom4.position);
    }

    /** returns the corresponding torsion angle */
    public double getDihedralAngle(Molecule molecule)
    {
        if ( molecule != this.molecule )
            throw new IllegalArgumentException("molecule mismatch in AtomTorsion!");
        return getDihedralAngle();
    }

    /** returns an IndexTorsion version of this AtomTorsion */
    public IndexTorsion getIndexTorsion()
    {
        int index1 = molecule.contents.indexOf(atom1) + 1;
        int index2 = molecule.contents.indexOf(atom2) + 1;
        int index3 = molecule.contents.indexOf(atom3) + 1;
        int index4 = molecule.contents.indexOf(atom4) + 1;

        // check indices; if 0, then there's no such Atom
        if ( index1 == 0 || index2 == 0 || index3 == 0 || index4 == 0 )
            throw new IllegalArgumentException("Atom index error in AtomTorsion!");

        List<Integer> atomNumbersToRotate = new LinkedList<>();
        for (Atom a : atomsToRotate)
            atomNumbersToRotate.add(molecule.contents.indexOf(a)+1);
        atomNumbersToRotate = ImmutableList.copyOf(atomNumbersToRotate);

        return new IndexTorsion(index1, index2, index3, index4, atomNumbersToRotate);
    }

    /**
     * Returns the hash code for this AtomTorsion
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(molecule, atom1, atom2, atom3, atom4, atomsToRotate);
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
        if ( !(obj instanceof AtomTorsion) )
            return false;

        AtomTorsion t = (AtomTorsion)obj;
        if ( Objects.equals(molecule, t.molecule) &&
             Objects.equals(atom1, t.atom1) &&
             Objects.equals(atom2, t.atom2) &&
             Objects.equals(atom3, t.atom3) &&
             Objects.equals(atom4, t.atom4) &&
             Objects.equals(atomsToRotate, t.atomsToRotate) )
            return true;
        return false;
    }
}
