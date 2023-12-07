import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 *  Represents an atom.  This class is immutable.
 */
public class Atom implements Immutable, Serializable, Comparable<Atom>
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** the Element enum corresponding to this Atom */
    public final Element element;

    /** the location of this Atom */
    public final Vector3D position;

    /** atom type in TINKER for this Atom */
    public final int tinkerAtomType;

    /** constructs a new Atom */
    public Atom(String symbol, Vector3D position, int tinkerAtomType)
    {
        this.element        = Element.getElement(symbol);
        this.position       = position;
        if ( tinkerAtomType < 0 )
            throw new IllegalArgumentException("Tinker atom type must be positive!");
        this.tinkerAtomType = tinkerAtomType;
    }

    /** constructs a new Atom */
    public Atom(Element element, Vector3D position, int tinkerAtomType)
    {
        this.element = element;
        this.position = position;
        if ( tinkerAtomType < 0 )
            throw new IllegalArgumentException("Tinker atom type must be positive!");
        this.tinkerAtomType = tinkerAtomType;
    }

    public int compareTo(Atom a)
    {
        double x1 =   position.getX();
        double x2 = a.position.getX();
        double y1 =   position.getY();
        double y2 = a.position.getY();
        double z1 =   position.getZ();
        double z2 = a.position.getZ();
        return ComparisonChain.start().compare(x1, x2).compare(y1, y2).compare(z1, z2).result();
    }

    /**
     * Returns a copy of this atom with a new atom type.
     */
    public Atom setAtomType(int atomType)
    {
        return new Atom(element, position, atomType);
    }

    /**
     * Returns a copy of this atom with a new position.
     */
    public Atom moveAtom(Vector3D newPosition)
    {
        return new Atom(element.symbol, newPosition, tinkerAtomType);
    }

    /**
     * Factory method to create a new Atom by transforming this one.
     * The rotation is applied before the translation.  
     * @param rot a three-dimensional rotation.
     * @param shift a vector that we add to the position of this after it has been rotated.
     */
    public Atom transform(Rotation rot, Vector3D shift)
    {
        return new Atom(element.symbol, rot.applyTo(position).add(shift), tinkerAtomType);
    }

    /** convenience method for moving an atom based on a map of old atoms and new atoms */
    public Atom moveAtom(Map<Atom,Atom> atomMap)
    {
        if (atomMap.containsKey(this))
            return atomMap.get(this);
        else
            return this;
    }

    /**
     * Return a string representation of this Atom.
     * @return the string
     */
    @Override
    public String toString()
    {
        return String.format("%-2s %10.6f %10.6f %10.6f", element.symbol, position.getX(), position.getY(), position.getZ());
        //return String.format("%-2s (%3d) %10.6f %10.6f %10.6f", element.symbol, tinkerAtomType, position.getX(), position.getY(), position.getZ());
    }

    /**
     * Returns the hash code.
     * @return the hash code
     */
    public int hashCode()
    {
        return Objects.hash(element, position, tinkerAtomType);
    }

    // tolerance can be added, but will destroy consistency with hash code
    // public static final double TOLERANCE = 0.000000;

    /** 
     * Tests for object equality
     * @param obj another object
     * @return true if equal
     */
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        else if ( obj == this )
            return true;
        else if ( !(obj instanceof Atom) )
            return false;

        Atom anotherAtom = (Atom)obj;
        if ( element == anotherAtom.element &&
             tinkerAtomType == anotherAtom.tinkerAtomType )
            {
              /*  if ( Math.abs(position.getX() - anotherAtom.position.getX()) < TOLERANCE &&
                     Math.abs(position.getY() - anotherAtom.position.getY()) < TOLERANCE &&
                     Math.abs(position.getZ() - anotherAtom.position.getZ()) < TOLERANCE    )
             */
            if ( Objects.equals(position, anotherAtom.position) )
                 return true;
            }
        return false;
    }

    /** for testing */
    public static void main(String[] args)
    {
        Atom atom1 = new Atom("C", new Vector3D(1.0, 2.0, 3.0), 103);
        Atom atom2 = new Atom("C", new Vector3D(2.0, 2.0, 3.0), 103);
        Atom atom3 = new Atom("C", new Vector3D(3.0, 3.0, 3.0), 103);
        Atom atom4 = new Atom("C", new Vector3D(3.0, 2.0, 3.0), 103);
        Atom atom5 = new Atom("C", new Vector3D(3.0, 2.0, 4.0), 103);
        List<Atom> list = new LinkedList<>();
        list.add(atom2);
        list.add(atom1);
        list.add(atom3);
        list.add(atom4);
        list.add(atom5);
        Collections.shuffle(list);
        Collections.sort(list);
        for (Atom a : list)
            System.out.println(a);
    }
} // end of class Atom
