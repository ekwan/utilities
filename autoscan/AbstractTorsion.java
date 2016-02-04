import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class represents a torsion angle.
 */
public abstract class AbstractTorsion implements Immutable, Serializable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /**
     * Returns the torsion angle in degrees.  Atoms are explicitly defined already.
     * @return the angle in degrees
     */
    public abstract double getDihedralAngle();

    /**
     * Returns the torsion angle in degrees.  Atoms are implicitly defined by atom number,
     * so a Molecule must be provided.
     * @param molecule the Molecule associated with this torsion angle
     * @return the angle in degrees
     */
    public abstract double getDihedralAngle(Molecule molecule);

    /**
     * Returns the corresponding vector2-vector3 dihedral angle in degrees.
     * @param vector1 position of atom1
     * @param vector2 position of atom2
     * @param vector3 position of atom3
     * @param vector4 position of atom4
     * @return the dihedral angle in degrees
     */
    public static double getDihedralAngle(Vector3D vector1, Vector3D vector2, Vector3D vector3, Vector3D vector4)
    {
        Vector3D b1 = vector2.add(-1.0, vector1);
        Vector3D b2 = vector3.add(-1.0, vector2);
        Vector3D b3 = vector4.add(-1.0, vector3);

        // make sure the vectors are not collinear
        if (Vector3D.angle(b1,b2) == 0 || Vector3D.angle(b2,b3) == 0)
            {
                System.out.println("Warning! Collinear dihedral angle!");
                return 0.0;
            }

        // compute dihedral angle
        double term1 = Vector3D.dotProduct(b1.scalarMultiply( b2.getNorm() ), Vector3D.crossProduct(b2, b3) );
        double term2 = Vector3D.dotProduct(Vector3D.crossProduct(b1, b2), Vector3D.crossProduct(b2, b3) );
        return Math.toDegrees( Math.atan2(term1, term2) );
    }

    // torsions shuld be rotated through the Molecule class
}
