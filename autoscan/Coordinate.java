import java.util.*;

public class Coordinate implements Immutable
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

