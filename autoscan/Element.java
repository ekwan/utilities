/**
 * This class represents an element.
 */
public enum Element
{
    CARBON  ("C", 3.75, 0.105),
    HYDROGEN("H", 2.50, 0.03),
    NITROGEN("N", 3.25, 0.17),
    OXYGEN  ("O", 2.96, 0.21),
    SULFUR  ("S", 3.60, 0.355),
    FLUORINE ("F", 2.94, 0.061),
    CHLORINE ("Cl", 3.40, 0.300), 
    PHOSPHORUS ("P", 3.72, 0.20),
    BROMINE ("Br", 3.15, 0.15),
    TITANIUM ("Ti", 3.75, 0.105),
    DUMMY ("Q", 0, 0);

    /** the atomic symbol for this element */
    public String symbol;

    /** rough vdW parameter sigma from OPLS */
    public double sigma;
    
    /** rough vdW paramater epsilon from OPLS */
    public double epsilon;

    /** Constructor. */
    Element(String symbol, double sigma, double epsilon)
    {
        this.symbol = symbol;
        this.sigma = sigma;
        this.epsilon = epsilon;
    }

    /**
     * Returns the Element corresponding to the inputted String.
     * @param symbol the symbol for the requested element
     * @return the corresponding Element enum element
     */
    public static Element getElement(String symbol)
    {
        if ( symbol.equals("C") || symbol.equals("CA") )
            return CARBON;
        else if ( symbol.equals("H") || symbol.equals("HN") || symbol.equals("HS") || symbol.equals("HO") )
            return HYDROGEN;
        else if ( symbol.equals("N") )
            return NITROGEN;
        else if ( symbol.equals("O") || symbol.equals("O-") || symbol.equals("OH") )
            return OXYGEN;
        else if ( symbol.equals("S") || symbol.equals("SH") || symbol.equals("SS") )
            return SULFUR;
        else if ( symbol.equals("P") )
            return PHOSPHORUS;
        else if ( symbol.equals("F") )
            return FLUORINE;
        else if ( symbol.equals("Cl") )
            return CHLORINE;
        else if ( symbol.equals("Q") )
            return DUMMY;
        else if ( symbol.equals("Br") )
            return BROMINE;
        else if ( symbol.equals("Ti") )
            return TITANIUM;
        else
            throw new IllegalArgumentException("Unrecognized atom symbol (" + symbol + ")!");
    }

    /**
     * Returns a string representation of this Element.
     * @return the String
     */
    public String toString()
    {
        return String.format("%s (%s, sigma=%.3f, epsilon=%.3f)", name(), symbol, sigma, epsilon);
    }
}
