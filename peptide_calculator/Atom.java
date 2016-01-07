import java.util.*;

public class Atom
{
    private final String symbol;
    private final Map<Isotope,Double> abundances;  // isotope --> relative abundances
    private final Double standardWeight;

    // data taken from NIST
    // http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    public static final Atom HYDROGEN  = new Atom.Builder("H", 1.00794).
                                         addIsotope(1, 1.00782503207, 0.999885).
                                         addIsotope(2, 2.0141017778, 0.000115).
                                         build();

    public static final Atom CARBON    = new Atom.Builder("C", 12.0107).
                                         addIsotope(12, 12.0, 0.9893).
                                         addIsotope(13, 13.0033548378, 0.0107).
                                         build();

    public static final Atom NITROGEN  = new Atom.Builder("N", 14.0067).
                                         addIsotope(14, 14.0030740048, 0.99636).
                                         addIsotope(15, 15.0001088982, 0.00364).
                                         build();

    public static final Atom OXYGEN    = new Atom.Builder("O", 15.9994).
                                         addIsotope(16, 15.99491461956, 0.99757).
                                         addIsotope(17, 16.99913170, 0.00038).
                                         addIsotope(18, 17.9991610, 0.00205).
                                         build();

    public static final Atom SULFUR    = new Atom.Builder("S", 32.065).
                                         addIsotope(32, 31.97207100, 0.9499).
                                         addIsotope(33, 32.97145876, 0.0075).
                                         addIsotope(34, 33.96786690, 0.0425).
                                         addIsotope(36, 35.96708076, 0.0001).
                                         build();

    public static final Atom SODIUM    = new Atom.Builder("Na", 22.9897692809).
                                         addIsotope(23, 22.9897692809, 1.0000).
                                         build();

    public static final Atom POTASSIUM = new Atom.Builder("K", 39.0983).
                                         addIsotope(39, 38.96370668, 0.932581).
                                         addIsotope(40, 39.96399848, 0.000117).
                                         addIsotope(41, 40.96182576, 0.067302).
                                         build();

    public static final Set<Atom> KNOWN_ELEMENTS;

    static
    {
        HashSet<Atom> tempSet = new HashSet<Atom>();
        tempSet.add(CARBON);
        tempSet.add(HYDROGEN);
        tempSet.add(NITROGEN);
        tempSet.add(OXYGEN);
        tempSet.add(SULFUR);
        tempSet.add(SODIUM);
        tempSet.add(POTASSIUM);
        KNOWN_ELEMENTS = Collections.unmodifiableSet(tempSet);
    }

    private Atom(Builder builder)
    {
        symbol = builder.symbol;
        abundances = Collections.unmodifiableMap(builder.abundances);
        standardWeight = builder.standardWeight;
    }

    public String toString()
    {
        return symbol + "\n" + abundances.toString();
    }

    public String getSymbol()
    {
        return symbol;
    }

    public Double getMostProbableMass()
    {
        Double mostProbableMass = null;
        Double biggestAbundance = null;
        for (Isotope i : abundances.keySet())
            {
                Double thisMass = i.getExactMass();
                Double thisAbundance = abundances.get(i);
                if ( biggestAbundance == null || thisAbundance > biggestAbundance )
                    {
                        mostProbableMass = thisMass;
                        biggestAbundance = thisAbundance;
                    }
            }
        return mostProbableMass;
    }

    public Map<Isotope,Double> getAbundances()
    {
        return abundances;
    }

    public Double getStandardWeight()
    {
        return standardWeight;
    }

    public static class Builder
    {
        private final String symbol;
        private HashMap<Isotope,Double> abundances;
        private final Double standardWeight;

        public Builder(String symbol, Double standardWeight)
        {
            this.symbol = symbol;
            abundances = new HashMap<Isotope,Double>();
            this.standardWeight = standardWeight;
        }

        public Builder addIsotope(int nominalMass, Double exactMass, Double abundance)
        {
            abundances.put(new Isotope(nominalMass, exactMass), abundance);
            return this;
        }

        public Atom build()
        {
            Double abundanceSum = 0.0;
            for ( Isotope i : abundances.keySet() )
                abundanceSum += abundances.get(i);
            if ( abundances.size() == 0 || standardWeight == null ||
                 abundanceSum > 1.001 || abundanceSum < 0.999 || symbol == null )
                throw new IllegalStateException("error in building Atom");
            return new Atom(this);
        }
    }

    private static class Isotope
    {
        private int nominalMass;
        private Double exactMass;

        private Isotope(int nominalMass, Double exactMass)
        {
            this.nominalMass = nominalMass;
            this.exactMass = exactMass;
        }

        public int getNominalMass()
        {
            return nominalMass;
        }

        public Double getExactMass()
        {
            return exactMass;
        }

        public String toString()
        {
            return String.format("(%d, %.4f)", nominalMass, exactMass);
        }
    }
}
