import java.util.*;

public class Molecule
{
    private final String name;
    private final Map<Atom,Integer> contents;  // atom --> how many
    private final int charge;

    public static final Molecule LEUCINE_AZIRIDINE_ALDEHYDE = Molecule.parseFormula("C7H13NO", "Azl");
    public static final Molecule TERT_BUTYL_ISOCYANIDE      = Molecule.parseFormula("C5H9N", "tert-butyl isocyanide");
    public static final Molecule WATER                      = Molecule.parseFormula("H2O", "water");

    public enum Adducts
    {
        PROTON              (Molecule.parseFormula("H",     "[M+H]+",      1), true),
        PROTON_DOUBLE       (Molecule.parseFormula("H2",    "[M+2H]++",    2), true),
        PROTON_SODIUM       (Molecule.parseFormula("NaH",   "[M+H+Na]++",  2), false),
        PROTON_AMMONIUM     (Molecule.parseFormula("NH5",   "[M+H+NH4]++", 2), false),
        PROTON_METHANOL     (Molecule.parseFormula("CH5O",  "[M+MeOH+H]+", 1), false),
        SODIUM              (Molecule.parseFormula("Na",    "[M+Na]+",     1), true),
        POTASSIUM           (Molecule.parseFormula("K",     "[M+K]+",      1), false),
        PROTON_ACETONITRILE (Molecule.parseFormula("C2H4N", "[M+H+MeCN]+", 1), false),
        AMMONIUM            (Molecule.parseFormula("NH4",   "[M+NH4]+",    1), false);

        public Molecule molecule;
        public boolean common;

        private Adducts(Molecule molecule, boolean common)
        {
            this.molecule = molecule;
            this.common = common;
        }

        private static Set<Molecule> createCommonAdducts(Molecule inputMolecule)
        {
            LinkedHashSet<Molecule> returnSet = new LinkedHashSet<Molecule>();
            Map<Atom,Integer> inputContents = inputMolecule.getContents();
            int inputCharge = inputMolecule.getCharge();
            for (Adducts a : Molecule.Adducts.values())
                {
                    if (a.common==false)
                        continue;
                    Molecule m = a.molecule;
                    Map<Atom,Integer> thisContents = m.getContents();
                    int thisCharge = m.getCharge();
                    Map<Atom,Integer> newContents = addContents(inputContents, thisContents);
                    int newCharge = thisCharge + inputCharge;
                    String newName = m.getName();
                    Molecule newMolecule = new Molecule(newName, newContents, newCharge);
                    returnSet.add(newMolecule);
                }
            return returnSet;
        }

        private static Set<Molecule> createUncommonAdducts(Molecule inputMolecule)
        {
            LinkedHashSet<Molecule> returnSet = new LinkedHashSet<Molecule>();
            Map<Atom,Integer> inputContents = inputMolecule.getContents();
            int inputCharge = inputMolecule.getCharge();
            for (Adducts a : Molecule.Adducts.values())
                {
                    if (a.common==true)
                        continue;
                    Molecule m = a.molecule;
                    Map<Atom,Integer> thisContents = m.getContents();
                    int thisCharge = m.getCharge();
                    Map<Atom,Integer> newContents = addContents(inputContents, thisContents);
                    int newCharge = thisCharge + inputCharge;
                    String newName = m.getName();
                    Molecule newMolecule = new Molecule(newName, newContents, newCharge);
                    returnSet.add(newMolecule);
                }
            return returnSet;
        }

        public static String createAdductString(Molecule inputMolecule)
        {
            String returnString = "Parent Molecule:\n";
            returnString = returnString + inputMolecule.toString() + "\n\n";
            returnString = returnString + "   Common ESI Adducts:\n";
            for (Molecule m : createCommonAdducts(inputMolecule))
                returnString = returnString + m.toAdductString();
            returnString = returnString + "\n   Uncommon ESI Adducts:\n";
            for (Molecule m : createUncommonAdducts(inputMolecule))
                returnString = returnString + m.toAdductString();
            return returnString;
        }

        private static Map<Atom,Integer> addContents(Map<Atom,Integer> original, Map<Atom,Integer> toBeAdded)
        {
            LinkedHashMap<Atom,Integer> originalCopy = new LinkedHashMap<Atom,Integer>(original);
            for (Atom a : toBeAdded.keySet())
                {
                    Integer thisNumber = toBeAdded.get(a);
                    if (originalCopy.keySet().contains(a))
                        {
                            Integer oldNumber = originalCopy.get(a);
                            originalCopy.put(a, thisNumber + oldNumber);
                        }
                    else
                        originalCopy.put(a,thisNumber);
                }
            return originalCopy;
        }
    }

    private Molecule(Builder builder)
    {
        this.name = builder.name;
        this.contents = Collections.unmodifiableMap(builder.contents);
        this.charge = builder.charge;
    }

    private Molecule(String name, Map<Atom,Integer> contents, int charge)
    {
        this.name = name;
        this.contents = Collections.unmodifiableMap(contents);
        this.charge = charge;
    }

    public String getName()
    {
        return name;
    }

    public int getCharge()
    {
        return charge;
    }

    public String getMolecularFormula()
    {
        String molecularFormula = "";
        for (Atom a : contents.keySet())
            {
                molecularFormula = molecularFormula + a.getSymbol();
                if ( contents.get(a) != 1 )
                    molecularFormula = molecularFormula + contents.get(a);
            }
        return molecularFormula;
    }

    public String toString()
    {
        String returnString = "";
        String molecularFormula = getMolecularFormula();
        if ( ! molecularFormula.equals(name) )
            returnString = name + " / " + molecularFormula + " : ";
        else
            returnString = molecularFormula + " : ";
        returnString = returnString + String.format(" charge = %+d,", charge);
        returnString = returnString + String.format(" mol wt. = %.4f g/mol, exact mass = %.4f", 
                                                    getMolecularWeight(), getExactMass());
        if ( charge != 0 )
            returnString = returnString + String.format(", m/z = %.4f", getMassToChargeRatio());
        return returnString;
    }

    public String toAdductString()
    {
        return String.format("   %-20s %-20s m/z = %.4f\n", name, getMolecularFormula(), getMassToChargeRatio());
    }

    public Map<Atom,Integer> getContents()
    {
        return contents;
    }

    public Double getMolecularWeight()
    {
        Double molecularWeight = 0.0;
        for (Atom a : contents.keySet())
            molecularWeight += a.getStandardWeight() * (double)contents.get(a);
        return molecularWeight;
    }

    // returns the mass of the biggest peak in the isotopic distribution
    public Double getExactMass()
    {
        Double exactMass = 0.0;
        for (Atom a : contents.keySet())
            exactMass += a.getMostProbableMass() * (double)contents.get(a);
        return exactMass;
    }

    public Double getMassToChargeRatio()
    {
        return getExactMass() / (double)charge;
    }

    public static Molecule parseFormula(String formula)
    {
        return parseFormula(formula, formula, 0);
    }

    public static Molecule parseFormula(String formula, String name)
    {
        return parseFormula(formula, name, 0);
    }

    // only designed to handle C, H, N, O, S
    public static Molecule parseFormula(String formula, String name, int charge)
    {
        // parse characters
        ArrayList<String> fields = new ArrayList<String>();
        String currentField = "";
        String currentCharacter = "";
        String previousCharacter = "";
        for (int i=0; i < formula.length(); i++)
            {
                previousCharacter = currentCharacter;
                currentCharacter = formula.substring(i,i+1);
                
                // deal with the first character
                if ( i==0 )
                    {
                        currentField = currentCharacter;
                        continue;
                    }

                // a number can't be in the first field
                if ( fields.size() == 1 )
                    {
                        if ( fields.get(0).matches("\\d") )
                            {
                                System.out.println("Number can't be the first charcater in a molecular formula.");
                                return null;
                            }
                    }

                // uppercase-number ("C5") or number-uppercase ("5C")
                if ( ( previousCharacter.matches("[A-Z]") && currentCharacter.matches("\\d")   ) ||
                          ( previousCharacter.matches("\\d")   && currentCharacter.matches("[A-Z]") )    )
                    {
                        fields.add(currentField);
                        currentField = currentCharacter;
                        if ( i == formula.length() - 1 )
                            fields.add(currentField);
                    }
                
                // uppercase-uppercase ("NO") or lowercase-uppercase ("FeN")
                else if ( ( previousCharacter.matches("[A-Z]") && currentCharacter.matches("[A-Z]")   ) ||
                          ( previousCharacter.matches("[a-z]") && currentCharacter.matches("[A-Z]") )    )
                    {
                        fields.add(currentField);
                        fields.add("1");
                        currentField = currentCharacter;
                        if ( i == formula.length() - 1 )
                            fields.add(currentField);
                    }

                // uppercase-lowercase ("Fe") or number-number ("20")
                else if ( ( previousCharacter.matches("[A-Z]") && currentCharacter.matches("[a-z]")   ) ||
                          ( previousCharacter.matches("\\d")   && currentCharacter.matches("\\d") )    )
                    {
                        currentField = currentField + currentCharacter;
                        if ( i == formula.length() - 1 )
                            fields.add(currentField);
                    }
                
                // unexpected
                else
                    {
                        System.out.println("Unexpected formula!");
                        return null;
                    }
            }

        if ( fields.size() == 0 )
            fields.add(currentField);
        
        if ( fields.size() == 1 )
            fields.add("1");
        else if ( fields.get(fields.size()-1).matches("[A-Za-z]") )
            fields.add("1");

        // debug only
        //System.out.println(formula);
        //for ( String s : fields )
        //    System.out.println(s);

        // check for an even number of fields
        if (fields.size() % 2 != 0)
            {
                System.out.println("Not an even number of fields (" + fields.size() + ")!");
                return null;
            }

        // parse fields
        Molecule.Builder builder = new Molecule.Builder(name);
        for (int i=0; i < fields.size() - 1; i += 2)
            {
                String element = fields.get(i);
                Integer numberOfAtoms = Integer.valueOf(fields.get(i+1));
                Atom thisAtom = null;
                for (Atom a : Atom.KNOWN_ELEMENTS)
                    {
                        if (element.equals(a.getSymbol()))
                            thisAtom = a;
                    }
                if ( thisAtom == null )
                    {
                        System.out.println("Did not recognize element " + element + ".");
                        return null;
                    }
                builder = builder.addAtom(thisAtom, numberOfAtoms);
            }
        return builder.setCharge(charge).build();
    }

    public static class Builder
    {
        private final String name;
        private LinkedHashMap<Atom,Integer> contents;
        private int charge;

        public Builder(String name)
            {
                this.name = name;
                contents = new LinkedHashMap<Atom,Integer>();
                charge = 0;
            }

        public Builder addAtom(Atom atom, Integer howMany)
            {
                contents.put(atom,howMany);
                return this;
            }

        public Builder setCharge(int charge)
            {
                this.charge = charge;
                return this;
            }

        public Molecule build()
            {
                if ( contents.keySet().size() == 0 || name == null )
                    throw new IllegalStateException("error building molecule");
                return new Molecule(this);
            }
    }
}
