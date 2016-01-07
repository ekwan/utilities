public enum AminoAcid
{
    // amino acid name(three letter abbreviation, one letter abbreviation, protecting group name (abbreviated)
    //                 protected mass, not protected mass, amino acid loading)
    ARGININE      ("Arg", "R", "Pbf",  "C6H14N4O2",  "C19H30N4O5S", 0.50),
    HISTIDINE     ("His", "H", "Trt",  "C6H9N3O2",   "C25H23N3O2",  0.50),
    ASPARTIC_ACID ("Asp", "D", "OtBu", "C4H7NO4",    "C8H15NO4",    0.50),
    GLUTAMIC_ACID ("Glu", "E", "OtBu", "C5H9NO4",    "C9H17NO4",    0.50),
    SERINE        ("Ser", "S", "tBu",  "C3H7NO3",    "C7H15NO3",    0.50),
    THREONINE     ("Thr", "T", "tBu",  "C4H9NO3",    "C8H17NO3",    0.50),
    ASPARAGINE    ("Asn", "N", "Trt",  "C4H8N2O3",   "C23H22N2O3",  0.50),
    GLUTAMINE     ("Gln", "Q", "Trt",  "C5H10N2O3",  "C24H24N2O3",  0.50),
    GLYCINE       ("Gly", "G", "",     "C2H5NO2",    "C2H5NO2",     0.50),
    PROLINE       ("Pro", "P", "",     "C5H9NO2",    "C5H9NO2",     0.50),
    ALANINE       ("Ala", "A", "",     "C3H7NO2",    "C3H7NO2",     0.50),
    VALINE        ("Val", "V", "",     "C5H11NO2",   "C5H11NO2",    0.50),
    ISOLEUCINE    ("Ile", "I", "",     "C6H13NO2",   "C6H13NO2",    0.50),
    LEUCINE       ("Leu", "L", "",     "C6H13NO2",   "C6H13NO2",    0.50),
    PHENYLALANINE ("Phe", "F", "",     "C9H11NO2",   "C9H11NO2",    0.50),
    TYROSINE      ("Tyr", "Y", "tBu",  "C9H11NO3",   "C13H19NO3",   0.50),
    TRYPTOPHAN    ("Trp", "W", "Boc",  "C11H12N2O2", "C16H20N2O4",  0.50);

    private String longName;
    private String shortName;
    private String protectingGroupName;
    private Molecule deprotectedMolecule;
    private Molecule protectedMolecule;
    private Double resinLoading;           // mmol per gram

    private AminoAcid(String longName, String shortName, String protectingGroupName, String deprotectedFormula, String protectedFormula, Double resinLoading)
    {
        this.longName = longName;
        this.shortName = shortName;
        this.protectingGroupName = protectingGroupName;

        String protectedMoleculeName = longName;
        if ( protectingGroupName.length() > 0 )
            protectedMoleculeName = protectedMoleculeName + "(" + protectingGroupName + ")";

        this.deprotectedMolecule = Molecule.parseFormula(deprotectedFormula, longName); 
        this.protectedMolecule = Molecule.parseFormula(protectedFormula, protectedMoleculeName); 
        this.resinLoading = resinLoading;
    }

    public String toString()
    {
        String returnString = "";
        if ( ! protectingGroupName.equals("") )
            returnString = String.format("%s(%s) (%s) : %s\n", longName, protectingGroupName, shortName, protectedMolecule.toString());
        returnString = returnString + String.format("%s (%s) : %s", longName, shortName, deprotectedMolecule.toString());
        return returnString;
    }

    public String getLongName()
    {
        return longName;
    }

    public String getShortName()
    {
        return shortName;
    }

    public String getProtectingGroupName()
    {
        return protectingGroupName;
    }

    public Molecule getProtectedMolecule()
    {
        return protectedMolecule;
    }

    public Molecule getDeprotectedMolecule()
    {
        return deprotectedMolecule;
    }

    public Double getResinLoading()
    {
        return resinLoading;
    }
}
