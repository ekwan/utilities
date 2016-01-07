import java.util.*;

public class Peptide
{
    private Molecule linearPeptide;
    private Molecule protectedMacrocycle;
    private Molecule deprotectedMacrocycle;

    private Peptide(Molecule linearPeptide, Molecule protectedMacrocycle, Molecule deprotectedMacrocycle)
    {
        this.linearPeptide = linearPeptide;
        this.protectedMacrocycle = protectedMacrocycle;
        this.deprotectedMacrocycle = deprotectedMacrocycle;
    }

    public Molecule getLinearPeptide()
    {
        return linearPeptide;
    }

    public Molecule getProtectedMacrocycle()
    {
        return protectedMacrocycle;
    }

    public Molecule getDeprotectedMacrocycle()
    {
        return deprotectedMacrocycle;
    }

    public static Peptide makeSequence(String sequence)
    {
        if ( sequence.trim().length() == 0 )
            throw new IllegalArgumentException("Empty sequence.");

        // parse amino acid string
        ArrayList<AminoAcid> AAsequence = new ArrayList<AminoAcid>();
        for (String s : sequence.split("-"))
            {
                String AAstring = s.trim();
                AminoAcid thisAminoAcid = null;
                for ( AminoAcid a : AminoAcid.values() )
                    {
                        String longName = a.getLongName();
                        String shortName = a.getShortName();
                        if ( AAstring.equalsIgnoreCase(longName) || AAstring.equalsIgnoreCase(shortName) )
                            {
                                thisAminoAcid = a;
                                break;
                            }
                    }
                if ( thisAminoAcid == null )
                    throw new IllegalArgumentException("Invalid amino acid found: " + AAstring);
                AAsequence.add(thisAminoAcid);
            }

        if ( AAsequence.size() < 2 )
            throw new IllegalArgumentException("Peptide must have at least two residues.");

        // construct linear peptide
        Map<Atom,Integer> linearPeptideContents = new LinkedHashMap<Atom,Integer>();
        String linearPeptideLongName = "";
        String linearPeptideShortName = "";
        for (int i=0; i < AAsequence.size(); i++)
            {
                AminoAcid thisAminoAcid = AAsequence.get(i);
                Molecule  thisMolecule  = thisAminoAcid.getProtectedMolecule();
                Map<Atom,Integer> thisContents = thisMolecule.getContents();
                
                // add molecular formula to current one
                linearPeptideContents = addContents(linearPeptideContents,thisContents);

                // remove water, but only for residues that aren't the first residue
                if ( i > 0 )
                    linearPeptideContents = subtractContents(linearPeptideContents,Molecule.WATER);
                
                // add to name
                String thisLongName = thisAminoAcid.getLongName();
                String thisShortName = thisAminoAcid.getShortName();
                if ( linearPeptideLongName.length() > 0 )
                    linearPeptideLongName = linearPeptideLongName + "-" + thisLongName;
                else
                    linearPeptideLongName = thisLongName;
                linearPeptideShortName = linearPeptideShortName + thisShortName;
            }
        String linearPeptideName = linearPeptideLongName + " (" + linearPeptideShortName + ")";
        Molecule.Builder linearPeptideBuilder = new Molecule.Builder(linearPeptideName);
        for (Atom a : linearPeptideContents.keySet())
            linearPeptideBuilder.addAtom(a, linearPeptideContents.get(a));
        Molecule linearPeptide = linearPeptideBuilder.build();

        // construct protected macrocycle
        Map<Atom,Integer> protectedMacrocycleContents = new LinkedHashMap<Atom,Integer>();
        String protectedMacrocycleLongName = Molecule.LEUCINE_AZIRIDINE_ALDEHYDE.getName();
        String protectedMacrocycleShortName = "";

        for (int i=0; i < AAsequence.size(); i++)
            {
                AminoAcid thisAminoAcid = AAsequence.get(i);
                Molecule  thisMolecule  = thisAminoAcid.getProtectedMolecule();
                Map<Atom,Integer> thisContents = thisMolecule.getContents();
                
                // add molecular formula to current one
                protectedMacrocycleContents = addContents(protectedMacrocycleContents,thisContents);

                // remove water, but only for residues that aren't the first residue
                if ( i > 0 )
                    protectedMacrocycleContents = subtractContents(protectedMacrocycleContents,Molecule.WATER);
                
                // add to name
                String thisLongName = thisAminoAcid.getLongName();
                String thisShortName = thisAminoAcid.getShortName();
                if ( protectedMacrocycleLongName.length() > 0 )
                    protectedMacrocycleLongName = protectedMacrocycleLongName + "-" + thisLongName;
                else
                    protectedMacrocycleLongName = thisLongName;
                protectedMacrocycleShortName = protectedMacrocycleShortName + thisShortName;
            
                // macrocyclize at end
                if ( i== AAsequence.size() - 1 )
                    {
                        protectedMacrocycleContents = addContents(protectedMacrocycleContents,Molecule.LEUCINE_AZIRIDINE_ALDEHYDE);
                        protectedMacrocycleContents = addContents(protectedMacrocycleContents,Molecule.TERT_BUTYL_ISOCYANIDE);
                        protectedMacrocycleContents = subtractContents(protectedMacrocycleContents,Molecule.WATER);
                    }
            }
            
        String protectedMacrocycleName = protectedMacrocycleLongName + " (" + protectedMacrocycleShortName + ")";
        Molecule.Builder protectedMacrocycleBuilder = new Molecule.Builder(protectedMacrocycleName);
        for (Atom a : protectedMacrocycleContents.keySet())
            protectedMacrocycleBuilder.addAtom(a, protectedMacrocycleContents.get(a));
        Molecule protectedMacrocycle = protectedMacrocycleBuilder.build();

        // construct deprotected macrocycle
        Map<Atom,Integer> deprotectedMacrocycleContents = new LinkedHashMap<Atom,Integer>();
        String deprotectedMacrocycleLongName = Molecule.LEUCINE_AZIRIDINE_ALDEHYDE.getName();
        String deprotectedMacrocycleShortName = "";

        for (int i=0; i < AAsequence.size(); i++)
            {
                AminoAcid thisAminoAcid = AAsequence.get(i);
                Molecule  thisMolecule  = thisAminoAcid.getDeprotectedMolecule();
                Map<Atom,Integer> thisContents = thisMolecule.getContents();
                
                // add molecular formula to current one
                deprotectedMacrocycleContents = addContents(deprotectedMacrocycleContents,thisContents);

                // remove water, but only for residues that aren't the first residue
                if ( i > 0 )
                    deprotectedMacrocycleContents = subtractContents(deprotectedMacrocycleContents,Molecule.WATER);
                
                // add to name
                String thisLongName = thisAminoAcid.getLongName();
                String thisShortName = thisAminoAcid.getShortName();
                if ( deprotectedMacrocycleLongName.length() > 0 )
                    deprotectedMacrocycleLongName = deprotectedMacrocycleLongName + "-" + thisLongName;
                else
                    deprotectedMacrocycleLongName = thisLongName;
                deprotectedMacrocycleShortName = deprotectedMacrocycleShortName + thisShortName;
            
                // macrocyclize at end
                if ( i== AAsequence.size() - 1 )
                    {
                        deprotectedMacrocycleContents = addContents(deprotectedMacrocycleContents,Molecule.LEUCINE_AZIRIDINE_ALDEHYDE);
                        deprotectedMacrocycleContents = addContents(deprotectedMacrocycleContents,Molecule.TERT_BUTYL_ISOCYANIDE);
                        deprotectedMacrocycleContents = subtractContents(deprotectedMacrocycleContents,Molecule.WATER);
                    }
            }
            
        String deprotectedMacrocycleName = deprotectedMacrocycleLongName + " (" + deprotectedMacrocycleShortName + ")";
        Molecule.Builder deprotectedMacrocycleBuilder = new Molecule.Builder(deprotectedMacrocycleName);
        for (Atom a : deprotectedMacrocycleContents.keySet())
            deprotectedMacrocycleBuilder.addAtom(a, deprotectedMacrocycleContents.get(a));
        Molecule deprotectedMacrocycle = deprotectedMacrocycleBuilder.build();

        // return object
        return new Peptide(linearPeptide, protectedMacrocycle, deprotectedMacrocycle);
    }

    private static Map<Atom,Integer> addContents(Map<Atom,Integer> original, Molecule molecule)
    {
        return addContents(original, molecule.getContents());
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

    private static Map<Atom,Integer> subtractContents(Map<Atom,Integer> original, Molecule molecule)
    {
        return subtractContents(original, molecule.getContents());
    }

    private static Map<Atom,Integer> subtractContents(Map<Atom,Integer> original, Map<Atom,Integer> toBeSubtracted)
    {
        LinkedHashMap<Atom,Integer> originalCopy = new LinkedHashMap<Atom,Integer>(original);
        for (Atom a : toBeSubtracted.keySet())
            {
                Integer thisNumber = toBeSubtracted.get(a);
                if (originalCopy.keySet().contains(a))
                    {
                        Integer oldNumber = originalCopy.get(a);
                        originalCopy.put(a, oldNumber - thisNumber);
                        if ( oldNumber - thisNumber < 1 )
                            System.out.println("Warning, negative molecular formula!");
                    }
                else
                    originalCopy.put(a,thisNumber);
            }
        return originalCopy;
    }
}
