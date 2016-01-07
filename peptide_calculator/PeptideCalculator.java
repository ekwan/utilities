import java.util.*;
import java.io.*;

public class PeptideCalculator
{
    public static final String PEPTIDE_FILENAME = "peptideHistory.txt";

    public static void main (String[] args)
    {
        // print banner
        System.out.println("====================");
        System.out.println(" Peptide Calculator ");
        System.out.println("                    ");
        System.out.println(" Eugene Kwan        ");
        System.out.println(" November 2013      ");
        System.out.println(" Harvard University ");
        System.out.println("====================\n");

        // read peptide database
        File peptideFile = new File(PEPTIDE_FILENAME);
        System.out.print("Reading peptide database...");
        ArrayList<String> peptideDatabase = new ArrayList<String>();
        if ( !peptideFile.isFile() )
            System.out.println("warning: peptide database file not found!");
        else
            {
                try ( Scanner scanner = new Scanner(new FileReader(PEPTIDE_FILENAME)) )
                    {
                        while (scanner.hasNextLine())
                            {
                                String currentLine = scanner.nextLine();
                                if ( !currentLine.contains("Azl") )
                                    continue;
                                String modifiedLine = "Pro - " + currentLine.substring(11).replaceAll("_D|_L|_S|_R","");
                                peptideDatabase.add(modifiedLine);
                                //System.out.println("peptide # " + (peptideDatabase.size()-1) + ": " + modifiedLine);
                            }
                        System.out.println(peptideDatabase.size() + " peptide sequences have been read.");
                    }
                catch (IOException e)
                    {
                        System.out.println("error while reading file.");
                    }
            }

        // main loop
        while (true)
            {
                // ask for input
                Scanner console = new Scanner(System.in);
                System.out.println("Please enter a peptide sequence (e.g. Pro-Ala-Glu):");
                System.out.print("> ");
                String inputString = console.nextLine();

                // make peptide
                Peptide thisPeptide = null;
                try
                    {
                        // if input string contains a number, try to parse it
                        if (inputString.matches(".*\\d.*"))
                            {
                                try
                                    {
                                        int peptideNumber = Integer.parseInt(inputString);
                                        if ( peptideNumber <= peptideDatabase.size() + 1 && peptideNumber >= 0 )
                                            {
                                                inputString = peptideDatabase.get(peptideNumber);
                                                System.out.println("Peptide number " + peptideNumber + " is " + inputString + ".\n");
                                            }
                                        else
                                            {
                                                System.out.println("Requested peptide number out of range.\n");
                                                continue;
                                            }
                                    }
                                catch (NumberFormatException e)
                                    {
                                        System.out.println("Please enter a valid peptide sequence or numeric identifier.\n");
                                        continue;
                                    }
                            }

                        thisPeptide = Peptide.makeSequence(inputString);
                    }
                catch (Exception e)
                    {
                        System.out.println("Invalid peptide.\n");
                        continue;
                    }

                // print out result
                System.out.println("\n" + getPeptideString(thisPeptide));

                // ask whether the user wants more peptides
                System.out.print("\nCalculate another peptide? ");
                inputString = console.nextLine();
                if ( ! (inputString.equalsIgnoreCase("y") || inputString.equalsIgnoreCase("yes") ) )
                    break;
                else
                    System.out.println("\n\n\n");
            }
    }

    public static String getPeptideString(Peptide inputPeptide)
    {
        Molecule linearPeptide             = inputPeptide.getLinearPeptide();
        String linearPeptideString         = Molecule.Adducts.createAdductString(linearPeptide);
        
        Molecule protectedMacrocycle       = inputPeptide.getProtectedMacrocycle();
        String protectedMacrocycleString   = Molecule.Adducts.createAdductString(protectedMacrocycle);
        
        Molecule deprotectedMacrocycle     = inputPeptide.getDeprotectedMacrocycle();
        String deprotectedMacrocycleString = Molecule.Adducts.createAdductString(deprotectedMacrocycle);

        String returnString                = "=== Linear Peptide ===\n" + linearPeptideString + "\n";
        returnString                       = returnString + "=== Protected Macrocycle ===\n" + protectedMacrocycleString + "\n";
        returnString                       = returnString + "=== Deprotected Macrocycle ===\n" + deprotectedMacrocycleString;

        return returnString;
    }
}
