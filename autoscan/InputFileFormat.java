import java.io.*;
import java.util.*;

/**
 * Represents communication from Mandor to the outside world.
 * This class is immutable.
 */
public abstract class InputFileFormat implements FileFormat, Immutable
{
    /** The string that will be written to a file */
    public final String stringRepresentation;

    public InputFileFormat(String stringRepresentation)
    {
        this.stringRepresentation = stringRepresentation;
    }

    /**
     * Writes the file to disk.
     */
    public void write(String filename)
    {
        writeStringToDisk(stringRepresentation, filename);
	}

    /**
     * Convenience method that writes a string to a file.
     */
    public static void writeStringToDisk(String string, String filename)
    {
        try (PrintWriter outputFile = new PrintWriter(filename))
            {
                outputFile.print(string);
            }
        catch (IOException e)
            {
                // abort program if there's a problem
                System.out.println("Error writing to " + filename + "!");
                e.printStackTrace();
            }
    }

    /**
     * Convenience method that appends a string to a file.
     */
    public static void appendStringToDisk(String string, String filename)
    {
        try
            {
                File file = new File(filename);
                if ( ! file.exists() )
                    file.createNewFile();
                FileWriter fileWriter = new FileWriter(filename,true);
                BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
                bufferedWriter.write(string);
                bufferedWriter.close();
            }
        catch (IOException e)
            {
                System.out.println("Error appending to " + filename + "!");
                e.printStackTrace();
            }
    }

    public int hashCode()
    {
        return Objects.hash(stringRepresentation);
    }

    public String toString()
    {
        return stringRepresentation;
    }
}
