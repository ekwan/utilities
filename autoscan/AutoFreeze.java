import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;

/**
 * Reads in a series of output files.  For each output file, a series of
 * frozen bond jobs is generated.  Reads the bonding connectivity of
 * all files from one template file.
 */
public class AutoFreeze
{
    /**
     * Changes the connectivity in inputMolecule to match that of templateMolecule.
     * @param templateMolecule the bonds will be read from this Molecule
     * @param inputMolecule the bonds in this molecule will be replaced by those from templateMolecule
     * @return the molecule with the new connectivity graph
     */
    public static Molecule copyConnectivity(Molecule templateMolecule, Molecule inputMolecule) {
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> oldConnectivity = templateMolecule.connectivity;
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
        List<Atom> newContents = new ArrayList<>(inputMolecule.contents);
        for (Atom a : newContents)
            newConnectivity.addVertex(a);
        for (DefaultWeightedEdge e : oldConnectivity.edgeSet())
            {
                // get old edge data
                Double bondOrder      = oldConnectivity.getEdgeWeight(e);
                Atom templateFromAtom = oldConnectivity.getEdgeSource(e);
                Atom templateToAtom   = oldConnectivity.getEdgeTarget(e);
				int templateFromIndex = templateMolecule.contents.indexOf(templateFromAtom);
                int templateToIndex   = templateMolecule.contents.indexOf(templateToAtom);
                if ( templateFromIndex < 0 || templateToIndex < 0 )
                    throw new IllegalArgumentException("check indices");
                Atom inputFromAtom    = newContents.get(templateFromIndex);
                Atom inputToAtom      = newContents.get(templateToIndex);

                // create new edge
                if (! newContents.contains(inputFromAtom) || ! newContents.contains(inputToAtom))
                    System.out.println("FromAtom: " + inputMolecule.getAtomString(inputFromAtom) +  " ToAtom: " + inputMolecule.getAtomString(inputToAtom));
                DefaultWeightedEdge newEdge = newConnectivity.addEdge(inputFromAtom,inputToAtom);
                newConnectivity.setEdgeWeight(newEdge, bondOrder);
            }
        return new Molecule(inputMolecule.name, newContents, newConnectivity, inputMolecule.energy);
    }

    public static void main(String[] args)
        {
            // read connectivity from template file
            GaussianOutputFile connectivityTemplateFile = new GaussianOutputFile("clay_new-glu_gluC6-AADL-Na_1Me2O-phenylate-ts-b3lyp_d3bj-631gd-pcm_thf.out");
            Molecule connectivityTemplateMolecule = connectivityTemplateFile.molecule;
            
            // add the forming bond (atom numbers)
            // ensure there are no loops
            connectivityTemplateMolecule = connectivityTemplateMolecule.addBond(34,19);
            connectivityTemplateMolecule = connectivityTemplateMolecule.removeBond(11,37);
            connectivityTemplateMolecule = connectivityTemplateMolecule.removeBond(6,43);
            /*GaussianInputFile g = new GaussianInputFile(connectivityTemplateMolecule, "title", "#p geom=connect", "", true, 0, 1);
            g.write("temp.gjf");
            
            Atom a1 = connectivityTemplateMolecule.contents.get(33);
            Atom a2 = connectivityTemplateMolecule.contents.get(18);
            Set<Atom> set = connectivityTemplateMolecule.getHalfGraph(a1, a2);
            for (Atom a : set)
                System.out.println(connectivityTemplateMolecule.getAtomString(a));
            System.exit(1);*/

            // parameters for new input files
            String keywords = "%mem=12GB\n%nprocshared=8\n#p b3lyp 6-31g* scrf=(pcm,solvent=tetrahydrofuran) pop=none opt=modredundant freq=noraman";
            String tail = "B 19 34 F\n";
 
            // read each file
			File directory = new File("/Users/ekwan/research/clay/templates/");
			File[] files = directory.listFiles((d, name) -> name.endsWith(".out"));
            if ( files == null )
                throw new NullPointerException("null files");
            for (File f : files) {
                // get filename
                String filename = f.getName();
                System.out.println("> " + filename);

                // copy connectivity
                Molecule inputMolecule = new GaussianOutputFile(f.getAbsolutePath()).molecule;
                Molecule newMoleculeTemplate = copyConnectivity(connectivityTemplateMolecule, inputMolecule); 

                // create input files
				for (double bondDistance = 2.20; bondDistance <= 2.70; bondDistance += 0.05) {
                    // set desired bond distance
                    Molecule newMolecule = newMoleculeTemplate.setDistance(34, 19, bondDistance); // indices
                    /*Atom atom1 = newMolecule.contents.get(18);
                    Atom atom2 = newMolecule.contents.get(33);
                    System.out.println(bondDistance);
                    System.out.println(newMolecule.getDistance(atom1,atom2));
                    System.out.println(newMolecule.directlyConnected(atom1,atom2));*/
                    
                    // write out result
                    String outputFilename = "output/" + filename.replaceAll(".out",".gjf");
                    outputFilename = outputFilename.replaceAll("-ts-", String.format("-%3.0f-", bondDistance*100.0));
                    GaussianInputFile outputGJF = new GaussianInputFile(newMolecule, "title", keywords, tail, false, 0, 1);
                    outputGJF.write(outputFilename);
                    System.out.println("   " + outputFilename);
                }
            }
        }
}
