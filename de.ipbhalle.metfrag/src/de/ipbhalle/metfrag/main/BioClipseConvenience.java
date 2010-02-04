package de.ipbhalle.metfrag.main;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import de.ipbhalle.metfrag.chemspiderClient.ChemSpider;
import de.ipbhalle.metfrag.fragmenter.Fragmenter;
import de.ipbhalle.metfrag.keggWebservice.KeggWebservice;
import de.ipbhalle.metfrag.main.AssignFragmentPeak;
import de.ipbhalle.metfrag.main.CleanUpPeakList;
import de.ipbhalle.metfrag.main.PeakMolPair;
import de.ipbhalle.metfrag.main.WrapperSpectrum;
import de.ipbhalle.metfrag.massbankParser.Peak;
import de.ipbhalle.metfrag.molDatabase.PubChemLocal;
import de.ipbhalle.metfrag.pubchem.PubChemWebService;
import de.ipbhalle.metfrag.scoring.Scoring;
import de.ipbhalle.metfrag.tools.MolecularFormulaTools;
import de.ipbhalle.metfrag.tools.PPMTool;


public class BioClipseConvenience {
	
	
	private boolean sumFormulaRedundancyCheck = false;

	private String sumFormula = "";
	private int mode = 1;
	private double exactMass = 272.06847;


	private double mzabs = 0.01;
	private double mzppm = 50;
	private int count = 0;
	private String peaks = "119.051 467.616 45\n" +
	   "123.044 370.662 36\n" +
	   "147.044 6078.145 606\n" +
	   "153.019 10000.0 999\n" +
	   "179.036 141.192 13\n" +
	   "189.058 176.358 16\n" +
	   "273.076 10000.000 999\n" +
	   "274.083 318.003 30\n";	
	private List<IAtomContainer> molecules;	
	
	/**
	 * Instantiates a new database retrieval just for testing with the default values.
	 */
	public BioClipseConvenience()
	{
		
	}
	
	/**
	 * Instantiate and set all parameters needed for MetFrag.
	 * 
	 * @param mzabs the mzabs
	 * @param mzppm the mzppm
	 * @param peaks the peaks
	 * @param database the database
	 * @param sumFormulaRedundancyCheck the sum formula redundancy check
	 * @param sumFormula the sum formula
	 * @param mode the mode
	 * @param exactMass the exact mass
	 */
	public BioClipseConvenience(double mzabs, double mzppm, String peaks, List<IAtomContainer> molecules, 
			boolean sumFormulaRedundancyCheck, String sumFormula, int mode, double exactMass)
	{
		this.mzabs = mzabs;
		this.mzppm = mzppm;
		this.peaks = peaks;
		this.molecules = molecules;
		this.sumFormulaRedundancyCheck = sumFormulaRedundancyCheck;
		this.sumFormula = sumFormula;
		this.mode = mode;
		this.exactMass = exactMass;
	}


	/**
	 * MetFrag convenience method!!!! :)
	 * 
	 * Retrieve candidate molecules from the selected database and fragment them.
	 * An array of scores is returned 
	 * 
	 * @throws Exception the exception
	 */
	public List<Double> metFrag() throws Exception
	{

		WrapperSpectrum spectrum = new WrapperSpectrum(this.peaks, mode, exactMass);
		HashMap<Integer, ArrayList<String>> scoreMap = new HashMap<Integer, ArrayList<String>>();
		HashMap<Double, Vector<String>> realScoreMap = new HashMap<Double, Vector<String>>();

		count = 0;

		for (int c = 0; c < molecules.size(); c++) {

			Vector<Peak> listOfPeaks = new Vector<Peak>();
			IAtomContainer molecule = molecules.get(c);

			//skip if molecule is not connected
			boolean isConnected = ConnectivityChecker.isConnected(molecule);
			if(!isConnected)
				continue;

			try
			{
				//add hydrogens
				CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(molecule.getBuilder());

				for (IAtom atom : molecule.atoms()) {
					IAtomType type = matcher.findMatchingAtomType(molecule, atom);
					AtomTypeManipulator.configure(atom, type);
				}
				CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
				hAdder.addImplicitHydrogens(molecule);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
			}
			//there is a bug in cdk?? error happens when there is a S or Ti in the molecule
			catch(IllegalArgumentException e)
			{
				//skip it
				continue;
			}


			//render original compound....thats the first picture in the list
			int countTemp = 0;


			Double massDoubleOrig = MolecularFormulaTools.getMonoisotopicMass(MolecularFormulaManipulator.getMolecularFormula(molecule));
			massDoubleOrig = (double)Math.round((massDoubleOrig)*10000)/10000;
			String massOrig = massDoubleOrig.toString();
			countTemp++;

			Fragmenter fragmenter = new Fragmenter((Vector<Peak>)spectrum.getPeakList().clone(), mzabs, mzppm, mode, true, this.sumFormulaRedundancyCheck, false, false);     
			List<IAtomContainer> l = null;
			try
			{
				l = fragmenter.generateFragmentsInMemory(molecule, true, 2);
			}
			catch(OutOfMemoryError e)
			{
				System.out.println("OUT OF MEMORY ERROR! " + molecule.getID());
				continue;
			}


			List<IAtomContainer> fragments = l;  

			//get the original peak list again
			Vector<Peak> peakListParsed = spectrum.getPeakList();


			//clean up peak list
			CleanUpPeakList cList = new CleanUpPeakList((Vector<Peak>) peakListParsed.clone());
			Vector<Peak> cleanedPeakList = cList.getCleanedPeakList(spectrum.getExactMass());


			//now find corresponding fragments to the mass
			AssignFragmentPeak afp = new AssignFragmentPeak();
			afp.setHydrogenTest(true);
			afp.AssignFragmentPeak(fragments, cleanedPeakList, mzabs, mzppm, spectrum.getMode(), false);

			Vector<PeakMolPair> hits = afp.getHits();			


			//now "real" scoring --> depends on intensities
			Scoring score = new Scoring(spectrum.getPeakList());
			double currentScore = score.computeScoring(afp.getHitsMZ());

			//save score in hashmap...if there are several hits with the same score --> vector of strings
			if(realScoreMap.containsKey(currentScore))
			{
				Vector<String> tempList = realScoreMap.get(currentScore);
				tempList.add(molecule.getID());
				realScoreMap.put(currentScore, tempList);
			}
			else
			{
				Vector<String> temp = new Vector<String>();
				temp.add(molecule.getID());
				realScoreMap.put(currentScore, temp);
			}


			//save score in hashmap...if there are several hits with the same
			//amount of identified peaks --> ArrayList
			if(scoreMap.containsKey(hits.size()))
			{
				ArrayList<String> tempList = scoreMap.get(hits.size());
				tempList.add(molecule.getID());
				scoreMap.put(hits.size(), tempList);
			}
			else
			{
				ArrayList<String> temp = new ArrayList<String>();
				temp.add(molecule.getID());
				scoreMap.put(hits.size(), temp);
			}

			Vector<Double> peaks = new Vector<Double>();
			Vector<Double> intensities = new Vector<Double>();

			//get all the identified peaks
			for (int i = 0; i < hits.size(); i++) {
				listOfPeaks.add(hits.get(i).getPeak());
				peaks.add(hits.get(i).getPeak().getMass());
				intensities.add(hits.get(i).getPeak().getRelIntensity());
				//all found peaks are later on marked in the spectrum
				//xyFound.add(hits.get(i).getPeak().getMass(), hits.get(i).getPeak().getRelIntensity());
			}


			List<IAtomContainer> hitsList = new ArrayList<IAtomContainer>();
			for (int i = 0; i < hits.size(); i++) {
				hitsList.add(AtomContainerManipulator.removeHydrogens(hits.get(i).getFragment()));
				//Render.Highlight(AtomContainerManipulator.removeHydrogens(molecule), hitsList , Double.toString(hits.get(i).getPeak()));
			}

			System.out.println(" " + molecule.getID() + " " + afp.getHits().size() + " " + currentScore);
			//done fragments
			count++;
		}

		HashMap<Double, Vector<String>> scoresNormalized = normalize(realScoreMap);
		Double[] scores = new Double[scoresNormalized.size()];
		scores = scoresNormalized.keySet().toArray(scores);
		//		Arrays.sort(scores);
		List<Double> result;
		result = Arrays.asList(scores);		
		return result;
	}

	
	/**
	 * Normalize the score between 0 and 1.
	 */
	private HashMap<Double, Vector<String>> normalize(HashMap<Double, Vector<String>> realScoreMap)
	{
		HashMap<Double, Vector<String>> ret = new HashMap<Double, Vector<String>>();
		
		double maxScore = 0; 
		for (Double score : realScoreMap.keySet()) {
			if(score > maxScore)
			{
				maxScore = score;
			}
		}
		
		for (Double score : realScoreMap.keySet()) {
			Vector<String> list = realScoreMap.get(score);
			ret.put((Math.round((score / maxScore) * 1000) / 1000.0), list);
		}
		
		return ret;
	}
	
	
	public static void main(String[] args) {
//		MetFlowConvenience test = new MetFlowConvenience();
//		try {
//			System.out.println(test.metFrag());
//		} catch (Exception e) {
//			e.printStackTrace();
//		}

		
		try {
			String peaks = "119.051 467.616 45\n123.044 370.662 36\n" +
			"147.044 6078.145 606\n153.019 10000.0 999\n179.036 141.192 13\n189.058 176.358 16\n" +
			"273.076 10000.000 999\n274.083 318.003 30\n";
			String smiles = "C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O";

			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			//parse smiles
			IAtomContainer molecule = sp.parseSmiles(smiles);
			//configure atoms
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
			//add all hydrogens explicitly
			CDKHydrogenAdder adder1 = CDKHydrogenAdder.getInstance(molecule.getBuilder());
	        adder1.addImplicitHydrogens(molecule);
	        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule); 

	        List<IAtomContainer> molecules = new ArrayList<IAtomContainer>();
	        molecules.add(molecule);
	        
			BioClipseConvenience bcc = new BioClipseConvenience(0.01, 50.0, peaks, molecules, true, "", 1, 272.06847);
			//System.out.println(bcc.metFrag());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
