package de.ipbhalle.metfrag.similarity;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmilesParser;

import de.ipbhalle.metfrag.molDatabase.PubChemLocal;
//import de.ipbhalle.metfrag.tools.Render;

public class Similarity {
	
	private Float[][] matrix = null;
	private Map<String, IAtomContainer> candidateToStructure = null;
	private Map<String, Integer> candidateToPosition = null;
	private StringBuilder allSimilarityValues = new StringBuilder();
	private float similarityThreshold;
	
	public Similarity(Map<String, IAtomContainer> candidateToStructure, float similarityThreshold) throws CDKException
	{
		matrix = new Float[candidateToStructure.size()][candidateToStructure.size()];
		this.candidateToStructure = candidateToStructure;
		initializePositions();
		this.similarityThreshold = similarityThreshold;
		calculateSimilarity(similarityThreshold);
	}
	
	//initialize the position of the candidates in the matrix
	private void initializePositions()
	{
		int i = 0;
		candidateToPosition = new HashMap<String, Integer>();
		for (String candidate : candidateToStructure.keySet()) {
			candidateToPosition.put(candidate, i);
			i++;
		}
	}
	
	
	/**
	 * Returns a similarity matrix based on the tanimoto distance
	 * of the fingerprints. (Upper triangular matrix)
	 * 
	 * @param molList the mol list
	 * @param similarityThreshold the similarity threshold
	 * 
	 * @return the float[][]
	 * 
	 * @throws CDKException the CDK exception
	 */
	private Float[][] calculateSimilarity(float similarityThreshold) throws CDKException
	{
		Map<String, BitSet> candidateToFingerprint = new HashMap<String, BitSet>();
		for (String strCandidate : candidateToStructure.keySet()) {
			Fingerprinter f = new Fingerprinter();
			BitSet fp = f.getFingerprint(candidateToStructure.get(strCandidate));
			candidateToFingerprint.put(strCandidate, fp);
		}

		int countJ = 0;
		int countI = 0;
		for (String candidate1 : candidateToStructure.keySet()) {
//			System.out.print(candidate1 + " ");
			for (String candidate2 : candidateToStructure.keySet()) {
				if(countJ < countI || candidate1.equals(candidate2))
				{
//					System.out.print("x ");
					countJ++;
					continue;
				}		
				Float similarity = compareFingerprints(candidateToFingerprint.get(candidate1), candidateToFingerprint.get(candidate2));
				matrix[countI][countJ] = similarity;
				allSimilarityValues.append(similarity + "\n");
				
				countJ++;
//				System.out.print(compareFingerprints(candidateToFingerprint.get(candidate1), candidateToFingerprint.get(candidate2)) + " ");
			}
			countJ = 0;
			countI++;
//			System.out.println("");
		}
		
		
		return matrix;
	}
	
	/**
	 * Gets the all similarity values.
	 * 
	 * @return the all similarity values
	 */
	public StringBuilder getAllSimilarityValues()
	{
		return allSimilarityValues;
	}
	
	/**
	 * Gets the tanimoto distance between two candidate structures.
	 * 
	 * @param candidate1 the candidate1
	 * @param candidate2 the candidate2
	 * 
	 * @return the tanimoto distance
	 */
	public float getTanimotoDistance(String candidate1, String candidate2)
	{
		int pos1 = 0;
		int pos2 = 0;
		try
		{
			pos1 = candidateToPosition.get(candidate1);
			pos2 = candidateToPosition.get(candidate2);
		}
		catch(NullPointerException e)
		{
			return 0;
		}
		if(matrix[pos1][pos2] != null)
			return matrix[pos1][pos2];
		else
			return matrix[pos2][pos1];
	}
	
	
	
	/**
	 * Gets the tanimoto distance from a list of candidates and groups them!.
	 * 
	 * @param candidateGroup the candidate group
	 * @param threshold the threshold
	 * 
	 * @return the tanimoto distance list
	 */
	public List<SimilarityGroup> getTanimotoDistanceList(List<String> candidateGroup)
	{
		
		
		List<SimilarityGroup> groupedCandidates = new ArrayList<SimilarityGroup>();
		for (String cand1 : candidateGroup) {
			SimilarityGroup simGroup = new SimilarityGroup(cand1);
			for (String cand2 : candidateGroup) {
				if(cand1.equals(cand2))
					continue;
				else if(cand1 == null || cand2 == null)
					continue;
				else
				{
					Float tanimoto = getTanimotoDistance(cand1, cand2);
					if(tanimoto > similarityThreshold)
						simGroup.addSimilarCompound(cand2, tanimoto);
				}
			}
			//now add similar compound to the group list
			//if(!isContainedInPreviousResults(groupedCandidates, simGroup))
			groupedCandidates.add(simGroup);
		}
		groupedCandidates = cleanList(groupedCandidates);
		return groupedCandidates;
	}
	
	
	/**
	 * Checks if is identical.
	 * 
	 * @param similar1 the similar1
	 * @param similar2 the similar2
	 * 
	 * @return true, if is identical
	 */
	private boolean isIdentical(List<String> similar1, List<String> similar2)
	{		
		String[] candidatesToCompare = new String[similar1.size()];
		candidatesToCompare = similar1.toArray(candidatesToCompare);
		Arrays.sort(candidatesToCompare);
		
		boolean isAlreadyContained = false;
		
		String[] candidateListTemp = new String[similar2.size()];
		candidateListTemp = similar2.toArray(candidateListTemp);
		Arrays.sort(candidateListTemp);
		
		if(Arrays.equals(candidatesToCompare, candidateListTemp))
		{
			isAlreadyContained = true;
		}

		return isAlreadyContained;
	}
	
	/**
	 * Removes the duplicates.
	 * 
	 * @param groupedCandidates the grouped candidates
	 * 
	 * @return the list< similarity group>
	 */
	private List<SimilarityGroup> removeDuplicates(List<SimilarityGroup> groupedCandidates)
	{		
		List<SimilarityGroup> toRemove = new ArrayList<SimilarityGroup>();
		
		
		if(groupedCandidates.size() == 1)
			return groupedCandidates;
		else
		{
			for (SimilarityGroup simGroup1 : groupedCandidates) {
				
				if(toRemove.contains(simGroup1))
					continue;
				
				String[] candidatesToCompare = new String[simGroup1.getSimilarCandidatesWithBase().size()];
				candidatesToCompare = simGroup1.getSimilarCandidatesWithBase().toArray(candidatesToCompare);
				Arrays.sort(candidatesToCompare);
				
				
				for (SimilarityGroup simGroup2 : groupedCandidates) {
										
					if(simGroup1.equals(simGroup2) || toRemove.contains(simGroup2))
						continue;
					
					List<String> temp = simGroup2.getSimilarCandidatesWithBase();
					String[] candidateListTemp = new String[temp.size()];
					candidateListTemp = simGroup2.getSimilarCandidatesWithBase().toArray(candidateListTemp);
					Arrays.sort(candidateListTemp);
					if(Arrays.equals(candidatesToCompare, candidateListTemp))
						toRemove.add(simGroup2);			
				}
			}
			
			for (SimilarityGroup simGroup : toRemove) {
				groupedCandidates.remove(simGroup);
			}
			return groupedCandidates;
		}
	}
	
	
	/**
	 * Cleans list: remove transitive relations
	 * Only the largest sets are used
	 * 
	 * @param groupedCandidates the grouped candidates
	 * 
	 * @return the list< similarity group>
	 */
	private List<SimilarityGroup> cleanList(List<SimilarityGroup> groupedCandidates)
	{
		List<SimilarityGroup> cleanedList = new ArrayList<SimilarityGroup>();		
		List<SimilarityGroup> alreadyRemoved = new ArrayList<SimilarityGroup>();
		
		Map<String, List<String>> candidatesToCompareMapBase = new HashMap<String, List<String>>();
		Map<String, SimilarityGroup> candidatesToCompareMap = new HashMap<String, SimilarityGroup>();
		for (SimilarityGroup simGroup : groupedCandidates) {
			candidatesToCompareMap.put(simGroup.getCandidateTocompare(), simGroup);
			candidatesToCompareMapBase.put(simGroup.getCandidateTocompare(), simGroup.getSimilarCandidatesWithBase());
		}
		
		if(groupedCandidates.size() == 1)
			cleanedList = groupedCandidates;
		else
		{	
			//initialize with one value
			cleanedList.add(groupedCandidates.get(0));
			for (String key1 : candidatesToCompareMap.keySet()) {
				for (String key2 : candidatesToCompareMap.keySet()) {
					if(key1.equals(key2))
						continue;
					
					List<String> similar1 = candidatesToCompareMapBase.get(key1);
					List<String> similar2 = candidatesToCompareMapBase.get(key2);
					
					if(isIdentical(similar1, similar2))
					{
						continue;
					}
					
					if(similar1.size() > similar2.size() && similar1.containsAll(similar2))
					{
						cleanedList.remove(candidatesToCompareMap.get(key2));
						alreadyRemoved.add(candidatesToCompareMap.get(key2));
					}
					else if(similar2.containsAll(similar1))
					{
						cleanedList.remove(candidatesToCompareMap.get(key1));
						alreadyRemoved.add(candidatesToCompareMap.get(key1));
					}
					else
					{
						if(!cleanedList.contains(candidatesToCompareMap.get(key2)) && !alreadyRemoved.contains(candidatesToCompareMap.get(key2)))
							cleanedList.add(candidatesToCompareMap.get(key2));
						if(!cleanedList.contains(candidatesToCompareMap.get(key1)) && !alreadyRemoved.contains(candidatesToCompareMap.get(key1)))
							cleanedList.add(candidatesToCompareMap.get(key1));
					}
						
				}
			}
			cleanedList = removeDuplicates(cleanedList);
		}
		return cleanedList;
	}
	
	
	/**
	 * Checks if the given structures are isomorph.
	 * 
	 * @param candidate1 the candidate1
	 * @param candidate2 the candidate2
	 * 
	 * @return true, if is isomorph
	 * 
	 * @throws CDKException the CDK exception
	 */
	public boolean isIsomorph(String candidate1, String candidate2) throws CDKException
	{
		IAtomContainer cand1 = this.candidateToStructure.get(candidate1);
		IAtomContainer cand2 = this.candidateToStructure.get(candidate2);
		return UniversalIsomorphismTester.isIsomorph(cand1, cand2);
	}
	
	
	
	
	
	/**
	 * Compare fingerprints by returning the Tanimoto distance.
	 * 
	 * @param bitSet1 the bit set1
	 * @param bitSet2 the bit set2
	 * 
	 * @return the float
	 * 
	 * @throws CDKException the CDK exception
	 */
	private float compareFingerprints(BitSet bitSet1, BitSet bitSet2) throws CDKException
	{
		return Tanimoto.calculate(bitSet1, bitSet2);
	}
	
	
	
	public static void main(String[] args) {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		try {
			IMolecule mol = sp.parseSmiles("C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=CC=C3O");
			IMolecule mol1 = sp.parseSmiles("C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O");
			IMolecule mol2 = sp.parseSmiles("C1C(OC(=O)C2=C(C=C(C=C21)O)O)C3=CC=C(C=C3)O");
			IMolecule mol3 = sp.parseSmiles("C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O");
			IMolecule mol4 = sp.parseSmiles("C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O");
			IMolecule mol5 = sp.parseSmiles("C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=CC=C3)O");
			IMolecule mol6 = sp.parseSmiles("C1C(OC2=CC(=C(C=C2C1=O)O)O)C3=CC=C(C=C3)O");
			IMolecule mol7 = sp.parseSmiles("C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=CC=C3O");
			IMolecule mol8 = sp.parseSmiles("C1C(C(=O)C2=CC(=C(C=C2O1)O)O)C3=CC=C(C=C3)O");
			Map<String, IAtomContainer> testMap = new HashMap<String, IAtomContainer>();
			testMap.put("179999", mol);
			testMap.put("932", mol1);
			testMap.put("10333412", mol2);
			testMap.put("439246", mol3);
			testMap.put("667495", mol4);
			testMap.put("113638", mol5);
			testMap.put("23724670", mol6);
			testMap.put("13889010", mol7);
			testMap.put("125100", mol8);
			
			List<String> cands = new ArrayList<String>();
			cands.add("179999");
			cands.add("932");
			cands.add("10333412");
			cands.add("439246");
			cands.add("667495");
			cands.add("113638");
			cands.add("23724670");
			cands.add("13889010");
			cands.add("125100");
			
			Similarity sim = new Similarity(testMap, (float)0.90);
//			System.out.println(sim.getTanimotoDistance("cand3", "cand1"));
			
			List<SimilarityGroup> groupedCandidates = sim.getTanimotoDistanceList(cands);
			for (SimilarityGroup similarityGroup : groupedCandidates) {
				if(similarityGroup.getSimilarCompounds().size() == 0)
					System.out.print("Single: " + similarityGroup.getCandidateTocompare() + "\n");
				else
				{
					System.out.print("Group of " + similarityGroup.getSimilarCompounds().size() + " " + similarityGroup.getCandidateTocompare() +  ": ");
					for (int i = 0; i < similarityGroup.getSimilarCompounds().size(); i++) {
						System.out.print(similarityGroup.getSimilarCompounds().get(i) + "(" + similarityGroup.getSimilarCompoundsTanimoto().get(i) + ") ");
					}
					System.out.println("");
				}
			}
			
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
