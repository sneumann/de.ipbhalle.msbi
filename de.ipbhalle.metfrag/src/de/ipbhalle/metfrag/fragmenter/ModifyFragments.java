package de.ipbhalle.metfrag.fragmenter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;


import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.aromaticity.AromaticityCalculator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;

import de.ipbhalle.metfrag.massbankParser.Peak;
//import de.ipbhalle.metfrag.tools.Render;


/**
 * The Class ModifyFragments.
 */
public class ModifyFragments {
	
	/** The all rings. */
	private List<IBond> aromaticBonds;
	private IRingSet allRings;
	
	/**
	 * Post process a fragment. --> find rearrangements with oxygen. (water)
	 * 
	 * @param original the original
	 * @param hydrogenUsed the hydrogen used
	 * @param aromaticBonds the aromatic bonds
	 * @param allRings the all rings
	 * 
	 * @return the i atom container set
	 * 
	 * @throws CDKException the CDK exception
	 * @throws CloneNotSupportedException the clone not supported exception
	 */
	public List<IAtomContainer> postProcess(IAtomContainer original, List<IBond> aromaticBonds, IRingSet allRingsOrig, boolean hydrogenUsed) throws CDKException, CloneNotSupportedException
	{
		//Render.Draw(original, "Original Main");
		
		List<IAtomContainer> ret = new ArrayList<IAtomContainer>();
		this.aromaticBonds = aromaticBonds;
		
		if(allRingsOrig.getAtomContainerCount() > 0)
		{
			//get the rings which are not broken up yet
			HashMap<IBond, Integer> bondMap = new HashMap<IBond, Integer>();
			int count = 0;
			
			for (IBond bondOrig : original.bonds()) {
				bondMap.put(bondOrig, count);
				count++;
			}
			
			
			//check for rings which are not broken up!
			IRingSet validRings = new RingSet();
			for (int i = 0; i < allRingsOrig.getAtomContainerCount(); i++) {
				int bondcount = 0;
				
				for (IBond bondRing : allRingsOrig.getAtomContainer(i).bonds()) {
					if(bondMap.containsKey(bondRing))
						bondcount++;
				}
				if(bondcount == allRingsOrig.getAtomContainer(i).getBondCount())
					validRings.addAtomContainer(allRingsOrig.getAtomContainer(i));
			}
			//rings which are not split up
			this.allRings = validRings;
		}
		else
		{
			//empty RingSet
			this.allRings = new RingSet();
		}
		
		//AllRingsFinder allRingsFinder = new AllRingsFinder();
        //allRingsFinder.setTimeout(100000);
        //allRings = allRingsFinder.findAllRings(original);
		
		IAtomContainer temp = new AtomContainerMetFrag();
		List<IAtom> doneAtoms = new ArrayList<IAtom>();
		
        for (IBond bond : original.bonds()) {
          
            // lets see if it is a terminal bond...check for oxygen
			for (IAtom atom : bond.atoms()) {
				if(doneAtoms.contains(atom))
				{
					continue;
				}
				else
				{
					doneAtoms.add(atom);
				}
				//a possible hit...oxygen is not contained in a ring
				if (atom.getSymbol().startsWith("O") && !allRings.contains(atom))
				{
					temp = checkForWater(original, atom, hydrogenUsed);
					if(temp.getAtomCount() > 0)
					{
						ret.add(temp);
						//create a atom container
						temp = new AtomContainerMetFrag();
					}
				}
				//a possible hit...N is not contained in a ring
				if (atom.getSymbol().startsWith("N") && !allRings.contains(atom))
				{
					temp = splitN(original, atom);
					if(temp.getAtomCount() > 0)
					{
						ret.add(temp);
						//create a atom container
						temp = new AtomContainerMetFrag();
					}
				}
			}

        }
				
    	return ret;
	}
	
	
	/**
	 * Split up NH2 + H. If there is a another hydrogen on the next C atom.
	 * 
	 * @param frag the frag
	 * @param candidateNitrogen the candidate nitrogen
	 * 
	 * @return the i atom container
	 */
	private IAtomContainer splitN(IAtomContainer frag, IAtom candidateNitrogen)
	{
		IAtomContainer ret = new AtomContainerMetFrag();
		//create a copy from the original fragment
    	List<IBond> part = new ArrayList<IBond>();
    	part = traverse(frag, candidateNitrogen, part);
    	IAtomContainer fragCopy = makeAtomContainer(candidateNitrogen, part);
    	
    	IBond bondToRemove = null;
    	IAtom hydrogenToRemove = null;
    	IAtom firstH = null;
    	IAtom secondH = null;
    	int count = 0;
    	
    	//NH2
		//there should only be one connected atom (execpt the hydrogens)
    	List<IBond> bondList = fragCopy.getConnectedBondsList(candidateNitrogen);
    	for (int i = 0; i < bondList.size(); i++) {
    		for (IAtom atom : bondList.get(i).atoms()) {
        		//thats the nitrogen
        		if(atom.getSymbol().startsWith("N"))
        			continue;
        		//found hydrogen
        		if(atom.getSymbol().startsWith("H"))
        		{
        			if(firstH == null)
        				firstH = atom;
        			else
        				secondH = atom;
        			count++;
        		}
        		//found C atom
        		if(atom.getSymbol().startsWith("C"))
        		{
        			//this carbon should have a hydrogen connected to it
    				bondToRemove = new Bond(candidateNitrogen, atom);
    				List<IBond> bondList1 = fragCopy.getConnectedBondsList(atom);
    				for (int j = 0; j < bondList1.size(); j++) {
    	        		for (IAtom atom1 : bondList1.get(j).atoms()) {
                    		//thats the hydrogen
                    		if(atom1.getSymbol().startsWith("H"))
                    		{
                    			hydrogenToRemove = atom1;
                    		}
        				}
    				}
        		}
    		}
    	}
    	//now split up NH3
    	if(hydrogenToRemove != null && count == 2)
    	{
    		//remove all the atoms
        	fragCopy.removeAtomAndConnectedElectronContainers(hydrogenToRemove);
        	fragCopy.removeAtomAndConnectedElectronContainers(firstH);
        	fragCopy.removeAtomAndConnectedElectronContainers(secondH);
        	fragCopy.removeAtomAndConnectedElectronContainers(candidateNitrogen);
        	//and now the bonds
			fragCopy.removeBond(bondToRemove);
			fragCopy.removeBond(candidateNitrogen, firstH);
			fragCopy.removeBond(candidateNitrogen, secondH);
			return fragCopy;
    	}
    	else
    		return ret;
    	
	}
    
    
    /**
     * Check for 2 hydrogens in the next 2 Atoms.
     * 
     * @param candidateOxygen the candidate oxygen atom
     * @param frag the frag
     * @param proton the proton
     * 
     * @return true, if successful
     */
    private IAtomContainer checkForWater(IAtomContainer frag, IAtom candidateOxygen, boolean proton)
    {

    	IAtomContainer ret = new AtomContainerMetFrag();
    	
    	//create a copy from the original fragment
    	//Render.Draw(frag, "Original");
    	List<IBond> part = new ArrayList<IBond>();
    	part = traverse(frag, candidateOxygen, part);
    	IAtomContainer fragCopy = makeAtomContainer(candidateOxygen, part);
    	//Render.Draw(fragCopy, "Copy");
    	
    	//there should only be one connected atom
    	List<IBond> bondList = fragCopy.getConnectedBondsList(candidateOxygen);
    	
    	//first C atom
    	IAtom tempAtom = new Atom();
    	
    	//atom and bonds to be removed
    	IAtom firstAtom = new Atom();
    	IBond firstBond = new Bond();
    	
    	//atom and bonds to be removed
    	IAtom secondAtom = new Atom();
    	IBond secondBond = new Bond();
    	
    	boolean firstHydrogen = false;
    	boolean secondHydrogen = false;
    	
    	boolean tempAtomCheck = false;
    	
    	//at most 2 bonds between the oxygen and other atoms (at most 1 H and 2 C)
    	for (int i = 0; i < bondList.size(); i++) {
    		for (IAtom atom : bondList.get(i).atoms()) {
        		//thats the oxygen
        		if(atom.getSymbol().startsWith("O"))
        			continue;
        		//save atom...make sure O is not between 2 C
        		if(atom.getSymbol().startsWith("C") && !tempAtomCheck)
        		{
        			tempAtom = atom;
        			tempAtomCheck = true;
        		}
        		else if(atom.getSymbol().startsWith("H"))
        		{
        			//one hydrogen found
        			firstHydrogen = true;
        			//this will be removed with the oxygen
        			firstAtom = atom;
        			firstBond = bondList.get(i);
        		}
        		else
        		{
        			//not a possible hit
        			return ret;
        		}
        	}
		}
    	
    	
    	
    	//bond between oxygen and carbon
    	//IBond tempBond = fragCopy.getBond(candidateOxygen, tempAtom);
    	
    	Vector<IAtom> carbonList = new Vector<IAtom>();
    	
    	
    	
    	//all bonds from carbon...check for hydrogen
    	bondList = fragCopy.getConnectedBondsList(tempAtom);
    	for (int i = 0; i < bondList.size(); i++) {
			
    		if(bondList.get(i).getConnectedAtom(tempAtom).getSymbol().startsWith("H") && !firstHydrogen)
    		{
    			//one hydrogen found
    			firstHydrogen = true;
    			//the atom to be removed
    			firstAtom = bondList.get(i).getConnectedAtom(tempAtom);
    			//the bond to be removed
    			firstBond = fragCopy.getBond(tempAtom, bondList.get(i).getConnectedAtom(tempAtom));	
    		}
			else if(firstHydrogen && bondList.get(i).getConnectedAtom(tempAtom).getSymbol().startsWith("H"))
			{
				secondHydrogen = true;
				//the atom to be removed
    			secondAtom = bondList.get(i).getConnectedAtom(tempAtom);
    			//the bond to be removed
    			secondBond = fragCopy.getBond(tempAtom, bondList.get(i).getConnectedAtom(tempAtom));
    			break;
			}
    			
    		
    		if(bondList.get(i).getConnectedAtom(tempAtom).getSymbol().startsWith("C"))
    		{
    			//add carbon to list to find more hydrogen
    			carbonList.add(bondList.get(i).getConnectedAtom(tempAtom));
    		}
		}
    	
    	
    	//now check those carbons in the list for another hydrogen atom
    	IBond bondToBeDouble = new Bond();
    	IAtom tempAtom1 = new Atom();
    	if(!secondHydrogen)
    	{
	    	for (int i = 0; i < carbonList.size(); i++) 
	    	{
	    		//get next carbon
	    		bondList = fragCopy.getConnectedBondsList(carbonList.get(i));
	    		//now check carbon from list for another hydrogen
	    		for (int j = 0; j < bondList.size(); j++) {
	    			if(bondList.get(j).getConnectedAtom(carbonList.get(i)).getSymbol().startsWith("H"))
	        		{
	        			//the atom to be removed
	        			secondAtom = bondList.get(j).getConnectedAtom(carbonList.get(i));
	        			//the bond to be removed
	        			secondBond = fragCopy.getBond(carbonList.get(i), bondList.get(j).getConnectedAtom(carbonList.get(i)));
	        			
	        			secondHydrogen = true;
	        			
	        			bondToBeDouble = fragCopy.getBond(tempAtom, carbonList.get(i));
	        			tempAtom1 = carbonList.get(i);
	        			break;
	        		}
	    		}
	    		
			}
    	}  	
    	
    	
    	//rearrangement
    	if(firstHydrogen && secondHydrogen)
    	{
	    	//Render.Draw(frag, "Hit");
			//remove the hydrogen and the bond to the hydrogen 
			fragCopy.removeAtomAndConnectedElectronContainers(secondAtom);
			fragCopy.removeBond(secondBond);
			//Render.Draw(fragCopy, "Rearranged!!!1");
			//from the first hydrogen too
			fragCopy.removeAtomAndConnectedElectronContainers(firstAtom);
			fragCopy.removeBond(firstBond);
			//Render.Draw(fragCopy, "Rearranged!!!2");
			
			//now remove the oxygen too
			fragCopy.removeAtomAndConnectedElectronContainers(candidateOxygen);
			fragCopy.removeBond(tempAtom,candidateOxygen);
			
			
			//create a double bond from the bond between the carbon
			//IBond doubleBond = frag.getBond(tempAtom, carbonList.get(i));
			
			//only change to double bond if the ring is not aromatic and there is a valid bond
			if(!this.aromaticBonds.contains(bondToBeDouble) && bondToBeDouble.getID() != null)
			{
				String temp = bondToBeDouble.getID();
				fragCopy.removeBond(bondToBeDouble);
				IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
				IBond b = builder.newBond(tempAtom, tempAtom1, IBond.Order.DOUBLE);
				b.setID(temp);
				//b.setID(Fragmenter.getBondNumber() + "");
				//increase bond number
				Fragmenter.setBondNumber(Fragmenter.getBondNumber() + 1);
		        
				fragCopy.addBond(b);
			}
			
	        //Render.Draw(fragCopy, "Rearranged");
	        return fragCopy;
    	}
    	
    	return ret;
    }
    
    
    /**
     * Traverse like in a the real fragmenter. The main purpose is to create a copy of
     * of the old fragment
     * 
     * @param atomContainer the atom container
     * @param atom the atom
     * @param bondList the bond list
     * 
     * @return the list< i bond>
     */
    private List<IBond> traverse(IAtomContainer atomContainer, IAtom atom, List<IBond> bondList) {
        List<IBond> connectedBonds = atomContainer.getConnectedBondsList(atom);
        for (IBond aBond : connectedBonds) {
            if (bondList.contains(aBond))
                continue;
            bondList.add(aBond);
            IAtom nextAtom = aBond.getConnectedAtom(atom);
            if (atomContainer.getConnectedAtomsCount(nextAtom) == 1)
                continue;
            traverse(atomContainer, nextAtom, bondList);
        }
        return bondList;
    }
    
    private IAtomContainer makeAtomContainer(IAtom atom, List<IBond> parts) {
        IAtomContainer partContainer = new AtomContainerMetFrag();
        partContainer.addAtom(atom);
        for (IBond aBond : parts) {
            for (IAtom bondedAtom : aBond.atoms()) {
                if (!partContainer.contains(bondedAtom))
                    partContainer.addAtom(bondedAtom);
            }
            partContainer.addBond(aBond);
        }
        return partContainer;
    }

}
