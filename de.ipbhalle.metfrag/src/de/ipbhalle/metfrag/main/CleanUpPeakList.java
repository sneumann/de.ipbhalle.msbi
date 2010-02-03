package de.ipbhalle.metfrag.main;

import java.util.Vector;

import de.ipbhalle.metfrag.massbankParser.Peak;

public class CleanUpPeakList {
	
	private Vector<Peak> peakList;
	
	/**
	 * Instantiates a new clean up peak list.
	 * 
	 * @param peakList the peak list
	 */
	public CleanUpPeakList(Vector<Peak> peakList)
	{
		this.peakList = peakList;
	}
	
	/**
	 * Gets the cleaned peak list.
	 * 
	 * @return the cleaned peak list
	 */
	public Vector<Peak> getCleanedPeakList(double mass)
	{
		removeHeavyPeaks(mass - 1);
		return this.peakList;
	}
	
	/**
	 * Removes the peaks which correspond to the original molecule. (molecule peak)
	 * 
	 * @param mass the mass
	 */
	private void removeHeavyPeaks(double mass)
	{
		Vector<Peak> tempList = (Vector<Peak>)this.peakList.clone();
		for (int i = 0; i < tempList.size(); i++) {
			if(tempList.get(i).getMass() >= mass)
			{
				//get last element
				int temp = peakList.size();
				peakList.remove(temp - 1);
			}
		}
	}

}
