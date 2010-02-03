package de.ipbhalle.metfrag.tools;

public class PPMTool {
	
	/**
	 * Gets the PPM for a specified peak (m/z).
	 * 
	 * @param peak the peak
	 * @param ppm the ppm
	 * 
	 * @return the allowed error with the given ppm
	 */
	public static double getPPM(double peak, double ppm)
	{
		double ret = 0.0;
		//calculate the allowed error for the given peak m/z in Th
		ret = (peak/1000000.0) * ppm;
		return ret;
	}

}
