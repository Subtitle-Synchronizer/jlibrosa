package com.java.jlibrosa.test;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class FFTTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub

		 double[] xValues = new double[]{0.5,1.5,1.2,0.5,0.8,1.5,1.9,2.9};

	    double[] tempConversion = new double[xValues.length];
	    double[] tempImag = new double[xValues.length];

	    FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.UNITARY);
	    try {           
	        Complex[] complx = transformer.transform(xValues, TransformType.FORWARD);

	        for (int i = 0; i < complx.length; i++) {               
	            double rr = (complx[i].getReal());
	            
	            double ri = (complx[i].getImaginary());

	            tempConversion[i] = rr; // Math.sqrt((rr * rr) + (ri * ri));
	            tempImag[i] = ri;
	        }

	    } catch (IllegalArgumentException e) {
	        System.out.println(e);
	    }
		
	    System.out.println(tempConversion);
	}

}
