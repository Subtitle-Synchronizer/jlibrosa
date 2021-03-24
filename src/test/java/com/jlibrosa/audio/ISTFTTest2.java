package com.jlibrosa.audio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.math3.complex.Complex;

public class ISTFTTest2 {

	public static void main(String[] args) throws NumberFormatException, IOException {
		// TODO Auto-generated method stub

		Complex [][] compArray = readFromFile();
		JLibrosa jLibrosa = new JLibrosa();
		
		float[] magValues = jLibrosa.generateInvSTFTFeatures(compArray,
                44100, 40,256, 128, 64);
		System.out.println("test");
	}
	
	public static Complex[][] readFromFile() throws NumberFormatException, IOException {
		String savedGameFile = "/Users/vishrud/Desktop/Vasanth/Technology/Mobile-ML/Spleeter_TF2.0/local/output2darray.csv";
		Complex[][] board = new Complex[129][6881];
		BufferedReader reader = new BufferedReader(new FileReader(savedGameFile));
		String line = "";
		int row = 0;
		while((line = reader.readLine()) != null)
		{
		   String[] cols = line.split(","); //note that if you have used space as separator you have to split on " "
		   for(int i=0;i<cols.length;i=i+1) {
			   String procStr = cols[i].replaceAll("[(]", "");
			   procStr = procStr.replaceAll("[)]", "");
			   String [] splStr = null;
			   
			   splStr = procStr.split("xx");
			   String realStr = splStr[0];
			   String imgStr = splStr[1];
			   
			   double real = Double.parseDouble(realStr);
			   double imag = Double.parseDouble(imgStr);
			   board[row][i] = new Complex(real, imag);
		   }
		   
		   
		   row++;
		}
		reader.close();
		return board;
	}

}
