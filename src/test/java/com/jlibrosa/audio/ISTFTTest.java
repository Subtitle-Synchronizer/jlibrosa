package com.jlibrosa.audio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.math3.complex.Complex;

public class ISTFTTest {

	public static void main(String[] args) throws NumberFormatException, IOException {
		// TODO Auto-generated method stub

		Complex [][] compArray = readFromFile();
		JLibrosa jLibrosa = new JLibrosa();
		
		float[] magValues = jLibrosa.generateInvSTFTFeatures(compArray,
                44100, 40,4096, 128, 1024);
	}
	
	public static Complex[][] readFromFile() throws NumberFormatException, IOException {
		String savedGameFile = "/Users/vishrud/Downloads/twodarray.txt";
		Complex[][] board = new Complex[2049][212];
		BufferedReader reader = new BufferedReader(new FileReader(savedGameFile));
		String line = "";
		int row = 0;
		while((line = reader.readLine()) != null)
		{
		   String[] cols = line.split(","); //note that if you have used space as separator you have to split on " "
		   int col = 0;
		   int counter =0;
		   for(int i=0;i<cols.length;i=i+2) {
			   String realStr = cols[i].replaceAll("[(]", "");
			   String imagStr = cols[i+1].replaceAll("[)]", "");
			   
			   double real = Double.parseDouble(realStr);
			   double imag = Double.parseDouble(imagStr);
			   board[row][counter] = new Complex(real, imag);
		   }
		   
		   
		   row++;
		}
		reader.close();
		return board;
	}

}
