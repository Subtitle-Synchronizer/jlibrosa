package com.java.jlibrosa.test;

import java.io.IOException;

import com.java.audio.librosafeatures.JLibrosa;
import com.java.audio.librosafeatures.WavFileException;

public class JLibrosaTest_1 {

	public static void main(String[] args) throws IOException, WavFileException {
		// TODO Auto-generated method stub

		JLibrosa jLibrosa = new JLibrosa();
		double [] loadValues = jLibrosa.loadAndRead("audioFiles/100263-2-0-126.wav", 44100, -1);
		double [][] loadValuesAcrossChannels = jLibrosa.loadAndReadAcrossChannels("audioFiles/100263-2-0-126.wav", 44100, -1);
		double [][] stftValues = jLibrosa.generateSTFTFeatures(loadValues, 44100,40);
		
		System.out.println("Test....");
		
		
	}

}
