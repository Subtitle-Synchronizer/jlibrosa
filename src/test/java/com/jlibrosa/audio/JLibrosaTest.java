package com.jlibrosa.audio;

import java.io.IOException;

import com.jlibrosa.audio.wavFile.WavFileException;

/**
 * 
 * This class tests the JLibrosa functionality for extracting MFCC and STFT Audio features for given Wav file.
 * 
 * @author abhi-rawat1
 *
 */
public class JLibrosaTest {

	public static void main(String[] args) throws IOException, WavFileException {
		String audioFilePath = "audioFiles/100263-2-0-126.wav";
		int defaultSampleRate = 44100;
		int defaultAudioDuration = -1;
		
		JLibrosa jLibrosa = new JLibrosa();

		double[] audioFeatureValues = jLibrosa.loadAndRead(audioFilePath, defaultSampleRate, defaultAudioDuration);
		double[][] stftValues = jLibrosa.generateSTFTFeatures(audioFeatureValues, 40);
		double[][] mfccValues = jLibrosa.generateMFCCFeatures(audioFeatureValues, 40);
		
		// double [][] loadValuesAcrossChannels = jLibrosa.loadAndReadAcrossChannels("audioFiles/100263-2-0-126.wav", 44100,-1);

		System.out.println(".......");
		System.out.println("Size of STFT Feature Values: (" + stftValues.length + " , " + stftValues[0].length + " )");

		System.out.println(".......");
		System.out.println("Size of MFCC Feature Values: (" + mfccValues.length + " , " + mfccValues[0].length + " )");
	}

}
