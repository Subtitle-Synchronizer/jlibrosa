package com.jlibrosa.audio;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import com.jlibrosa.audio.exception.FileFormatNotSupportedException;
import com.jlibrosa.audio.wavFile.WavFileException;

/**
 * 
 * This class tests the JLibrosa functionality for extracting MFCC and STFT Audio features for given Wav file.
 * 
 * @author abhi-rawat1
 *
 */
public class JLibrosaTest {

	public static void main(String[] args) throws IOException, WavFileException, FileFormatNotSupportedException {
		String audioFilePath = "audioFiles/001_children_playing.wav";
		int defaultSampleRate = -1;		//-1 value implies the method to use default sample rate
		int defaultAudioDuration = -1;	//-1 value implies the method to process complete audio duration
		
		JLibrosa jLibrosa = new JLibrosa();

		/* To read the magnitude values of audio files - equivalent to librosa.load('../audioFiles/1995-1826-0003.wav', sr=None) function */
		
		float audioFeatureValues [] = jLibrosa.loadAndRead(audioFilePath, defaultSampleRate, defaultAudioDuration);
		
		ArrayList<Float> audioFeatureValuesList = jLibrosa.loadAndReadAsList(audioFilePath, defaultSampleRate, defaultAudioDuration);
		
		
		for(int i=0;i<10;i++) {
			System.out.printf("%.6f%n", audioFeatureValues[i]);
		}
		
		
		/* To read the no of frames present in audio file*/
		int nNoOfFrames = jLibrosa.getNoOfFrames();
		
		
		/* To read sample rate of audio file */
		int sampleRate = jLibrosa.getSampleRate();
		
		/* To read number of channels in audio file */
		int noOfChannels = jLibrosa.getNoOfChannels();
		
		

		float [][] melSpectrogram = jLibrosa.generateMelSpectroGram(audioFeatureValues, sampleRate, 2048, 128, 256);
		
		System.out.println("/n/n");
		System.out.println("***************************************");
		System.out.println("***************************************");
		System.out.println("***************************************");
		System.out.println("/n/n");

		
		/* To read the MFCC values of an audio file 
		 *equivalent to librosa.feature.mfcc(x, sr, n_mfcc=40) in python
		 * */
		
		float[][] mfccValues = jLibrosa.generateMFCCFeatures(audioFeatureValues, sampleRate, 40);
		
		float[] meanMFCCValues = jLibrosa.generateMeanMFCCFeatures(mfccValues, mfccValues.length, mfccValues[0].length);
		
		System.out.println(".......");
		System.out.println("Size of MFCC Feature Values: (" + mfccValues.length + " , " + mfccValues[0].length + " )");

		for(int i=0;i<1;i++) {
			for(int j=0;j<10;j++) {
				System.out.printf("%.6f%n", mfccValues[i][j]);
			}
		}
		
		
		
		/* To read the STFT values of an audio file 
		 *equivalent to librosa.core.stft(x, sr, n_mfcc=40) in python
		 *Note STFT values return would be complex in nature with real and imaginary values.
		 * */
		
		Complex[][] stftComplexValues = jLibrosa.generateSTFTFeatures(audioFeatureValues, sampleRate, 40);
		
		
		System.out.println(".......");
		System.out.println("Size of STFT Feature Values: (" + stftComplexValues.length + " , " + stftComplexValues[0].length + " )");

		
		
		for(int i =0; i<1;i++) {
			for(int j=0;j<10;j++) {
				double realValue = stftComplexValues[i][j].getReal();
				double imagValue = stftComplexValues[i][j].getImaginary();
				System.out.println("Real and Imag values of STFT are " + realValue + "," + imagValue);
			}
			 
		}
		
			}
	
	
	
	

}
