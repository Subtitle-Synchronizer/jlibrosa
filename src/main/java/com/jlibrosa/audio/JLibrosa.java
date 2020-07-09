package com.jlibrosa.audio;

import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.apache.commons.math3.complex.Complex;

import com.jlibrosa.audio.process.AudioFeatureExtraction;
import com.jlibrosa.audio.wavFile.WavFile;
import com.jlibrosa.audio.wavFile.WavFileException;

/**
 * 
 * This Class is an equivalent of Python Librosa utility used to extract the Audio features from given Wav file.
 * 
 * @author abhi-rawat1
 *
 */
public class JLibrosa {
	private int BUFFER_SIZE = 4096;
	/*private int mNumFrames;
	private int mSampleRate;
	private int mChannels;
*/
	
	public JLibrosa() {

	}

	/**
	 * This function is used to load the audio file and read its Numeric Magnitude
	 * Feature Values.
	 * 
	 * @param path
	 * @param sr
	 * @param readDurationInSec
	 * @return
	 * @throws IOException
	 * @throws WavFileException
	 */
	public double[][] loadAndReadAcrossChannels(String path, int sr, int readDurationInSec)
			throws IOException, WavFileException {
		double[][] magValues = readMagnitudeValuesFromFile(path, sr, readDurationInSec);
		return magValues;
	}

	/**
	 * This function is used to load the audio file and read its Numeric Magnitude
	 * Feature Values.
	 * 
	 * @param path
	 * @param sr
	 * @param readDurationInSeconds
	 * @return
	 * @throws IOException
	 * @throws WavFileException
	 */
	private double[][] readMagnitudeValuesFromFile(String path, int sampleRate, int readDurationInSeconds)
			throws IOException, WavFileException {

		File sourceFile = new File(path);
		WavFile wavFile = null;

		wavFile = WavFile.openWavFile(sourceFile);
		int mNumFrames = (int) (wavFile.getNumFrames());
		int mSampleRate = (int) wavFile.getSampleRate();
		int mChannels = wavFile.getNumChannels();

		if (readDurationInSeconds != -1) {
			mNumFrames = readDurationInSeconds * mSampleRate;
		}

		if (sampleRate != -1) {
			mSampleRate = sampleRate;
		}

		// Read the magnitude values across both the channels and save them as part of
		// multi-dimensional array
		double[][] buffer = new double[mChannels][mNumFrames];
		int frameOffset = 0;
		int loopCounter = ((mNumFrames * mChannels) / BUFFER_SIZE) + 1;
		for (int i = 0; i < loopCounter; i++) {
			frameOffset = wavFile.readFrames(buffer, mNumFrames, frameOffset);
		}

		return buffer;

	}

	
	
	/**
	 * This method is used to read number of frames present in the audio file.
	 * 
	 * @param filePath
	 * @return no of frames present in audio file
	 * @throws IOException
	 * @throws WavFileException
	 */
	
	public int getNoOfFrames(String filePath) throws IOException, WavFileException {
		File sourceFile = new File(filePath);
		WavFile wavFile = null;
		int nNoOfFrames = 0;
		
		wavFile = WavFile.openWavFile(sourceFile);
		nNoOfFrames = (int) (wavFile.getNumFrames());
		return nNoOfFrames;

	}
	

	/**
	 * This method is used to read number of frames present in the audio file.
	 * 
	 * @param filePath
	 * @return no of frames present in audio file
	 * @throws IOException
	 * @throws WavFileException
	 */
	
	public int getSampleRate(String filePath) throws IOException, WavFileException {
		File sourceFile = new File(filePath);
		WavFile wavFile = null;
		int nSampleRate = 0;
		
		wavFile = WavFile.openWavFile(sourceFile);
		nSampleRate = (int) (wavFile.getSampleRate());
		return nSampleRate;

	}

	
	/**
	 * This method is used to read number of frames present in the audio file.
	 * 
	 * @param filePath
	 * @return no of frames present in audio file
	 * @throws IOException
	 * @throws WavFileException
	 */
	
	public int getNoOfChannels(String filePath) throws IOException, WavFileException {
		File sourceFile = new File(filePath);
		WavFile wavFile = null;
		int nChannels = 0;
		
		wavFile = WavFile.openWavFile(sourceFile);
		nChannels = (int) (wavFile.getNumChannels());
		return nChannels;

	}
	
	
	
	/**
	 * This function calculates and returns the MFCC values of given Audio Sample
	 * values.
	 * 
	 * @param magValues
	 * @param nMFCC
	 * @return
	 */
	public double[][] generateMFCCFeatures(double[] magValues, int mSampleRate, int nMFCC) {

		AudioFeatureExtraction mfccConvert = new AudioFeatureExtraction();
		mfccConvert.setSampleRate(mSampleRate);
		mfccConvert.setN_mfcc(nMFCC);
		float [] mfccInput = mfccConvert.extractMFCCFeatures(magValues); //extractMFCCFeatures(magValues);
		
		int nFFT = mfccInput.length / nMFCC;
		double[][] mfccValues = new double[nMFCC][nFFT];

		// loop to convert the mfcc values into multi-dimensional array
		for (int i = 0; i < nFFT; i++) {
			int indexCounter = i * nMFCC;
			int rowIndexValue = i % nFFT;
			for (int j = 0; j < nMFCC; j++) {
				mfccValues[j][rowIndexValue] = mfccInput[indexCounter];
				indexCounter++;
			}
		}

		return mfccValues;
		
	}

	/**
	 * This function calculates and return the Mean MFCC values.
	 * 
	 * @param mfccValues
	 * @param nMFCC
	 * @param nFFT
	 * @return
	 */
	public float[] generateMeanMFCCFeatures(double[][] mfccValues, int nMFCC, int nFFT) {
		// code to take the mean of mfcc values across the rows such that
		// [nMFCC x nFFT] matrix would be converted into
		// [nMFCC x 1] dimension - which would act as an input to tflite model

		float[] meanMFCCValues = new float[nMFCC];
		for (int p = 0; p < nMFCC; p++) {
			double fftValAcrossRow = 0;
			for (int q = 0; q < nFFT; q++) {
				fftValAcrossRow = fftValAcrossRow + mfccValues[p][q];
			}
			double fftMeanValAcrossRow = fftValAcrossRow / nFFT;
			meanMFCCValues[p] = (float) fftMeanValAcrossRow;
		}
		return meanMFCCValues;
	}

	/**
	 * This function calculates and returns the STFT values of given Audio Sample
	 * values. STFT stands for Short Term Fourier Transform
	 * 
	 * @param magValues
	 * @param nMFCC
	 * @return
	 */
	public Complex [][] generateSTFTFeatures(double[] magValues, int mSampleRate, int nMFCC) {
		AudioFeatureExtraction featureExtractor = new AudioFeatureExtraction();
		featureExtractor.setSampleRate(mSampleRate);
		featureExtractor.setN_mfcc(nMFCC);
		Complex [][] stftValues = featureExtractor.extractSTFTFeaturesAsComplexValues(magValues);
		return stftValues;
	}

	/**
	 * This function loads the audio file, reads its Numeric Magnitude Feature
	 * values and then takes the mean of amplitude values across all the channels and
	 * convert the signal to mono mode
	 * 
	 * @param path
	 * @param sampleRate
	 * @param readDurationInSeconds
	 * @return
	 * @throws IOException
	 * @throws WavFileException
	 */
	public double[] loadAndRead(String path, int sampleRate, int readDurationInSeconds)
			throws IOException, WavFileException {

		double[][] magValueArray = readMagnitudeValuesFromFile(path, sampleRate, readDurationInSeconds);

		DecimalFormat df = new DecimalFormat("#.#####");
		df.setRoundingMode(RoundingMode.CEILING);

		int mNumFrames = this.getNoOfFrames(path);
		int mChannels = this.getNoOfChannels(path);
		
		// take the mean of amplitude values across all the channels and convert the
		// signal to mono mode
		double[] meanBuffer = new double[mNumFrames];
		for (int q = 0; q < mNumFrames; q++) {
			double frameVal = 0;
			for (int p = 0; p < mChannels; p++) {
				frameVal = frameVal + magValueArray[p][q];
			}
			meanBuffer[q] = Double.parseDouble(df.format(frameVal / mChannels));
		}

		return meanBuffer;

	}

}
