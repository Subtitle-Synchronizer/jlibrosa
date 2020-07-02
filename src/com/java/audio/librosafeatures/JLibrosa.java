package com.java.audio.librosafeatures;

import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import com.java.audio.process.MFCC;

public class JLibrosa {
	
	private String path;
	private long sr;
	private boolean mono;
	
	
	private int BUFFER_SIZE = 4096;
	private int mNumFrames;
	private int mSampleRate;
	private int mChannels;
	
	public JLibrosa() {
		
	}
	
	public double [][] loadAndReadAcrossChannels(String path, int sr, int readDurationInSec) throws IOException, WavFileException{
		double[][] magValues = readMagnitudeValuesFromFile(path, sr, readDurationInSec);
		return magValues;
	}
	
	
	private double[][] readMagnitudeValuesFromFile(String path, int sr, int readDurationInSeconds) throws IOException, WavFileException{
		
		File sourceFile = new File(path);
        WavFile wavFile = null;

        wavFile = WavFile.openWavFile(sourceFile);
        mNumFrames = (int) (wavFile.getNumFrames());
        mSampleRate = (int) wavFile.getSampleRate();
        mChannels = wavFile.getNumChannels();

        if(readDurationInSeconds != -1) {
        	mNumFrames = readDurationInSeconds * mSampleRate;
        }
        
        if(sr!=-1) {
        	mSampleRate = sr; 
        }
        
        
        //Read the magnitude values across both the channels and save them as part of 
        //multi-dimensional array
        double[][] buffer = new double[mChannels][mNumFrames];
        int frameOffset = 0;
        int loopCounter = ((mNumFrames * mChannels)/BUFFER_SIZE) + 1;
        for (int i = 0; i < loopCounter; i++) {
            frameOffset = wavFile.readFrames(buffer, mNumFrames, frameOffset);
        }
        
        return buffer;
        
	}
	
	
	public double [][] generateMFCCFeatures(double[] magValues, int sr, int nMFCC){
		
		MFCC mfccConvert = new MFCC();
        mfccConvert.setSampleRate(sr);
        mfccConvert.setN_mfcc(nMFCC);
        float[] mfccInput = mfccConvert.process(magValues);

        int nFFT = mfccInput.length/nMFCC;
        double [][] mfccValues = new double[nMFCC][nFFT];

        //loop to convert the mfcc values into multi-dimensional array
        for(int i=0;i<nFFT;i++){
            int indexCounter = i * nMFCC;
            int rowIndexValue = i%nFFT;
            for(int j=0;j<nMFCC;j++){
                mfccValues[j][rowIndexValue]=mfccInput[indexCounter];
                indexCounter++;
            }
        }
		
		return mfccValues;
		
	}
	
	
	public double [][] generateSTFTFeatures(double[] magValues, int sr, int nMFCC){
		MFCC mfccConvert = new MFCC();
		mfccConvert.setSampleRate(sr);
		mfccConvert.setN_mfcc(nMFCC);
		double [][] stftValues = mfccConvert.stftMagSpec(magValues);
		return stftValues;
		
	}
	
	public double [] loadAndRead(String path, int sr, int readDurationInSeconds) throws IOException, WavFileException {
	
		double[][] magValueArray = readMagnitudeValuesFromFile(path, sr, readDurationInSeconds);
		
        DecimalFormat df = new DecimalFormat("#.#####");
        df.setRoundingMode(RoundingMode.CEILING);

        //take the mean of amplitude values across all the channels and convert the signal to mono mode
        double [] meanBuffer = new double[mNumFrames];
        for(int q=0;q<mNumFrames;q++){
            double frameVal = 0;
            for(int p=0;p<mChannels;p++){
                frameVal = frameVal + magValueArray[p][q];
            }
                meanBuffer[q]=Double.parseDouble(df.format(frameVal/mChannels));
        }

        return meanBuffer;
		
	}
	
	
	

}
