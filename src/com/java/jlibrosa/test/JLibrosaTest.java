package com.java.jlibrosa.test;

import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;

import com.java.audio.librosafeatures.WavFile;
import com.java.audio.librosafeatures.WavFileException;
import com.java.audio.process.MFCC;

public class JLibrosaTest {

	public static void main(String[] args) throws IOException, WavFileException {
		// TODO Auto-generated method stub
		
		int mNumFrames;
		int mSampleRate;
		int mChannels;
		
		File sourceFile = new File("audioFiles/100263-2-0-126.wav");
        WavFile wavFile = null;
        
            wavFile = WavFile.openWavFile(sourceFile);
            mNumFrames = (int) (wavFile.getNumFrames());
            mSampleRate = (int) wavFile.getSampleRate();
            mChannels = wavFile.getNumChannels();

            double[][] buffer = new double[mChannels][mNumFrames];
            int frameOffset = 0;
            int loopCounter = ((mNumFrames * mChannels)/4096) + 1;
            for (int i = 0; i < loopCounter; i++) {
                frameOffset = wavFile.readFrames(buffer, mNumFrames, frameOffset);
            }


            DecimalFormat df = new DecimalFormat("#.#####");
            df.setRoundingMode(RoundingMode.CEILING);

            double [] meanBuffer = new double[mNumFrames];
            for(int q=0;q<mNumFrames;q++){
                double frameVal = 0;
                for(int p=0;p<mChannels;p++){
                    frameVal = frameVal + buffer[p][q];
                }
                    meanBuffer[q]=Double.parseDouble(df.format(frameVal/mChannels));
            }


            //MFCC java library.
            MFCC mfccConvert = new MFCC();
            mfccConvert.setSampleRate(mSampleRate);
            int nMFCC = 40;
            mfccConvert.setN_mfcc(nMFCC);
            float[] mfccInput = mfccConvert.process(meanBuffer);

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

            //code to take the mean of mfcc values across the rows such that
            //[nMFCC x nFFT] matrix would be converted into
            //[nMFCC x 1] dimension - which would act as an input to tflite model
            float [] meanMFCCValues = new float[nMFCC];
            for(int p=0;p<nMFCC;p++){
                double fftValAcrossRow = 0;
                for(int q=0;q<nFFT;q++){
                    fftValAcrossRow = fftValAcrossRow + mfccValues[p][q];
                }
                double fftMeanValAcrossRow = fftValAcrossRow/nFFT;
                meanMFCCValues[p] = (float) fftMeanValAcrossRow;
            }

       
        
	}

}
