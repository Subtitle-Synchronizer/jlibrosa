package com.jlibrosa.audio;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.analysis.function.Log;
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
public class JLibrosaTest_1 {

	public static void main(String[] args) throws IOException, WavFileException, FileFormatNotSupportedException {
		String audioFilePath = "audioFiles/AClassicEducation.wav";
		int defaultSampleRate = -1;		//-1 value implies the method to use default sample rate
		int defaultAudioDuration = 5;	//-1 value implies the method to process complete audio duration
		
		JLibrosa jLibrosa = new JLibrosa();

		/* To read the magnitude values of audio files - equivalent to librosa.load('../audioFiles/1995-1826-0003.wav', sr=None) function */
		
		float audioFeatureValues [] = jLibrosa.loadAndReadWithOffset(audioFilePath, defaultSampleRate, defaultAudioDuration, 10);
		
		float audioFeatureValues1 [] = jLibrosa.loadAndReadWithOffset(audioFilePath, defaultSampleRate, 15, 0);
		
		ArrayList<Float> audioFeatureValuesList = jLibrosa.loadAndReadAsList(audioFilePath, defaultSampleRate, defaultAudioDuration);
		
		//writeToFile(audioFeatureValues);
		
		
		for(int i=0;i<10;i++) {
			System.out.printf("%.6f%n", audioFeatureValues[i]);
		}
		
		
		/* To read the no of frames present in audio file*/
		int nNoOfFrames = jLibrosa.getNoOfFrames();
		
		
		/* To read sample rate of audio file */
		int sampleRate = jLibrosa.getSampleRate();
		
		/* To read number of channels in audio file */
		int noOfChannels = jLibrosa.getNoOfChannels();
		
		Complex[][] stftComplexValues = jLibrosa.generateSTFTFeaturesWithPadOption(audioFeatureValues, sampleRate, 40, true);
		
		
		
		float[] invSTFTValues = jLibrosa.generateInvSTFTFeatures(stftComplexValues, sampleRate, 40);
		

		
		
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
		
		Complex[][] stftComplexValues1 = jLibrosa.generateSTFTFeatures(audioFeatureValues, sampleRate, 40);
		
		
		
		float[] invSTFTValues1 = jLibrosa.generateInvSTFTFeatures(stftComplexValues, sampleRate, 40);
		
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
	
	
	
	public static void encodePcmToMp3(byte[] pcm) {
//      LameEncoder encoder = new LameEncoder(new javax.sound.sampled.AudioFormat(44100.0f, 16, 2, true, false), 256, MPEGMode.STEREO, Lame.QUALITY_HIGHEST, false);

  //fast  gmm2  LameEncoder encoder = new LameEncoder(new javax.sound.sampled.AudioFormat(88200.0f, 16, 2, true, false), 256, MPEGMode.STEREO, Lame.QUALITY_HIGHEST, false);

		
		/*
     LameEncoder encoder = new LameEncoder(new javax.sound.sampled.AudioFormat(44100.0f, 16, 2, true, false), 256, MPEGMode.STEREO, Lame.QUALITY_HIGHEST, false);


      ByteArrayOutputStream mp3 = new ByteArrayOutputStream();
      byte[] buffer = new byte[encoder.getPCMBufferSize()];

      int bytesToTransfer = Math.min(buffer.length, pcm.length);
      int bytesWritten;
      int currentPcmPosition = 0;

       while (0 < (bytesWritten = encoder.encodeBuffer(pcm, currentPcmPosition, bytesToTransfer, buffer))) {
          currentPcmPosition += bytesToTransfer;
          bytesToTransfer = Math.min(buffer.length, pcm.length - currentPcmPosition);
          //Log.e("logmessage", "current position: " + currentPcmPosition);
          mp3.write(buffer, 0, bytesWritten);
      }

      encoder.close();

      File file = new File("gtmm5.mp3");
      if (!file.exists()) {
          try {
              file.createNewFile();
          } catch (IOException e) {
              e.printStackTrace();
             // Log.e("logmessage", "cannot create file");
          }

          FileOutputStream stream = null;
          try {
              stream = new FileOutputStream("gtmm5.mp3");
              stream.write(mp3.toByteArray());
          } catch (FileNotFoundException e) {
              e.printStackTrace();
          } catch (IOException e) {
              e.printStackTrace();
          }

      }
      
      */
   //   return mp3.toByteArray();
  }
	
	
	public static void writeToFile(float[] array) throws IOException 
	{ 
	    
	        byte [] consByteArray = new byte[4*array.length];
	        
	        
		        for(int i=0;i<array.length;i++) {
		        	byte[] byteArray = ByteBuffer.allocate(4).putFloat(array[i] * 32767).array();
		        	for(int k=0;k<byteArray.length;k++) {
		        		consByteArray[i*4+k]=byteArray[k];
		        	}
		        }
		        
		        encodePcmToMp3(consByteArray);
	        
	} 
	

}
