package com.jlibrosa.audio;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

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
public class SpleeterTest {
	
	static JLibrosa jLibrosa = null;

	public static void main(String[] args) throws IOException, WavFileException, FileFormatNotSupportedException {
		String audioFilePath = "audioFiles/AClassicEducation.wav";
		int defaultSampleRate = -1;		//-1 value implies the method to use default sample rate
		int defaultAudioDuration = 20;	//-1 value implies the method to process complete audio duration
		
		jLibrosa = new JLibrosa();

		/* To read the magnitude values of audio files - equivalent to librosa.load('../audioFiles/1995-1826-0003.wav', sr=None) function */
		
		float [][] stereoFeatureValues = jLibrosa.loadAndReadStereo(audioFilePath, defaultSampleRate, defaultAudioDuration);
		
		float [][] stereoTransposeFeatValues = transposeMatrix(stereoFeatureValues);
		
		float[][][][] array4D = new float[10][10][9][8];
		
		Double[][][][] stftValues = _stft(stereoTransposeFeatValues);
		
			}
	
	
	private static Double[][][][] _stft(float[][] stereoMatrix) {
		
		int N=4096;
		int H=1024;
		int sampleRate = 44100;
		
		ArrayList<Complex [][]> stftValuesList = new ArrayList<Complex [][]>();
		
		for (int i=0;i<stereoMatrix[0].length;i++) {
			float[] doubleStream = getColumnFromMatrix(stereoMatrix, i);
			
			Complex[][] stftComplexValues = jLibrosa.generateSTFTFeatures(doubleStream, sampleRate, 40,4096,128,1024);
			Complex[][] transposedSTFTComplexValues = transposeMatrix(stftComplexValues);
			stftValuesList.add(transposedSTFTComplexValues);
		}
		
		int segmentLen = 512;
		
		Double [][][] stft3DMatrixValues = gen3DMatrixFrom2D(stftValuesList, segmentLen);
		
		
		int splitValue = (stft3DMatrixValues.length + segmentLen - 1)/segmentLen;
		
		Double [][][][] stft4DMatrixValues = gen4DMatrixFrom3D(stft3DMatrixValues, splitValue, segmentLen);
		
		return stft4DMatrixValues;
		
	}
	
	
	private static Double [][][][] gen4DMatrixFrom3D(Double [][][] stft3DMatrixValues, int splitValue, int segmentLen){
		
		
		int yVal = 1024;
		int zVal = stft3DMatrixValues[0][0].length;
		
		Double[][][][] stft4DMatrixValues = new Double [splitValue][segmentLen][yVal][zVal];
		
		
		
		for(int p=0;p<splitValue;p++) {
			for(int q=0;q<segmentLen;q++) {
				int retInd = (p*segmentLen)+q;
				for(int r=0;r<yVal;r++) {
					for(int s=0;s<zVal;s++) {
						stft4DMatrixValues[p][q][r][s] = stft3DMatrixValues[retInd][r][s];
					}
				}
			}
		}
		
		return stft4DMatrixValues;
		
	}
	
	
	private static Double[][][] gen3DMatrixFrom2D(ArrayList<Complex[][]> mat2DValuesList, int segmentLen){
		
		int padSize = computePadSize(mat2DValuesList.get(0).length, segmentLen);
		int matrixXLen = mat2DValuesList.get(0).length + padSize;
		
		Double [][][] complex3DMatrix = new Double[matrixXLen][mat2DValuesList.get(0)[0].length][mat2DValuesList.size()];
		
		
		for(int k=0;k<mat2DValuesList.size();k++) {
			Complex[][] mat2DValues = mat2DValuesList.get(k);
			for(int i=0;i<matrixXLen;i++) {
				for(int j=0;j<mat2DValues[0].length;j++) {
					double value = 0.0;
					if(i<mat2DValues.length) {
						value = mat2DValues[i][j].abs();
					}
					
					complex3DMatrix[i][j][k]=value;
				}
			}
		}
		return complex3DMatrix;
		
	}
	
	
	private static int computePadSize(int currentMatrixLen, int segmentLen) {
		int tensorSize = currentMatrixLen % segmentLen;
		int padSize = segmentLen - tensorSize;
		return padSize;
	}
	
	
	private static float [] getColumnFromMatrix(float[][] floatMatrix, int column) {
		double[] doubleStream =  IntStream.range(0, floatMatrix.length).mapToDouble(i -> floatMatrix[i][column]).toArray();
		float[] floatArray = new float[doubleStream.length];
		for (int i = 0 ; i < doubleStream.length; i++)
		{
		    floatArray[i] = (float) doubleStream[i];
		}
		return floatArray;
	}
	
	
	public static Complex [][] transposeMatrix(Complex[][] matrix){
	    int m = matrix.length;
	    int n = matrix[0].length;

	    Complex[][] transposedMatrix = new Complex[n][m];

	    for(int x = 0; x < n; x++) {
	        for(int y = 0; y < m; y++) {
	            transposedMatrix[x][y] = matrix[y][x];
	        }
	    }

	    return transposedMatrix;
	}
	
	
	public static float[][] transposeMatrix(float[][] matrix){
	    int m = matrix.length;
	    int n = matrix[0].length;

	    float[][] transposedMatrix = new float[n][m];

	    for(int x = 0; x < n; x++) {
	        for(int y = 0; y < m; y++) {
	            transposedMatrix[x][y] = matrix[y][x];
	        }
	    }

	    return transposedMatrix;
	}
	

}
