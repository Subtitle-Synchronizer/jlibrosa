package com.jlibrosa.audio;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;


public class JLibrosaMatteoTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub

		JLibrosa jLibrosa = new JLibrosa();

		JSONParser jsonParser = new JSONParser();
		 try (FileReader reader = new FileReader("/Users/vishrud/Downloads/example_data_experiment.json"))
	        {
	            //Read JSON file
	            Object obj = jsonParser.parse(reader);
	 
	            JSONObject signalObject = (JSONObject) obj;
	            
	            JSONArray sigArr = (JSONArray) signalObject.get("raw_signal");
	            
	            float [] sigFltArr = new float[sigArr.size()];
	            
	            for(int i=0;i<sigArr.size();i++) {
	            	float val = Float.parseFloat((String.valueOf(sigArr.get(i))));
	            	sigFltArr[i] =  val/32768.0f;
	            }
	            
	        	
	    		float [][] melSpectrogram = jLibrosa.generateMelSpectroGram(sigFltArr, 22050, 1024, 128, 128);
	    		
	            
	            System.out.println(1000);
	             
	           
	        } catch (FileNotFoundException e) {
	            e.printStackTrace();
	        } catch (IOException e) {
	            e.printStackTrace();
	        } catch (ParseException e) {
	            e.printStackTrace();
	        }
	}

}
