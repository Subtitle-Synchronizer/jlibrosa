package com.jlibrosa.audio.exception;


/**
 * 
 * This Class is an custom exception class to throw exceptions when unsupported files are provided as input for processing.
 * 
 * @author vvasanth
 *
 */

public class FileFormatNotSupportedException extends Exception {
	
	public FileFormatNotSupportedException(String message) {
		super(message);
	}

}
