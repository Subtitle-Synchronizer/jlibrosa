# jlibrosa

jLibrosa has been conceptualized to build as an equivalent of Python's librosa library. 

Librosa library is widely used to process audio files to generate various values such as magnitude, stft, istft, mfcc etc. These concepts are widely employed in building prediction systems associated with audio form of data.

jLibrosa helps in performing similar kind of processing and to generate the aforementioned features in java. 

jLibrosa could be used in Android environment, so Deep Learning models built in Tensorflow and ported to TFLite could be leveraged in Android environment for building mobile apps.

Refer this blogpost for understanding how to integrate jlibrosa with Android apps for making predictions using TFLite models.

Processed values of audio files generated from jLibrosa would be very similar to the respective values from Python librosa files. Refer the JlibrosaTest.java and the corresponding python notebook to compare the values from both the libraries.