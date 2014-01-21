/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package simulatorGUI;

import java.io.ByteArrayOutputStream;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;
import org.apache.commons.codec.binary.Base64;

public class BinaryReaderWriter {
	public static void main(String [ ] args) throws UnsupportedEncodingException, DataFormatException
	{
		// this input will break the fixed buffer method
		double[] centroids = {123.1212234143345453223123123, 28464632322456781.23, 31232342342341234121.0};
		
		// this input will break the fixed buffer method
		//double[] centroids = {123.1212234143345453223123123, 28464632322456781.23, 31.0};
		BinaryReaderWriter brw = new BinaryReaderWriter();
		String output = brw.compressCentroids(centroids);
		brw.decompressCentroids(output);
	}
	void decompressCentroids(String encoded) throws DataFormatException{
		byte[] binArray = Base64.decodeBase64(encoded);
		
		/*
		// This block of code is the fixed buffer version
		//
System.out.println("binArray length " + binArray.length);
		Inflater deCompressor = new Inflater();
        deCompressor.setInput(binArray, 0, binArray.length);
        byte[] decompressed = new byte[1024];
		int decompressedLength = deCompressor.inflate(decompressed);
        deCompressor.end();
System.out.println("decompressedLength = " + decompressedLength);
		byte[] decompressedData = new byte[decompressedLength];
		for(int i=0;i<decompressedLength;i++){
			decompressedData[i] = decompressed[i];
		}
		*/
		
		
		// This block of code is the variable buffer version
		//
		ByteArrayOutputStream bos = new ByteArrayOutputStream(binArray.length);
		Inflater deCompressor = new Inflater();
        deCompressor.setInput(binArray, 0, binArray.length);
        byte[] decompressed = new byte[1024];
		while (!deCompressor.finished()) {
			int decompressedLength = deCompressor.inflate(decompressed);
			bos.write(decompressed, 0, decompressedLength);
		}
        deCompressor.end();
		byte[] decompressedData = bos.toByteArray();
		
		
		ByteBuffer bb = ByteBuffer.wrap(decompressedData);
		bb.order(ByteOrder.LITTLE_ENDIAN);
System.out.println("decompressedData length = " + decompressedData.length);
		double[] doubleValues = new double[decompressedData.length / 8];
		for (int i = 0; i< doubleValues.length; i++){
			doubleValues[i] = bb.getDouble(i * 8);
		}

		for(double dbl : doubleValues){
			System.out.println(dbl);
		}	
	}
	
	String compressCentroids(double[] centroids){
		byte[] cinput = new byte[centroids.length * 8];
		ByteBuffer buf = ByteBuffer.wrap(cinput);
		buf.order(ByteOrder.LITTLE_ENDIAN);
		for (double cent : centroids){
			buf.putDouble(cent);
		}
		
		byte[] input = buf.array();
System.out.println("raw length = " + input.length);
		byte[] output = new byte[input.length * 2];
		Deflater compresser = new Deflater();
		compresser.setInput(input);
		compresser.finish();
		int compressedLength = compresser.deflate(output);
		compresser.end();
System.out.println("Compressed length = " + compressedLength);
		byte[] compressed = new byte[compressedLength];
		for(int i = 0; i < compressedLength; i++){
			compressed[i] = output[i];
		}
		
		String decrypted = Base64.encodeBase64String(compressed);
		return decrypted;
	}
}
