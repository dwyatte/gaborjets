import java.io.File;
import java.util.HashSet;
import java.util.Iterator;

public class MathUtils {
	
	// construct an array 0:len-1
	public static int[] arrayConstructor(int len) {
		int[] A = new int[len];
		for(int i = 0; i < len; i++) {
			A[i] = i;
		}
		return A;
	}
	
	// tile an an array numTimes
	public static int[] repeatingArray(int[] baseArray, int numTimes) {
		int[] newArray = new int[baseArray.length * numTimes];
		
		for(int i = 0; i < newArray.length; i++) {
			newArray[i] = baseArray[i % baseArray.length];
		}
		
		return newArray;
	}
	
	
	// shuffle an array
	public static int[] shuffle(int[] A) {
		for(int i = 0; i < A.length; i++) {
			int swapInd = (int)(Math.random()*(double)A.length);
			int tmp = A[i];
			A[i] = A[swapInd];
			A[swapInd] = tmp;
		}
		
		return A;
	}
	public static File[] shuffle(File[] A) {
		for(int i = 0; i < A.length; i++) {
			int swapInd = (int)(Math.random()*(double)A.length);
			File tmp = A[i];
			A[i] = A[swapInd];
			A[swapInd] = tmp;
		}
		
		return A;
	}
	
	// reshape an array into a wxh matrix
	public static int[][] reshape(int[] A, int w, int h) {
		int[][] A2D = new int[w][h];
		
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				A2D[i][j] = A[j*w + i];
			}
		}
		
		return A2D;
	}
	
	// reshape an array into a wxh matrix
	public static double[][] reshape(double[] A, int w, int h) {
		double[][] A2D = new double[w][h];
		
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				A2D[i][j] = A[j*w + i];
			}
		}
		
		return A2D;
	}
	
	// creates a new vector by interleaving the two that are passed in
	public static double[] interleave(double[] v1, double[] v2) {
		double[] newV = new double[v1.length + v2.length];
		
		int oldVIdx = 0;
		for(int i = 0; i < newV.length; i += 2) {
			newV[i] = v1[oldVIdx];
			newV[i+1] = v2[oldVIdx];
			oldVIdx++;
		}
		
		return newV;
	}
	
	
	// this method does complex multiplication across 2 vectors, that are organized
	// in the manner v = [real1 imag1 real2 imag2 ... ...]
	public static double[] complexMultiply(double[] complexV1, double[] complexV2) {
        //return new Complex(x*w.real()-y*w.imag(),x*w.imag()+y*w.real());
		double[] productVec = new double[complexV1.length];
		
		for(int i = 0; i < complexV1.length-1; i+=2) {
			// get the real numbers
			double real1 = complexV1[i];
			double real2 = complexV2[i];
			// get the imaginary numbers
			double imag1 = complexV1[i+1];
			double imag2 = complexV2[i+1];
			
			double realProd = real1*real2-imag1*imag2;
			double imagProd = real1*imag2+imag1*real2;
			// do the math
			productVec[i] = realProd;
			productVec[i+1] = imagProd;
		}
		
		return productVec;
	}
	
	
	// translation of GWT similarity function
	public static double compareTwoJets(double[] magnitudeResps1, double[] magnitudeResps2, double[] phaseResps1, double[] phaseResps2) {
		double[] magnitudeProds = new double[magnitudeResps1.length];
		double[] phaseDiffs = new double[phaseResps1.length];
		double[] magnitudeResps1Squared = new double[magnitudeResps1.length];
		double[] magnitudeResps2Squared = new double[magnitudeResps2.length];
		double numerator = 0;
		double sumMagnitudeResps1 = 0;
		double sumMagnitudeResps2 = 0;
		double denominator = 0;
		
		// product of mag resps
		// diff of phase resps
		// squared mag resps
		for(int i = 0; i < magnitudeResps1.length; i++) {
			magnitudeProds[i] = magnitudeResps1[i]*magnitudeResps2[i];
			phaseDiffs[i] = phaseResps1[i]-phaseResps2[i];
			magnitudeResps1Squared[i] = Math.pow(magnitudeResps1[i], 2.0);
			magnitudeResps2Squared[i] = Math.pow(magnitudeResps2[i], 2.0);
		}
		
		// for the numerator, sum of magnitudes * cos of phases
		// for the denominator, sum the magnitude resps
		for(int i = 0; i < magnitudeResps1.length; i++) {
			numerator += magnitudeProds[i]*Math.cos(phaseDiffs[i]);
			sumMagnitudeResps1 += magnitudeResps1Squared[i];
			sumMagnitudeResps2 += magnitudeResps2Squared[i];
		}
		denominator = sumMagnitudeResps1 * sumMagnitudeResps2;
		
		return numerator/Math.sqrt(denominator);
	}	
	
  // Method to shift the wavenumber origin and
  // place it at the center for a more visually
  // pleasing display.  Must be applied
  // separately to the real part, the imaginary
  // part, and the amplitude spectrum for a wave-
  // number spectrum.
  // Courtesy of R.G.Baldwin
  public static double[][] shiftOrigin(double[][] data){
    int numberOfRows = data.length;
    int numberOfCols = data[0].length;
    int newRows;
    int newCols;
   
    double[][] output =
          new double[numberOfRows][numberOfCols];
         
    //Must treat the data differently when the
    // dimension is odd than when it is even.
   
    if(numberOfRows != 0){//odd
      newRows = numberOfRows +
                            (numberOfRows + 1)/2;
    }else{//even
      newRows = numberOfRows + numberOfRows/2;
    }//end else
   
    if(numberOfCols != 0){//odd
      newCols = numberOfCols +
                            (numberOfCols + 1)/2;
    }else{//even
      newCols = numberOfCols + numberOfCols/2;
    }//end else
   
    //Create a temporary working array.
    double[][] temp =
                    new double[newRows][newCols];
                   
    //Copy input data into the working array.
    for(int row = 0;row < numberOfRows;row++){
      for(int col = 0;col < numberOfCols;col++){
        temp[row][col] = data[row][col];
      }//col loop
    }//row loop
   
    //Do the horizontal shift first
    if(numberOfCols != 0){//shift for odd
      //Slide leftmost (numberOfCols+1)/2 columns
      // to the right by numberOfCols columns
      for(int row = 0;row < numberOfRows;row++){
        for(int col = 0;
                 col < (numberOfCols+1)/2;col++){
          temp[row][col + numberOfCols] =
                                  temp[row][col];
        }//col loop
      }//row loop
      //Now slide everything back to the left by
      // (numberOfCols+1)/2 columns
      for(int row = 0;row < numberOfRows;row++){
        for(int col = 0;
                       col < numberOfCols;col++){
          temp[row][col] =
             temp[row][col+(numberOfCols + 1)/2];
        }//col loop
      }//row loop
     
    }else{//shift for even
      //Slide leftmost (numberOfCols/2) columns
      // to the right by numberOfCols columns.
      for(int row = 0;row < numberOfRows;row++){
        for(int col = 0;
                     col < numberOfCols/2;col++){
          temp[row][col + numberOfCols] =
                                  temp[row][col];
        }//col loop
      }//row loop
     
      //Now slide everything back to the left by
      // numberOfCols/2 columns
      for(int row = 0;row < numberOfRows;row++){
        for(int col = 0;
                       col < numberOfCols;col++){
          temp[row][col] =
                 temp[row][col + numberOfCols/2];
        }//col loop
      }//row loop
    }//end else
    //Now do the vertical shift
    if(numberOfRows != 0){//shift for odd
      //Slide topmost (numberOfRows+1)/2 rows
      // down by numberOfRows rows.
      for(int col = 0;col < numberOfCols;col++){
        for(int row = 0;
                 row < (numberOfRows+1)/2;row++){
          temp[row + numberOfRows][col] =
                                  temp[row][col];
        }//row loop
      }//col loop
      //Now slide everything back up by
      // (numberOfRows+1)/2 rows.
      for(int col = 0;col < numberOfCols;col++){
        for(int row = 0;
                       row < numberOfRows;row++){
          temp[row][col] =
             temp[row+(numberOfRows + 1)/2][col];
        }//row loop
      }//col loop
     
    }else{//shift for even
      //Slide topmost (numberOfRows/2) rows down
      // by numberOfRows rows
      for(int col = 0;col < numberOfCols;col++){
        for(int row = 0;
                     row < numberOfRows/2;row++){
          temp[row + numberOfRows][col] =
                                  temp[row][col];
        }//row loop
      }//col loop
     
      //Now slide everything back up by
      // numberOfRows/2 rows.
      for(int col = 0;col < numberOfCols;col++){
        for(int row = 0;
                       row < numberOfRows;row++){
          temp[row][col] =
                 temp[row + numberOfRows/2][col];
        }//row loop
      }//col loop
    }//end else
   
    //Shifting of the origin is complete.  Copy
    // the rearranged data from temp to output
    // array.
    for(int row = 0;row < numberOfRows;row++){
      for(int col = 0;col < numberOfCols;col++){
        output[row][col] = temp[row][col];
      }//col loop
    }//row loop
    return output;
  }//end shiftOrigin method
}
