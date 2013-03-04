import java.awt.Point;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Serializable;
import java.util.Arrays;

import javax.imageio.ImageIO;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;

// GWTGrid implements serializable so that it can be written out
public class GWTGrid {

	int size, numScales, numOrientations;
	double[][][] freqKernels;
	DoubleFFT_2D fft;

	double[][][][] magnitudeResps;
	double[][][][] phaseResps;

	// construct a gwt grid from an image
	/*
	 * 4d freq kernels
	 * first dimension: scale
	 * second dimension: orientation
	 * third dimension: x
	 * fourth dimension: y
	 */
	public GWTGrid(double[] pixels, int size, double[][][] freqKernels) {
		numScales = freqKernels.length;
		numOrientations = freqKernels[0].length;
		this.size = size;
		this.freqKernels = freqKernels;

		// fft shift the image
		fft = new DoubleFFT_2D(size, size);
		// interleave our vectors since that is what this implementation of fft needs
		double[] realVector = pixels;
		double[] imaginaryVector = new double[pixels.length];
		Arrays.fill(imaginaryVector, 0.0);
		double[] complexVector = MathUtils.interleave(realVector, imaginaryVector);
		fft.complexForward(complexVector);
		// get phase and magnitude responses using fft shifted image
		getPhaseAndMagResps(complexVector);
	}


	// setup the freq kernels
	public static double[][][] genFreqKernel(int size, int numScales, int numOrientations, double sigma) {
		double[][][] freqKernels = new double[numScales][numOrientations][size*size*2];

		int xyResL = size;
		int xHalfResL = (int)((double)size/2.0);
		int yHalfResL = xHalfResL;
		double kxFactor = 2.0*Math.PI/(double)xyResL;
		double kyFactor = 2.0*Math.PI/(double)xyResL;

//		[tx,ty] = meshgrid(-xHalfResL:xHalfResL-1,-yHalfResL:yHalfResL-1);
//		tx = kxFactor*tx;
//		ty = kyFactor*(-ty);

		for(int scale = 0; scale < numScales; scale++) {
			double k0 = (Math.PI/2.0)*Math.pow((1.0/Math.sqrt(2.0)), scale);
			for(int orientation = 0; orientation < numOrientations; orientation++) {
				int curX = 0;
				double[][] freqKernel = new double[size][size];
				double kA = Math.PI*(double)orientation/(double)numOrientations;
				double k0X = k0*Math.cos(kA);
				double k0Y = k0*Math.sin(kA);
				for(int x = -xHalfResL; x <= xHalfResL-1; x++) {
					int curY = 0;
					for(int y = -yHalfResL; y <= yHalfResL-1; y++) {
						double tx = kxFactor*x;
						double ty = kyFactor*(-y);
						//              2*pi*(exp(-(Sigma/k0)^2/2*((k0X-tx).^2+(k0Y-ty).^2))-exp(-(Sigma/k0)^2/2*(k0^2+tx.^2+ty.^2)));
						freqKernel[curX][curY] = 2.0*Math.PI*(Math.exp(-Math.pow((sigma/k0),2.0)/2.0*(Math.pow((k0X-tx),2.0)+Math.pow((k0Y-ty),2.0)))-Math.exp(-Math.pow((sigma/k0),2.0)/2.0*(Math.pow(k0,2.0)+Math.pow(tx,2.0)+Math.pow(ty,2.0))));
						curY++;
					}
					curX++;
				}
				// fft shift freqkernels
				freqKernel = MathUtils.shiftOrigin(freqKernel);

				// put freqKernel into complex representation
				// first make freqKernel 1D. then interleave with imag 0s
				double[] zeros = new double[size*size];
				Arrays.fill(zeros, 0.0);

				double[] freqKernel1D = new double[size*size];
				int freqKernel1DIdx = 0;
				for(int tmpY = 0; tmpY < size; tmpY++) {
					for(int tmpX = 0; tmpX < size; tmpX++) {
						freqKernel1D[freqKernel1DIdx++] = freqKernel[tmpX][tmpY];
					}
				}

				double[] complexFreqKernel = MathUtils.interleave(freqKernel1D, zeros);

				// copy freqKernel into freqKernels
				for(int tmpX = 0; tmpX < size*size*2; tmpX++) {
					freqKernels[scale][orientation][tmpX] = complexFreqKernel[tmpX];
				}
			}
		}
		return freqKernels;
	}

	private void getPhaseAndMagResps(double[] imgFreq) {
		magnitudeResps = new double[numScales][numOrientations][size][size];
		phaseResps = new double[numScales][numOrientations][size][size];

		// for all orientations and scales, multiply imgFreq with freqKernel
		for(int scale = 0; scale < numScales; scale++) {
			for(int orientation = 0; orientation < numOrientations; orientation++) {
				double[] product = MathUtils.complexMultiply(freqKernels[scale][orientation], imgFreq);
				// now inverse transform the product
				fft.complexInverse(product, true);

				// get rid of our complex numbers and compute magnitude (abs) and phase (angle) repsones
				int curX = 0;
				int curY = 0;				
				for(int i = 0; i < size*size*2; i += 2) {
					double realPart = product[i];
					double imagPart = product[i+1];

					magnitudeResps[scale][orientation][curX][curY] = Math.sqrt(Math.pow(realPart, 2.0)+Math.pow(imagPart, 2.0));
					// this is not a perfect translation. possible source of error
					phaseResps[scale][orientation][curX][curY] = Math.atan2(imagPart, realPart) + Math.PI;
					curX++;
					if(curX == size) {
						curX = 0;
						curY++;
					}
				}
			}
		}
	}

	// gets magnitude response for index x,y by stacking the scale/orientation responses
	public double[] getMagnitudeResp(int x, int y) {
		double[] resp = new double[numScales*numOrientations];

		int respIdx = 0;
		for(int scale = 0; scale < numScales; scale++) {
			for(int orientation = 0; orientation < numOrientations; orientation++) {
				resp[respIdx++] = magnitudeResps[scale][orientation][x][y];
			}
		}
		return resp;
	}

	// gets phase response for index x,y by stacking the scale/orientation responses
	public double[] getPhaseResp(int x, int y) {
		double[] resp = new double[numScales*numOrientations];

		int respIdx = 0;
		for(int scale = 0; scale < numScales; scale++) {
			for(int orientation = 0; orientation < numOrientations; orientation++) {
				resp[respIdx++] = phaseResps[scale][orientation][x][y];
			}
		}

		return resp;
	}

	public static void main(String[] args) {
		double sigma = 1.0*Math.PI;
		int numOrientations = 8;
		int numScales = 5;	

		try {
			// read in images -- must be same size
			BufferedImage img = ImageIO.read(new File("corner1.jpg"));
			BufferedImage img2 = ImageIO.read(new File("corner2.jpg"));			
			int size = img.getWidth();
			double[] pixels = ImageUtils.RGBtoGrayDouble(ImageUtils.getPixels(img));
			double[] pixels2 = ImageUtils.RGBtoGrayDouble(ImageUtils.getPixels(img2));							

			long startTime = System.currentTimeMillis();
			// generate frequency kernel
			double[][][] freqKernels = GWTGrid.genFreqKernel(size, numScales, numOrientations, sigma);

			// create gwt grids
			GWTGrid gwt = new GWTGrid(pixels, size, freqKernels);
			GWTGrid gwt2 = new GWTGrid(pixels2, size, freqKernels);

			// point in img1 to find similarity of in img2
			int xCmpPoint = 64;
			int yCmpPoint = 58;

			double maxSim = 0;
			Point maxSimIdx = null; 
			for(int i = 0; i < img2.getWidth(); i++) {
				for(int j = 0; j < img2.getHeight(); j++) {
					// compare img1's mag/phase @ xCmpPoint/yCmpPoint to every location on img2
					double similarity = MathUtils.compareTwoJets(gwt.getMagnitudeResp(xCmpPoint, yCmpPoint), gwt2.getMagnitudeResp(i,j), 
					gwt.getPhaseResp(xCmpPoint, yCmpPoint), gwt2.getPhaseResp(i,j));
					if(similarity > maxSim) {
						maxSim = similarity;
						maxSimIdx = new Point(i, j);
					}
				}
			}
			long elapsedTime = System.currentTimeMillis() - startTime;
			System.out.println("maxSim=" + maxSim + ", maxSimIdx=" + maxSimIdx.x + "," + maxSimIdx.y);
			System.out.println("elapsed time=" + elapsedTime + " ms");
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}

}
