import java.awt.Color;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.MemoryImageSource;
import java.awt.image.PixelGrabber;
import java.util.Vector;

public class ImageUtils {

	/* scale image to preferred dimension
	 * first scale width to preferred width
	 * then scale scaledHeight further to preferred height
	 */
	public static Image scaleImage(Image I, int preferredWidth, int preferredHeight) {
		Image scaledI = I;
		int width = I.getWidth(null);
		int height = I.getHeight(null);
		double scaledWidth, scaledHeight;
		double scaleFactor;

		// scale width
		scaleFactor = (double)preferredWidth/(double)width;
		scaledWidth = width*scaleFactor;
		scaledHeight = height*scaleFactor;
		// scale height from width scaling
		scaleFactor = (double)preferredHeight/(double)scaledHeight;
		scaledWidth = scaledWidth*scaleFactor;
		scaledHeight = scaledHeight*scaleFactor;

		System.out.println("Scaling image to " + (int)scaledWidth + "x" + (int)scaledHeight);
		scaledI = I.getScaledInstance((int)scaledWidth, (int)scaledHeight, Image.SCALE_DEFAULT);

		return scaledI;
	}


	// stitch images
	public static Image stitchImages(Image left, Image right) {
		int leftWidth = left.getWidth(null);
		int leftHeight = left.getHeight(null);
		int rightWidth = right.getWidth(null);
		int rightHeight = right.getHeight(null);
		int taller;
		BufferedImage stitched;
		Graphics stitchedGraphics;

		// find taller image 
		if(leftHeight > rightHeight) {
			taller = leftHeight;
		}
		else {
			taller = rightHeight;
		}

		stitched = new BufferedImage(leftWidth+rightWidth, taller, BufferedImage.TYPE_INT_RGB);
		stitchedGraphics = stitched.getGraphics();
		// left image
		stitchedGraphics.drawImage(left,0,0,null);
		// right image
		stitchedGraphics.drawImage(right, leftWidth, 0, null);

		return stitched;
	}

	// stitches images with padding between
	public static Image stitchImages(Image left, Image right, int padWidth, Color padColor) {
		int leftWidth = left.getWidth(null);
		int leftHeight = left.getHeight(null);
		int rightWidth = right.getWidth(null);
		int rightHeight = right.getHeight(null);
		int taller;
		BufferedImage stitched;
		Graphics stitchedGraphics;

		// find taller image 
		if(leftHeight > rightHeight) {
			taller = leftHeight;
		}
		else {
			taller = rightHeight;
		}

		stitched = new BufferedImage(leftWidth+padWidth+rightWidth, taller, BufferedImage.TYPE_INT_RGB);
		stitchedGraphics = stitched.getGraphics();
		// left image
		stitchedGraphics.drawImage(left,0,0,null);
		// padding
		stitchedGraphics.setColor(padColor);
		stitchedGraphics.fillRect(leftWidth, 0, leftWidth+padWidth, taller);
		// right image
		stitchedGraphics.drawImage(right, leftWidth+padWidth, 0, null);

		return stitched;
	}

	// get pixels from an image
	public static int[] getPixels(Image I) {
		int width = I.getWidth(null); 
		int height = I.getHeight(null); 
		int pixels[] = new int[width * height];
		PixelGrabber pg = new PixelGrabber(I, 0, 0, width, height, pixels, 0, width);
		try { pg.grabPixels(); }
		catch(InterruptedException ie) { System.out.println("Interrupted"); }

		return pixels;
	}
	
	// get pixels from an image
	public static int[] getPixels(Image I, int x, int y, int width, int height) {
		int pixels[] = new int[width * height];
		PixelGrabber pg = new PixelGrabber(I, x,y , width, height, pixels, 0, width);
		try { pg.grabPixels(); }
		catch(InterruptedException ie) { System.out.println("Interrupted"); }
		return pixels;
	}
	
	public static double[] RGBtoGrayDouble(int[] RGBpixels) {
		double[] grayDoublePixels = new double[RGBpixels.length];
		
		for(int i = 0; i < RGBpixels.length; i++) {
			// we are working with grayscale images, so just mask off a channel.
			// blue is the fastest since there is no shifting involved
			int blue  = RGBpixels[i] & 0xff;
			grayDoublePixels[i] = (double)blue;
		}
		
		return grayDoublePixels;
	}
	
	public static int[] RGBtoGray(int[] RGBpixels) {
		int[] grayPixels = new int[RGBpixels.length];
		
		for(int i = 0; i < RGBpixels.length; i++) {
			// we are working with grayscale images, so just mask off a channel.
			// blue is the fastest since there is no shifting involved
			int blue  = RGBpixels[i] & 0xff;
			grayPixels[i] = blue;
		}
		
		return grayPixels;
	}
	

	public static int[] grayDoubleToRGB(double[] grayDoublePixels) {
		double pixel;
		int numClipneg = 0;
		int numClip = 0;
		int[] RGBpixels = new int[grayDoublePixels.length];
		int red, green, blue;
		int alpha = 255;

		for(int i = 0; i < grayDoublePixels.length; i++) {
			pixel = grayDoublePixels[i];
			// clip
			if(pixel < 0) {
				pixel = 0;
				numClipneg++;
			}
			if(pixel > 255) {
				pixel = 255;
				numClip++;
			}
			red = (int)pixel;
			green = (int)pixel;
			blue = (int)pixel;
			// make rgb
			RGBpixels[i] = (alpha << 24) | (red << 16) | (green << 8 ) | blue;
		}

		return RGBpixels;
	}
	
	public static int[] grayToRGB(int[] grayPixels) {
		double pixel;
		int numClipneg = 0;
		int numClip = 0;
		int[] RGBpixels = new int[grayPixels.length];
		int red, green, blue;
		int alpha = 255;

		for(int i = 0; i < grayPixels.length; i++) {
			pixel = grayPixels[i];
			// clip
			if(pixel < 0) {
				pixel = 0;
				numClipneg++;
			}
			if(pixel > 255) {
				pixel = 255;
				numClip++;
			}
			red = (int)pixel;
			green = (int)pixel;
			blue = (int)pixel;
			// make rgb
			RGBpixels[i] = (alpha << 24) | (red << 16) | (green << 8 ) | blue;
		}

		return RGBpixels;
	}


	// create a new java image with pixels
	public static Image makeImageFromPixels(int[] pixels, int width, int height) {
		return Toolkit.getDefaultToolkit().createImage(new MemoryImageSource (width, height, pixels, 0, width));
	}


	public static BufferedImage plotXtoImage(BufferedImage I, int x, int y, Color c) {
		int alpha = 255;
		int red = c.getRed();
		int green = c.getGreen();
		int blue = c.getBlue();
		int pixel = (alpha << 24) | (red << 16) | (green << 8 ) | blue;

		if(x > 0 && y > 0 && x < I.getWidth() && y < I.getHeight()) {
			I.setRGB(x-1, y, pixel);
			I.setRGB(x, y, pixel);
			I.setRGB(x+1, y, pixel);
			I.setRGB(x, y-1, pixel);
			I.setRGB(x, y+1, pixel);
		}

		return I;
	}
	
	/*
	 * this method takes a binary pixel array and inverts it
	 * 0 -> 255
	 * 255 -> 0
	 */
	public static int[] invertPixels(int[] pixels) {
		int[] inverted = new int[pixels.length];
		
		for(int i = 0; i < pixels.length; i++) {
			inverted[i] = 255-pixels[i];
		}
		
		return inverted;
	}
	
	public static int[] threshold(int[] pixels, int t) {
		int[] threshPixels = new int[pixels.length];
		
		for(int i = 0; i < pixels.length; i++) {
			if(pixels[i] < t) {
				threshPixels[i] = 0;
			}
			else {
				threshPixels[i] = 255;
			}
		}
		return threshPixels;
	}

	/***********************************************************************
	 * all functions below here work on grayscale double pixels, not RGB int
	 ***********************************************************************/

	// rotate pixels by angle
	public static double[] rotate(double[] pixels, int width, int height, double angle) {
		double newPixels[] = null;  

		if ((width*height) == pixels.length) 
		{ 
			angle = angle / 180 * Math.PI;	//program needs angle in radians	

			double midWidth = (double)width/2.0;
			double midHeight = (double)height/2.0;

			double CosAngle = Math.cos(angle);
			double SinAngle = Math.sin(angle);
			double mSinAngle = -Math.sin(angle);
			double transX, transY;
			int newX, newY;
			double newXD, newYD;
			double a, b;
			double newPixel;

			newPixels = new double[width*height]; 
			for (int y = 0; y < height; y++)
			{
				for (int x=0; x<width; x++)
				{
					transX = x-midWidth;
					transY = y-midHeight;
					newXD = (transX*CosAngle + transY * SinAngle) + midWidth;
					newYD = (transX*mSinAngle + transY * CosAngle) + midHeight;
					newX = (int)Math.floor(newXD);
					newY = (int)Math.floor(newYD);

					if ((newX<width-1)&(newX>=0)&(newY<height-1)&newY>=0)
					{
						a = newXD - newX;
						b = newYD - newY;
						newPixel = (1-a)*(1-b)*pixels[newY*width+newX] +
						a*(1-b)*pixels[newY*width+(newX+1)] +
						(1-a)*b*pixels[(newY+1)*width+newX] +
						a*b*pixels[(newY+1)*width+(newX+1)];
						newPixels[y*width+x]=newPixel;
					}
					else
					{
						newPixels[y*width+x]=pixels[0];//probaly white; fix
					}
				}
			}
		} 
		return newPixels; 
	}


	public static double[] adjustContrast (double pixels[], int width, int height, double contrastChange, double brightnessChange)  { 
		double newPixels[] = null; 

		int x, y;
		double pixel;

		if ((width*height) == pixels.length) 
		{ 
			newPixels = new double [width*height]; 
			int newIndex=0; 
			for (y = 0; y < height; y++)
			{
				for (x=0; x<width; x++)
				{
					pixel = pixels[y*width+x];
					pixel = ((pixel-128.0)*contrastChange)+(128.0+brightnessChange);
					newPixels[newIndex++]=pixel;
				}
			}
		} 
		return newPixels; 
	}


	// combine pixels in pix1 with those in pix2 using weights in pix1weight and pix2weight respectively
	public static double[] combinePixels(double[] pix1, double pix1weight, double[] pix2, double pix2weight, int width, int height) {
		double newPixels[] = null; 
		int newIndex;
		double p1, p2;
		int x,y;

		if ((width*height) == pix1.length && (width*height) == pix2.length) {
			newIndex = 0;
			newPixels = new double [width*height]; 
			for (y = 0; y < height; y++) {
				for (x=0; x<width; x++) {
					p1 = pix1[y*width+x];
					p2 = pix2[y*width+x];
					newPixels[newIndex++] = (p1*pix1weight)+(p2*pix2weight);
				}
			}
		}

		return newPixels;
	}


	// linearize pixels will linearize and clip to [0, 255]
	public static double[] linearizePixels(double pixels[], int width, int height, Vector abVector) { 
		double newPixels[] = null; 

		int x, y;
		double pixel, a, b, newLinear, linearVoltage;

		a = Double.parseDouble((String)abVector.elementAt(0));
		b = Double.parseDouble((String)abVector.elementAt(1));

		// double maxLinear = Math.pow(((1.0/a)*255.0),(1/b));
		double maxLinear =Math.pow((1.0/a)*255.0,(1.0/b));

//		System.out.println(" a = " + a + " b = " + b + " maxLinear = " + maxLinear);

		if ((width*height) == pixels.length) 
		{ 
			newPixels = new double [width*height]; 
			int newIndex=0; 
			for (y = 0; y < height; y++)
			{
				for (x=0; x<width; x++)
				{
					pixel = pixels[y*width+x];
					newLinear =(pixel/255.0)*(maxLinear-1.0)+1.0;
					linearVoltage=a*Math.pow(newLinear,b);
					//System.out.println("Old Linear " + pixel + " New Linear " + newLinear + " Linear Voltage " + linearVoltage);
					pixel = linearVoltage;
					newPixels[newIndex++]=pixel;
				}
			}
		} 

		return newPixels; 
	}

}
