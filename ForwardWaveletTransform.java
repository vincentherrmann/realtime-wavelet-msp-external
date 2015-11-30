import com.cycling74.max.*;
import com.cycling74.msp.*;
import java.lang.Math.*;
import java.util.*;

public class FWT6 extends MSPPerformer {
	
	// INSTANCE VARIABLES //
	String wavelet;			// name of the wavelet
	float[] lod; 			// lowpass deconstruction filter
	float[] hid; 			// highpass deconstruction filter
	
	float[][] overlap; 		// signal block overlap for convolution 
	int overlapSize;		// (may be very large, depending on transform level and filter size)
	
	int nSteps;				// number of transform steps
	
	float[][] transform;	// whole transform
	int transformSize;		// Samples of one complete transformation cycle
	int vecSteps;			// total number of signal vector in the whole transform
	
	float[][] delayVector;	// delay
	int nDelays;			// number of delay vectors = number of steps outside the signal vector
	
	int vecSize;			// signal vector size 
	int nF; 				// filter vector size
	int nS; 				// signal size for the transform
	
	float scaleFactor; 		// output scale
	
	int nTransform;			// number of the current signal vector in the whole transform
	int step;				// current step
	
	HashMap<String, float[]> wavelets = new HashMap<String, float[]>();
	
	// CONSTRUCTOR //
	public FWT6(int steps){
		
		declareInlets(new int[]{SIGNAL});
		declareOutlets(new int[]{SIGNAL, SIGNAL});
		wavelet = "DB2";
		declareAttribute("wavelet");
		
		nSteps = steps;
		if(nSteps<1){nSteps=1;}
		
		wavelets = Wavelets.writeWavelets();
		
		// Transform
		transformSize = 1 << nSteps; // 2^nSteps
		transform = new float[2][transformSize]; 
		Arrays.fill(transform[1], -1);
	}
	
	// DSP SETUP //
	public void dspsetup(MSPSignal[] in, MSPSignal[] out){
		
		// Filter
		lod = wavelets.get(wavelet);
		if(lod == null){
			post("No Wavelet '" + wavelet + "', using Haar wavelet instead");
			lod = wavelets.get("DB2");
		}
		hid = calcWaveletCoeffs(lod);
		nF = lod.length;
		scaleFactor = 1.f;
		
		// Transform
		transformSize = 1 << nSteps; // 2^nSteps
		transform = new float[2][transformSize]; 
		Arrays.fill(transform[1], -1);
				
		//Overlap
		overlap = new float[nSteps][0];
		for(int i = 0; i < nSteps; i++){
			overlap[i] = new float[(nF-1) * (1 << i)]; // = (nF-1)*2^i
		}
	
		nTransform = 0;
		
		vecSize = in[0].n;
		nS = vecSize;
		
		vecSteps = transformSize/vecSize;

		if(vecSteps <= 1){
			nDelays = 1;
		} else{
			nDelays = Integer.numberOfTrailingZeros(vecSteps)+1;
		}
		delayVector = new float[nDelays][vecSize];
		for(int i=1; i < nDelays; i++){
			delayVector[i] = new float[vecSize<<(i)];
		}
		
		post("FWT dsp setpup. vecSize: " + vecSize + " nSteps: " + nSteps + " transformSize: "+transformSize + " vecSteps: " +vecSteps + " nDelays: " + nDelays);
		for(int i = 0; i < nF; i++){
			post("lod["+i+"]: " + lod[i]);
			post("hid["+i+"]: " + hid[i]);
		}
		
	}
	
	// PERFORM //
	public void perform(MSPSignal[] in, MSPSignal[] out){

		float[] _in = in[0].vec;
		float[] _out1 = out[0].vec;
		float[] _out2 = out[1].vec;
		
		float[][] currentVector = new float[2][vecSize]; // signal + code
		System.arraycopy(_in, 0, currentVector[0], 0, vecSize);
		
		int regionCode = 1; // start with first step
		int level = 1;
		step = 0;
		
		while((step < nSteps) && (2*level <= vecSize)){
			
			System.arraycopy(currentVector[0], 0, delayVector[0], 0, vecSize);
			
			currentVector = transformStep(currentVector, lod, hid, regionCode);
			
			updateOverlap(delayVector[0]);
			
			regionCode *= 2; // go to next level;
			level = regionCode;
			step++;
		}
		
		if(transformSize > vecSize){
			
			int lastVectorStep = step; // last step that was calculated inside one signal vector
			int remainingSteps = Integer.numberOfTrailingZeros(nTransform);	// trailing zeros of nTransform give the count of further steps
																			// that will be calculated for this (nTransform)th signal vector
			int updateOverlapSteps = Integer.numberOfTrailingZeros(vecSteps - nTransform - 1);

			while(step < nSteps){ // for the remaining steps

				if(currentVector[0].length >= level){
					if(lastVectorStep + remainingSteps > step){ 
						System.arraycopy(currentVector[0], 0, delayVector[step-lastVectorStep+1], 0, level);
					} else{
						System.arraycopy(currentVector[0], 0, delayVector[step-lastVectorStep+1], level, level);
					}
				}
				
				if(lastVectorStep + remainingSteps > step){  
					
					
					float[][] calcVector = new float[2][2*level];
					System.arraycopy(currentVector[0], 0, calcVector[0], 0, currentVector.length);
					System.arraycopy(currentVector[1], 0, calcVector[1], 0, currentVector.length);
					
					calcVector = transformStep(calcVector, lod, hid, regionCode);
					
					currentVector = calcVector;
				}
				
				if(lastVectorStep + updateOverlapSteps > step){ 
					updateOverlap(delayVector[step-lastVectorStep+1]);
				}
				
				regionCode *= 2;
				level = regionCode;
				step++;
			}
			
			int tOffset = vecSize * nTransform;
			System.arraycopy(transform[0], tOffset, _out1, 0, vecSize);
			System.arraycopy(transform[1], tOffset, _out2, 0, vecSize);
	
		} else{
			
		System.arraycopy(currentVector[0], 0, _out1, 0, vecSize);
  		System.arraycopy(currentVector[1], 0, _out2, 0, vecSize);

		}
		
		
		if(nTransform < vecSteps - 1){ // update nTransform
			nTransform++;
		} else {
			nTransform = 0;
			transform = new float[2][transformSize];
			Arrays.fill(transform[1], -1);
		}
		
	}
	
	// TRANSFORM STEP //
	public float[][] transformStep(float[][] signal, float[] lod, float[] hid, int regionCode){
		
		// region code:	approximation = 0 (binary)
								// d1 = 1 		= 1
								// d2 = 10		= 2
								// d3 = 100		= 4
								// ...
		// for the transform step, give the region code of the detail you want to calculate (first step: 1)
		
		// (nTransform+1)th signal Vector of the whole transform
		
		nS = signal[0].length;
		float[][] newSignal = new float[2][nS];
		System.arraycopy(signal[0], 0, newSignal[0], 0, nS);
		System.arraycopy(signal[1], 0, newSignal[1], 0, nS);
		
		float signalValue = 0;
		float tmpA = 0;
		float tmpD = 0;
		
		int tOffset = nTransform*vecSize; // tOffset in transform vector
		
		int level = 1 << step;
		int everyNthSample = 2*level; // take every (2*level)th sample (downsampling!)
		for(int i = 0; i < level; i++){
			
			regionCode = level + i;	
			int offset = regionCode-level; // first sample (packet transform!)

			for(int n = offset; n < nS; n += everyNthSample){ // for every other sample of the current level
					
				// Convolution + Downsampling
				for(int k = 0; k < nF; k++){ // for every sample of the filter
					// find signalValue
					int p = n-k*level;
					if(p < 0){ 								// if it is before the signal vector
						int nLevel = overlap[step].length;
						signalValue = overlap[step][nLevel+p]; 	// take from overlap
					} else if(p > nS){ 						// if it is after the signal vector
						signalValue = 0;						// = 0
					} else {								// if it is on the signal vector
						signalValue = signal[0][p];				// take it from signal vector
					}
								
					// multiply accumulate
					tmpA += signalValue*lod[k]; // accumulate products for approximation
					tmpD += signalValue*hid[k]; // accumulate products for detail
				}
				
				tmpA *= scaleFactor; // scale output
				tmpD *= scaleFactor;
				
				// write values
				int m = n+level; // Detail position

				// write signal vector - not for the final output, but for intern calculation
				if(m < nS){			// if there is place in the signal vector for both values
					newSignal[0][n] = tmpA; 		// write Approx value
					newSignal[0][m] = tmpD; 		// write Detail value
					newSignal[1][m] = regionCode;	// write Detail Code
				} else{						// if there is only place for the Approx value
					newSignal[0][n] = tmpA; 		// write Approx
				}
				
				if(transformSize > vecSize){
					if(transform[1][tOffset+n] < level){
						transform[0][tOffset+n] = tmpA; 		// write Approx value
						transform[1][tOffset+n] = offset;
					} 
					
					if(transform[1][tOffset+m] < level){
						transform[0][tOffset+m] = tmpD;		 	// write Detail value
						transform[1][tOffset+m] = regionCode;	// write Detail code
					} 
				}
				
				tmpA = 0; // reset tmp
				tmpD = 0;
				}
		}
		return newSignal;		
	} 	
	
	// UPDATE OVERLAP //
	public void updateOverlap(float[] signal){

		int nLevel = overlap[step].length; // count of samples of the current level in the overlap
		int length = signal.length;
		
		if(length >= nLevel){ 	// if the signal vector is bigger than the overlap
			System.arraycopy(signal, length - nLevel, overlap[step], 0, nLevel); 		// write from the new vector into the overlap
		} else{					// if the signal vector is smaller than the overlap
			System.arraycopy(overlap[step], length, overlap[step], 0, nLevel-length); 	// shift the values in the overlap to the left
			System.arraycopy(signal, 0, overlap[step], nLevel - length, length);		// write the signal values to the right in the overlap
		}

	}
	
	/*public void writeWavelets(){
		// WAVELETS //
		float[] db2 = 			{0.70710677f, 0.70710677f}; // Haar scale coefficients
		float[] db4 = 			{0.4829629f, 0.8365163f, 0.22414386f, -0.12940952f}; // Daubechies 4 scale coefficients
		float[] db6 =			{0.33267054f, 0.8068915f, 0.4598775f, -0.13501102f, -0.08544128f, 0.035226293f};
		float[] db8 = 			{0.23037781f, 0.71484655f, 0.6308808f, -0.02798377f, -0.18703482f, 0.030841382f, 
					   			 0.03288301f, -0.010597402f};
		float[] db10 = 			{0.1601024f, 0.60382926f, 0.72430855f, 0.13842815f, -0.2422949f, -0.03224487f, 
					    		 0.0775715f, -0.00624149f, -0.012580752f, 0.0033357253f};
		float[] db12 = 			{0.11154074f, 0.4946239f, 0.7511339f, 0.31525034f, -0.2262647f, -0.12976687f, 
								 0.097501606f, 0.027522866f, -0.03158204f, 5.538422E-4f, 0.0047772573f, -0.0010773011f};
		float[] db14 = 			{0.077852055f, 0.39653933f, 0.7291321f, 0.4697823f, -0.143906f, -0.22403619f, 
							 	 0.07130922f, 0.08061261f, -0.038029935f, -0.016574541f, 0.0125509985f, 4.2957798E-4f, 
							 	 -0.0018016407f, 3.537138E-4f};
		
		float[] discreteMeyer = {0.0f, -1.0099999569414229e-12f, 8.519459636796214e-09f, -1.111944952595278e-08f, 
						 		 -1.0798819539621958e-08f, 6.0669757413511352e-08f, -1.0866516536735883e-07f,
								 8.2006806503864813e-08f, 1.1783004497663934e-07f, -5.5063405652522782e-07f,
								 1.1307947017916706e-06f, -1.4895492164971559e-06f, 7.367572885903746e-07f,
								 3.2054419133447798e-06f, -1.6312699734552807e-05f, 6.5543059305751491e-05f,
								 -0.00060115023435160925f, -0.002704672124643725f, 0.0022025341009110021f,
								 0.006045814097323304f, -0.0063877183184971563f, -0.011061496392513451f,
								 0.015270015130934803f, 0.017423434103729693f, -0.032130793990211758f,
								 -0.024348745906078023f, 0.063739024322801596f, 0.030655091960824263f,
								 -0.13284520043622938f, -0.035087555656258346f, 0.44459300275757724f, 
								 0.74458559231880628f, 0.44459300275757724f, -0.035087555656258346f,
								 -0.13284520043622938f, 0.030655091960824263f, 0.063739024322801596f,
								 -0.024348745906078023f, -0.032130793990211758f, 0.017423434103729693f, 
								 0.015270015130934803f, -0.011061496392513451f, -0.0063877183184971563f,
								 0.006045814097323304f, 0.0022025341009110021f, -0.002704672124643725f, 
								 -0.00060115023435160925f, 6.5543059305751491e-05f, -1.6312699734552807e-05f,
								 3.2054419133447798e-06f, 7.367572885903746e-07f, -1.4895492164971559e-06f,
								 1.1307947017916706e-06f, -5.5063405652522782e-07f, 1.1783004497663934e-07f,
								 8.2006806503864813e-08f, -1.0866516536735883e-07f, 6.0669757413511352e-08f,
								 -1.0798819539621958e-08f, -1.111944952595278e-08f, 8.519459636796214e-09f,
								 -1.0099999569414229e-12f};
		
		float[] realN12K4L2 = 	{-0.0017853301f, 0.013358873f, 0.036090743f, -0.03472219f, 0.041525062f, 0.56035835f, 
								0.77458614f, 0.22752075f, -0.16040927f, -0.06169435f, 0.017099408f, 0.0022852293f};
		float[] imgN12K4L2 = 	{-3.5706602E-4f, -1.847535E-4f, 0.032591484f, 0.013449902f, -0.058466725f, 0.27464306f,
								0.7795662f, 0.5409738f, -0.04031501f, -0.13320138f, -0.00591213f, 0.011426146f};
		
		wavelets.put("DB2", db2);
		wavelets.put("DB4", db4);
		wavelets.put("DB6", db6);
		wavelets.put("DB8", db8);
		wavelets.put("DB10", db10);
		wavelets.put("DB12", db12);
		wavelets.put("DB14", db14);
		wavelets.put("Haar", db2);
		wavelets.put("DiscreteMeyer", discreteMeyer);
		wavelets.put("RealN12K4L2", realN12K4L2);
		wavelets.put("ImgN12K4L2", imgN12K4L2);
	}*/
	
	// CALCULATE WAVELET COEFFICIENTS //
	public float[] calcWaveletCoeffs(float[] f){
		int n = f.length;
		float[] s = new float[n];
		for(int i = n-1; i >= 0; i--){ // reverse order and change sign of every other value
			if(i % 2 == 0){
				s[i] = f[n-1 - i];
			} else {
				s[i] = - f[n-1 - i];
			}
		}
		return s;
	}
	
}

