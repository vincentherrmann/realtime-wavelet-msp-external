import com.cycling74.msp.*;

import java.util.*;

public class IFWT6 extends MSPPerformer {
	
	// INSTANCE VARIABLES //
	String wavelet;			// name of the wavelet
	float[] lor; 			// lowpass reconstruction filter
	float[] hir; 			// highpass reconstruction filter
	
	float[][] overlap; 		// signal block overlap for convolution 
	int overlapSize;		// (may be very large, depending on transform level and filter size)
	
	int nSteps;				// number of transform steps
	int step;				// current step
	
	float[][] transform;		// whole transform
	float[] inputTransform;
	int transformSize;		// Samples of one complete transformation cycle
	int nTransform;			// number of the current signal vector in the whole transform
	int vecSteps;			// total number of signal vector in the whole transform
	
	int vecSize;			// signal vector size 
	int nF; 				// filter vector size
	int nS; 				// signal size for the transform
	
	int nDelay;				// delay in samples (level 1)
	
	float scaleFactor; 		// output scale
	
	HashMap<String, float[]> wavelets = new HashMap<String, float[]>();
	
	// CONSTRUCTOR //
	public IFWT6(int steps){
		nSteps = steps;
		declareInlets(new int[]{SIGNAL, SIGNAL});
		declareOutlets(new int[]{SIGNAL});
		
		wavelet = "DB2";
		declareAttribute("wavelet");
		
		wavelets = Wavelets.writeWavelets();
	}
	
	// DSP SETUP //
	public void dspsetup(MSPSignal[] in, MSPSignal[] out){
		
		// Filter
		lor = wavelets.get(wavelet);
		if(lor == null){
			post("No Wavelet '" + wavelet + "', using Haar wavelet instead");
			lor = wavelets.get("DB2");
		}
		hir = reverseOrder(calcWaveletCoeffs(lor));
		lor = reverseOrder(lor);
		nF = lor.length;
		scaleFactor = 1.f;
		nDelay = nF-1;
	
		nTransform = 0;

		vecSize = in[0].n;
		nS = vecSize;
		
		if(nSteps == 0){
			nSteps = 1;
		} 
		
		transformSize = 1 << nSteps; // 2^nSteps
		transform = new float[2][transformSize]; 
		inputTransform = new float[transformSize];
		Arrays.fill(transform[0], 0);
		vecSteps = transformSize/vecSize;
		
		overlapSize = (nF-1) * (transformSize) + (transformSize - 2) * nDelay; // space for filter + delay
		overlap = new float[nSteps][0];
		for(int i = 0; i < nSteps; i++){
			if(i < nSteps-1){
				overlap[i] = new float[(nF-1) * (1<<i) + 2*(1<<i)*nDelay]; // = (nF-1)*level + 2*level*nDelay
			} else{
				overlap[i] = new float[(nF-1) * (1<<i)]; // = (nF-1)*2^i
			}		
		}
		
		post("iFWT dsp setpup. vecSize: " + vecSize + ", nSteps: " + nSteps + ", transformSize: "+transformSize);
		for(int i = 0; i < nF; i++){
			post("lor["+i+"]: " + lor[i]);
			post("hir["+i+"]: " + hir[i]);
		}
		
	}
	
	// PERFORM //
	public void perform(MSPSignal[] in, MSPSignal[] out){

		float[] _in1 = in[0].vec;
		float[] _in2 = in[1].vec;
		float[] _out = out[0].vec;
		
		float[] currentVector;// = new float[2][vecSize]; // signal + code
		float[] delayVector;
		currentVector = (float[]) _in1.clone(); // [0][] = signal / [1][] = code

		int regionCode = 1 << (nSteps-1); // start with first step
		step = nSteps-1;
		
		if(transformSize > vecSize){
			
			int preSteps = Integer.numberOfTrailingZeros(nTransform); // number of steps, that can't be calculated inside the signal vector
			int updateOverlapSteps = Integer.numberOfTrailingZeros(vecSteps - nTransform - 1);
			int stepsInVector = Integer.numberOfTrailingZeros(vecSize);
			
			while(2*regionCode > vecSize){
				
				int level = regionCode;
				
				//calculate transform step
				if(preSteps + stepsInVector >= step+1){
					currentVector = new float[2*level];
					System.arraycopy(transform[0], nTransform*vecSize, currentVector, 0, 2*level);
					
					delayVector = new float[2*level];
					System.arraycopy(currentVector, 0, delayVector, 0, 2*level);
					
					currentVector = inverseTransformStep(currentVector, lor, hir, regionCode);
					
					updateOverlap(delayVector);
					System.arraycopy(currentVector, 0, transform[0], nTransform*vecSize, 2*level);
				}
				
				regionCode /= 2;
				step--;
			}
			
			int level = regionCode;
			int totalStepDelay = (transformSize-2*level)*nDelay;
			
			currentVector = new float[vecSize];
			int nLevel = overlap[step].length;

			System.arraycopy(transform[0], nTransform*vecSize, currentVector, 0, vecSize);
		}
		

		delayVector = new float[vecSize];
		
		while(step >= 0){
			

			System.arraycopy(currentVector, 0, delayVector, 0, vecSize);
					
			currentVector = inverseTransformStep(currentVector, lor, hir, regionCode);
					
			updateOverlap(delayVector);
			
			regionCode /= 2; // go to next level;
			step--;
		}
		
		
		if(transformSize > vecSize){
			// write inputTransform
			System.arraycopy(_in1, 0, inputTransform, nTransform*vecSize, vecSize); 
		}
		
		System.arraycopy(currentVector, 0, _out, 0, vecSize);
			
		
		if(nTransform < vecSteps - 1){ // update nTransform
			nTransform++;
		} else {
			Arrays.fill(transform[1], 0);
			System.arraycopy(inputTransform, 0, transform[0], 0, transformSize);
			nTransform = 0;
		}
	}	
	
	public float[] inverseTransformStep(float[] signal, float[] lor, float[] hir, int regionCode){
		
		float signalValue = 0;
		float tmp1 = 0;
		float tmp2 = 0;
				
		int tOffset = nTransform*vecSize; // tOffset in transform vector
		
		int level = 1 << step;
		int everyNthSample = level << 1; // take every (2*level)th sample (downsampling!)
		
		nS = signal.length;
		int nLevel = overlap[step].length;
		int delay = 0;
		if(2*level < transformSize){
			delay = 2*level * nDelay;
		}
		
		// write newSignal (delayed)
		float[] newSignal = new float[nS];
		if(delay < nS){
			// source, srcPos, destination, destPos, length
			System.arraycopy(signal, 0, newSignal, delay, nS - delay);
			System.arraycopy(overlap[step], nLevel - delay, newSignal, 0, delay);
		} else{
			System.arraycopy(overlap[step], nLevel - delay, newSignal, 0, nS); // fill the whole array with values from the overlap
		}

		for(int i = 0; i < level; i++){
			
			regionCode = level + i;	
			int offset = regionCode-level; // first sample (packet transform!)
			
			// Approximation (the values calculated before)
			for(int n = offset; n < nS; n +=  everyNthSample){

				for(int k = 0; k < nF; k += 2){
					
					// find signalValue
					int p = n - k*level;
					if(p < 0){ 								// if it is before the signal vector
						signalValue = overlap[step][nLevel+p]; 	// take from overlap
					} else if(p > nS){ 						// if it is after the signal vector
						signalValue = 0;						// = 0
					} else {								// if it is on the signal vector
						signalValue = signal[p];			// take it from signal vector
					}
					
					// multiply accumulate
					tmp1 += signalValue*lor[k]; // accumulate products for first value
					tmp2 += signalValue*lor[k+1]; // accumulate products for second value
				}
				
				newSignal[n] = (tmp1*scaleFactor);
				newSignal[n+level] = (tmp2*scaleFactor);

				if(level*2 > vecSize){
					transform[0][n+tOffset] = tmp1*scaleFactor;
					transform[0][n+level+tOffset] = tmp2*scaleFactor;
				}

				tmp1 = 0;
				tmp2 = 0;
			}

			// Detail
			for(int n = offset+level; n < nS; n +=  everyNthSample){

				// first value
				for(int k = 0; k < nF; k += 2){
					
					// find signalValue
					int p = n - delay - k*level;
					if(p < 0){ 								// if it is before the signal vector
						signalValue = overlap[step][nLevel+p]; 	// take from overlap
					} else if(p > nS){ 						// if it is after the signal vector
						signalValue = 0;						// = 0
					} else {								// if it is on the signal vector
						signalValue = signal[p];				// take it from signal vector
					}
					
					// multiply accumulate
					tmp1 += signalValue*hir[k]; // accumulate products for first value
					tmp2 += signalValue*hir[k+1]; // accumulate products for second value
				}
				
				newSignal[n-level] += (tmp1*scaleFactor);
				newSignal[n] += (tmp2*scaleFactor);

				if(level*2 > vecSize){
					transform[0][n-level+tOffset] += tmp1*scaleFactor;
					transform[0][n+tOffset] += tmp2*scaleFactor;
				}
				
				tmp1 = 0;
				tmp2 = 0;
				
			}
		}
	
		return newSignal;
	}

	
	public void updateOverlap(float[] signal){
		
		int length = signal.length;
		
		int nLevel = overlap[step].length; // count of samples of the current level in the overlap
		
		if(length >= nLevel){
			// source, srcPos, destination, destPos, length
			System.arraycopy(signal, length - nLevel, overlap[step], 0, nLevel); // write the new vector into the overlap
		} else{
			// source, srcPos, destination, destPos, length
			System.arraycopy(overlap[step], length, overlap[step], 0, nLevel - length);
			System.arraycopy(signal, 0, overlap[step], nLevel - length, length);
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
	
	
	public float[] reverseOrder(float[] f){
		int n = f.length;
		float[] s = new float[n];
		for(int i = n-1; i >= 0; i--){ // reverse order
				s[i] = f[n-1 - i];
		}
		return s;
	}
	
	public static float[] calcWaveletCoeffs(float[] f){
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

