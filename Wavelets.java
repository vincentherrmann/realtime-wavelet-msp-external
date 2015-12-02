import java.util.HashMap;

public class Wavelets {
	public static HashMap<String, float[]> writeWavelets(){
		
		HashMap<String, float[]> wavelets = new HashMap<String, float[]>();
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
		
		float[] realK4L2 = 		{-0.0017853301f, 0.013358873f, 0.036090743f, -0.03472219f, 0.041525062f, 0.56035835f, 
								0.77458614f, 0.22752075f, -0.16040927f, -0.06169435f, 0.017099408f, 0.0022852293f};
		float[] imgK4L2 = 		{-3.5706602E-4f, -1.847535E-4f, 0.032591484f, 0.013449902f, -0.058466725f, 0.27464306f,
								0.7795662f, 0.5409738f, -0.04031501f, -0.13320138f, -0.00591213f, 0.011426146f};
		float[] realK8L1 = 		{0.0514121165485221f, 0.300597326079206f, 0.664328833794970f, 0.600105861352414f, 0.00921092271022342f, -0.289580516275856f, 
								-0.0162908610933243f, 0.135603740614867f, -0.00950003071493000f, -0.0502866760081416f, 0.0122668723387378f, 0.0119481746983848f, 
								-0.00528156281796491f, -0.00123601008563526f, 0.00102425884497389f, -5.60257041337798e-05f, -6.37684238030743e-05f, 1.09065163002018e-05f};
		float[] imgK8L1 = 		{0.0171373721828407f, 0.145898767847310f, 0.473407348061759f, 0.706561671558521f, 0.367655611736314f, -0.209867786627167f, 
								-0.225054874431242f, 0.103928677255379f, 0.0977941712914308f, -0.0588603131481208f, -0.0265776141791887f, 0.0251088019645956f, 
								0.00181805310419792f, -0.00629958388124263f, 0.00120527082310358f, 0.000603826669228098f, -0.000278557401810837f, 3.27195489006054e-05f};
		
		
		wavelets.put("DB2", db2);
		wavelets.put("DB4", db4);
		wavelets.put("DB6", db6);
		wavelets.put("DB8", db8);
		wavelets.put("DB10", db10);
		wavelets.put("DB12", db12);
		wavelets.put("DB14", db14);
		wavelets.put("Haar", db2);
		wavelets.put("DiscreteMeyer", discreteMeyer);
		wavelets.put("RealK4L2", realK4L2);
		wavelets.put("ImgK4L2", imgK4L2);
		wavelets.put("RealK8L1", realK8L1);
		wavelets.put("ImgK8L1", imgK8L1);
		
		return wavelets;
	}
}
