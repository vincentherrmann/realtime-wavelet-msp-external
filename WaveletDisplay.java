import java.awt.*;
import java.awt.geom.*;
import java.util.Arrays;

import javax.swing.*;
import com.cycling74.msp.*;

public class WaveletDisplay extends MSPPerformer{
	
	float[][] buffer;
	float scale;
	int bufferSize;
	int nSteps;
	int nTransforms;
	int transformSize;
	int vecSize;			// signal vector size 
	Panel p;
	
	public WaveletDisplay(int steps, int transforms){
		declareInlets(new int[]{SIGNAL, SIGNAL});
		scale = 1.0f;
		declareAttribute("scale");
		
		nSteps = steps;
		if(nSteps<1){nSteps=1;}
		
		transformSize = 1 << nSteps; // 2^nSteps
		
		nTransforms = transforms;
		
		bufferSize = transformSize * nTransforms;
		buffer = new float[2][bufferSize];
		setupDisplay();
	}
	
	public void dspsetup(MSPSignal[] in, MSPSignal[] out){
		vecSize = in[0].n;
		
		bufferSize = transformSize * nTransforms;
		buffer = new float[2][bufferSize];
		post("bufferSize: " + bufferSize);
	}
	
	public void perform(MSPSignal[] in, MSPSignal[] out){
		
		float[] _in1 = in[0].vec;
		float[] _in2 = in[1].vec;
		
		if(transformSize <= vecSize){
			if(bufferSize > vecSize){
				System.arraycopy(buffer[0], vecSize, buffer[0], 0, bufferSize-vecSize);
				System.arraycopy(_in1, 0, buffer[0], bufferSize-vecSize, vecSize);
				
				System.arraycopy(buffer[1], vecSize, buffer[1], 0, bufferSize-vecSize);
				System.arraycopy(_in2, 0, buffer[1], bufferSize-vecSize, vecSize);
			} else{
				System.arraycopy(_in1, vecSize - bufferSize, buffer[0], 0, bufferSize);
				System.arraycopy(_in2, vecSize - bufferSize, buffer[1], 0, bufferSize);
			}
		}
	}
	
	public void bang(){
		//post("repaint?");
		//post(Arrays.toString(buffer[0]));
		p.repaint();
	}
	
	public class Panel extends JPanel{
		
		public void paintComponent(Graphics g) {    	
	        Graphics2D g2 = (Graphics2D) g; 
	        //post("repaint!");
	        //g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	        
	        int w = getWidth();
	        int h = getHeight();
	        
	        //post("width: " + w + ", height: " + h);
	        
	        //int nthTransform;
	        
	        int posX;
	        int posY;
	        int width;
	        int height;
	        int color;
	        
	        for(int i = 0; i < bufferSize; i++){
	        	//nthTransform = i/transformSize;
	        	float level = buffer[1][i];
	        	if(level > 0){
	        		posX = (int) ((i-level)*w)/bufferSize;
	        		posY = (int) ((h*(level-1)) / level);
	        		width = (int) ((2*level*w)/bufferSize);
	        		height = (int) (h/(2*level));
	        		//height = 10;
	        	} else{
	        		posX = (i*w)/bufferSize;
	        		posY =  (int) (h-(1f/transformSize)*h);
	        		width = (int) ((2*(transformSize/2)*w)/bufferSize);
	        		height = (int) (h/(transformSize));
	        	}
	        	color = (int) (Math.abs(buffer[0][i]) * 100);
	        	if(color > 255){
	        		color = 255;
	        	}
	        	g2.setColor(new Color(color, color, color));
	        	g2.fillRect(posX, posY, width+1, height+1);
	        	//post("level: " + level + ", color: " + color + ", posX: " + posX + ", posY: " + posY + ", width: " + width + ", height: " + height);
	        }
	    }
	}
	
	public void setupDisplay(){
    	JFrame f = new JFrame();
        //f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        p = new Panel();
        f.add(p);
        f.setSize(400,400);
        f.setLocation(200,200);
        f.setVisible(true);
    }

}
