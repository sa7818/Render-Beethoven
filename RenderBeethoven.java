import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.text.View;






// Main class
public class RenderBeethoven extends Frame implements ActionListener {
	RenderCanvas renderer;
	Vertex pos ;
	// Constructor
	public RenderBeethoven(String file) {
		super("Beethoven Bust");
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			int num = Integer.parseInt(line);
			Triangle triangles[] = new Triangle[num];
			for ( int i=0; i<num; i++ )
				triangles[i] = Triangle.parseTriangle(br.readLine());
			renderer = new RenderCanvas(triangles);
		} catch ( Exception e ) {
			System.out.println(e);
		}
		Panel controls = new Panel();
		Button button = new Button("Wireframe");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Depth");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Constant");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Gouraud");
		button.addActionListener(this);
		controls.add(button);
		button = new Button("Phong");
		button.addActionListener(this);
		controls.add(button);
		// add canvas and control panel
		add("Center", renderer);
		add("South", controls);
		addWindowListener(new ExitListener());
		// creat the canvas
		addWindowListener(new ExitListener());
	}
	class ExitListener extends WindowAdapter {
		public void windowClosing(WindowEvent e) {
			System.exit(0);
		}
	}

	// Action listener for buttons
	public void actionPerformed(ActionEvent e) {
		if ( ((Button)e.getSource()).getLabel().equals("Wireframe") )
			renderer.renderWireframe();
	//	else if ( ((Button)e.getSource()).getLabel().equals("Depth") )
	//		renderer.RenderZbuffer();
	//	else if ( ((Button)e.getSource()).getLabel().equals("Constant") )
	//		renderer.renderDemo(Color.black);;
	//	else if ( ((Button)e.getSource()).getLabel().equals("Gouraud") )
	//		renderer.RenderGouraud();
	//	else if ( ((Button)e.getSource()).getLabel().equals("Phong") )
//			renderer.RenderPhong();
		renderer.repaint();
	}
	public static void main(String[] args) {
		RenderBeethoven window = new RenderBeethoven("beethoven.mpt");
		window.setSize(400, 500);
		window.setVisible(true);
	}
}

class RenderCanvas extends Canvas {
	// image for display.
	Image image;
	// the triangles for rendering.
	Triangle triangles[];
	// canvas size
	int width, height;
	// pitch and paw angles for world to camera transform
		int pitch , yaw;// scale ;
		
		int size =3;double scale;
		Vertex pos = new Vertex(0,0,0);//camera position
		Vertex view = new Vertex(0,0,1); //view direction
		Vertex up = new Vertex(0,1,0);
		int focal = 50;
		
		float yrot=0, xrot=0;
		boolean camMove = false;
		boolean check = false;
		Vertex dx,dy,dz;
		 int viewLength =5;
		 int LR_angle=271, UD_angle=5;
		ArrowListener alis = new ArrowListener();
	// initialize the canvas
	public RenderCanvas(Triangle tri[]) {
		//focal = 50; size = 3;
		yaw = 0; pitch = 0;
		triangles = tri;
	  
		DragListener drag = new DragListener();
		addMouseListener(drag);
		addMouseMotionListener(drag);
		addKeyListener(new ArrowListener());
		addComponentListener(new ResizeListener());
	}
	class ResizeListener extends ComponentAdapter {
		public void componentResized(ComponentEvent e) {
			width = getWidth();
			height = getHeight();
			scale = Math.min(width/2, height/2);
			repaint();
		}
		
	}
	// render the wireframe of triangle mesh
	public void renderWireframe() {
		BufferedImage srcImage = new BufferedImage
			(width, height, BufferedImage.TYPE_INT_RGB);
		Triangle tri; //Vertex V = new Vertex(0,0,1);
		Triangle tri2 = null;
		Graphics2D g2D = srcImage.createGraphics();
		
		
		for ( int i=0; i<triangles.length; i++ ) {
			Matrix mat = setupCamera(pos, view, up);
			
		
			    tri = triangles[i].transformCam(width, height,scale,focal,size,yaw,pitch,mat,camMove,pos);
		
				g2D.draw(tri.toPolygon());
		
				
		}//camMove= false;
		image = srcImage;
	}
	// Demo -- generate a solid color image
	public void renderDemo(int color) {
		int pixels[] = new int [width*height];
		for ( int i=0; i<width*height; i++ )
			pixels[i] = color;
		image = createImage(new MemoryImageSource(width, height, pixels, 0, width));
	}
	
	private double[] GetWeights(int x,int y,int[] xCoordinates,int[] yCoordinates)
	{
	
		double[] areas = {0,0,0};
		double[] edgeLen = {0,0,0};
		double[] xyLen = {0,0,0};
		
		edgeLen[0] = Point.distance(xCoordinates[0], yCoordinates[0], xCoordinates[1], yCoordinates[1]);
		edgeLen[1] = Point.distance(xCoordinates[1], yCoordinates[1], xCoordinates[2], yCoordinates[2]);
		edgeLen[2] = Point.distance(xCoordinates[2], yCoordinates[2], xCoordinates[0], yCoordinates[0]);
		xyLen[0] = Point.distance(x, y, xCoordinates[0], yCoordinates[0]);
		xyLen[1] = Point.distance(x, y, xCoordinates[1], yCoordinates[1]);
		xyLen[2] = Point.distance(x, y, xCoordinates[2], yCoordinates[2]);
		
		double sum = 0;
		for(int i = 0;i<3;i++)
		{
			double s = 0.5*(edgeLen[(i+1)%3]+xyLen[(i+1)%3]+xyLen[(i+2)%3]);
			areas[i] = Math.sqrt(s*(s - edgeLen[(i+1)%3])*(s - xyLen[(i+1)%3])*(s - xyLen[(i+2)%3]));
			sum+=areas[i];
		}
		double minValue = 0;
		double[] rtn = {0,0,0};
		rtn[0] = Math.max(minValue, areas[0]/sum);
		rtn[1] = Math.max(minValue, areas[1]/sum);
		rtn[2] = Math.max(minValue, areas[2]/sum);
		
		return rtn;
	}


	private Vertex GetSurfaceNormal(Triangle t)
	{
		Vertex a=new Vertex(t.coords[1].x - t.coords[0].x,
				t.coords[1].y - t.coords[0].y , t.coords[1].z - t.coords[0].z);
		Vertex b=new Vertex(t.coords[2].x - t.coords[1].x,
				t.coords[2].y - t.coords[1].y , t.coords[2].z - t.coords[1].z);
		
		Vertex v = new Vertex(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
		v.Normalize();
		return v;
	}
	
   public double dot(Vertex v1, Vertex v2)
   {
	   double temp;
	   return temp = v1.x* v2.x+v1.y*v2.y+v1.z*v2.z;
   }
	public void RenderPhong()
	{
		
		int pixels[] = new int [width*height];
		Vertex V = new Vertex(0, 0, 1);
		Vertex temp = new Vertex(0,0,0);
		Triangle tri;
		for(int i = 0;i<(width*height);i++)
		{
			pixels[i] = Color.BLACK.getRGB();
		}
		
		int lastValue = 0;
		
		for (int i = 0;i<triangles.length;i++)
		{
			Matrix mat = setupCamera(pos, view, up);
	
			    tri = triangles[i].transform(width, height,scale,focal,size,yaw,pitch,mat,camMove,pos);
			
			Vertex NN = GetSurfaceNormal(tri);
			
			
			temp.x = tri.coords[0].x - V.x;
			temp.y = tri.coords[0].y - V.y;
			temp.z = tri.coords[0].z - V.z;
			if(dot(NN, temp)<0)
		//	if ( tri.normals[0].z<0 && tri.normals[1].z<0 && tri.normals[2].z<0 )
	{
				
				
				/// do transformation on normals here//////////////////////////////////////////////////////////////////////////
				Vertex n0 = tri.normals[0];  
				Vertex n1 = tri.normals[1];
				Vertex n2 = tri.normals[2];
				
				
				Polygon poly = tri.toPolygon();
				Rectangle rect = poly.getBounds();
				
				for(int h=rect.y;h<(rect.y+rect.height);h++)
					for(int w=rect.x;w<(rect.x+rect.width);w++)
					{
						if(poly.contains(w, h))
						{
						
							double[] weights = GetWeights(w, h, poly.xpoints, poly.ypoints);
							
							Vertex N = new Vertex(0, 0, 0);
							N.x = n0.x*weights[0]+n1.x*weights[1]+n2.x*weights[2];
							N.y = n0.y*weights[0]+n1.y*weights[1]+n2.y*weights[2];
							N.z = n0.z*weights[0]+n1.z*weights[1]+n2.z*weights[2];
						
							Vertex L = new Vertex(-0.5, -0.3, -1);
							L.Normalize();
							
							Vertex R = new Vertex(0,0,0);
							double NL = N.x*L.x+N.y*L.y+N.z*L.z;					
							R.x = 2*NL*N.x-L.x;
							R.y = 2*NL*N.y-L.y;
							R.z = 2*NL*N.z-L.z;
							R.Normalize();
							
							
						//	V = new Vertex(0, 0, 1);
						//	V.Normalize();
							
							double I = 0.5*(N.x*L.x+N.y*L.y+N.z*L.z)+
									0.5*Math.pow(V.x*R.x+V.y*R.y+V.z*R.z, 10);
							
							if(I<0)I*=-1;
							
							pixels[h*width+w] = (255<<24)|((int)(224*I)<<16)|((int)(144*I)<<8)|(int)(32*I);
							
							if(pixels[h*width+w] == Color.black.getRGB())
								pixels[h*width+w] = lastValue;
							else
								lastValue = pixels[h*width+w];
						
					}
					}
			}
		}

		image = createImage(new MemoryImageSource(width, height, pixels, 0, width));
	}
	
	public Matrix setupCamera(Vertex pos, Vertex view, Vertex up)
	{
		Vertex Pos = new Vertex(-pos.x, -pos.y, -pos.z);
		
		//Vertex Pos = new Vertex(0,0,0);
		Matrix mat;

		dz = view.normalize();
		dx = dz.cross(up);
		dx = dx.normalize();
		dy = dx.cross(dz);
	    mat = new Matrix(dx,dy,dz);
	    mat.translate(Pos);
		//}
		 
	
		return mat;
	}
	
	// End of demo
	public void paint(Graphics g) {
		// draw the rendering result.


		renderWireframe();
	//	RenderPhong();
		g.drawImage(image, 0, 0, this);
	}
	 
		// Action listener for mouse
		class DragListener extends MouseAdapter implements MouseMotionListener {
			int lastX, lastY;
			public void mousePressed(MouseEvent e) {
				lastX = e.getX();
				lastY = e.getY();
			}
			public void mouseMoved(MouseEvent e) {}
			// update pitch and yaw angles when the mouse is dragged.
			public void mouseDragged(MouseEvent e) {
				yaw += e.getX() - lastX;
				pitch += e.getY() - lastY;
				lastX = e.getX();
				lastY = e.getY();
				System.out.println("yaw "+yaw+" pitch"+pitch);
				repaint();
			}
		}
		// Action listener for keyboard
		class ArrowListener extends KeyAdapter {
			//public  boolean camMove = false;
			float xrotrad, yrotrad;
			Transformation trans = new Transformation();
			Triangle tri = new Triangle();
			double[][] r = new double[4][4]; double [] temp = new double[4];
			public void keyPressed(KeyEvent e) {
				if ( e.getKeyCode() == e.VK_DOWN )//&& focal>5 )
				{
					
				   UD_angle+=1;
				  //  viewLength +=2;
					pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
					pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
					pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
					
					//pos.x -= view.x;
					//pos.z -= view.z;
				}
				else if ( e.getKeyCode() == e.VK_UP)// && focal<50 )
				{
					
					 UD_angle-=1;
				//	viewLength -=2;
					pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
					pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
					pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
		
				//pos.x += view.x;
			//		pos.z += view.z;
			
				}
				 if ( e.getKeyCode() == e.VK_LEFT )//&& size>1 )
				{
					 LR_angle +=1;
					pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
					pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
					pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
					 
					camMove = true;
					
				//	pos.x -= up.cross(view).x;
				//	pos.z -= up.cross(view).z;
				}
					
					//pos.y -=1;//size --;
				else if ( e.getKeyCode() == e.VK_RIGHT)// && size<20 )
				{
					 LR_angle -=1;
					 pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
						
					
				//	camMove = true;
				//	pos.x += up.cross(view).x;
				//	pos.z += up.cross(view).z;
					
				}
				else if( e.getKeyCode() == e.VK_A)//left
				{
					    LR_angle += 1;
					    
					    pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
					      
				}
				else if( e.getKeyCode() == e.VK_D)//right
				{
					    LR_angle-=1;
					    pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
					
				 //  if (yrot >360) yrot -= 360;
				}
				else if( e.getKeyCode() == e.VK_W)//zoom in 
				{
					    viewLength +=2;
					    pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
						System.out.println(up.y);
				}
				else if( e.getKeyCode() == e.VK_S)//zoom out
				{
					    viewLength -=2;
					    pos.z = Math.sin(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.x = Math.cos(Math.toRadians(LR_angle))*Math.cos(Math.toRadians(UD_angle))*viewLength;
						pos.y = Math.sin(Math.toRadians(UD_angle))*viewLength;
		
				}
				 System.out.println("LRAngle "+ LR_angle + " UD"+ UD_angle+ " viewLenght" + viewLength);
				 repaint();
			}
			public boolean returnCamMove(){return camMove;}
		}
	   public  Vertex ArbitraryRotation(double yrot,Vertex axis,Vertex ver)
		{
		   double [] result = new double[4];double [] temp = new double[4];
		   double[][] translated1 = new double[4][4];  double[][] rotated = new double[4][4];  double[][] translated2 = new double[4][4];
		   double[][] last = new double[4][4]; double[][] temp1 = new double[4][4];
		   Vertex a = new Vertex(axis.x-pos.x, axis.y-pos.y,axis.z-pos.z);
		   Transformation trans = new Transformation();
			//a is up direction vertex
			double c = Math.cos(yrot); 
			double s = Math.sin(yrot);
			double t = 1-c;
			double[][] mat = new double[4][4];
			mat[0][0] = a.x*a.x*t+c ; 	  mat[0][1] = a.x*a.y*t-a.z*s; mat[0][2] = a.x*a.z*t+a.y*s; mat[0][3] =0;
			mat[1][0] = a.y*a.x*t+a.z*s;  mat[1][1] = a.y*a.y*t+c;     mat[1][2] = a.y*a.z*t-a.x*s;  mat[1][3] =0;
			mat[2][0] = a.z*a.x*t-a.y*s;  mat[2][1] = a.z*a.y*t+a.x*s; mat[2][2] = a.z*a.z*t+c; mat[2][3] =0;
			mat[3][0] = 0;                mat[3][1] = 0;               mat[3][2] = 0;    mat[3][3] =1;
		
		Triangle tri = new Triangle();
		temp[0] =ver.x;temp[1] =ver.y;temp[2] =ver.z;temp[3] =1;
		translated1 = trans.translationMatrix(new Vertex(-(a.x), -(a.y),-(a.z)));
		rotated = tri.multiplyTwoMatric(mat, translated1);
		translated2 = trans.translationMatrix(new Vertex(a.x, a.y,a.z));
		last = tri.multiplyTwoMatric( translated2,rotated);


		return new Vertex(result[0],result[1],result[2]);
		}
}

// Vertex definition: coordinates for a 3D point
class Vertex {
	double x, y, z;
	// constructors
	public Vertex(double a, double b, double c) {
		x = a; y = b; z = c;
	}
	
	public void Normalize()
	{
		double m = Math.sqrt(x*x+y*y+z*z);
		x/=m;y/=m;z/=m;
	}
	
	public static double Distance(Vertex v1,Vertex v2)
	{
		return Math.sqrt(Math.pow(v1.x-v2.x, 2)+Math.pow(v1.y-v2.y, 2)+Math.pow(v1.z-v2.z, 2));
	}
	
	public static Vertex parseVertex(String tokens[], int start) {
		double x = Double.parseDouble(tokens[start]);
		double y = Double.parseDouble(tokens[start+1]);
		double z = Double.parseDouble(tokens[start+2]);
		return new Vertex(x, y, z);
	}
	public Vertex normalize() {
		double l = length();
		return new Vertex(x/l, y/l, z/l);
	}
	public double dot(Vertex t) {
		return x*t.x + y*t.y + z*t.z;
	}
	public Vertex cross(Vertex t) {
		return new Vertex(y*t.z-z*t.y, z*t.x-x*t.z, x*t.y-y*t.x);
	}
	public double length() {
		return Math.sqrt(x*x+y*y+z*z);
	}
}

// Triangle definition
class Triangle {
	// positions and normals of the 3 coords
	Vertex coords[] = new Vertex[3], normals[] = new Vertex[3];
	// parse a triangle from a string
	public static Triangle parseTriangle(String input) {
		String tokens[] = input.split("\\s");
		Triangle tri = new Triangle();
		for ( int n=0 ; n<3 ; n++ ) {
			tri.coords[n] = Vertex.parseVertex(tokens, n*3);
			tri.normals[n] = Vertex.parseVertex(tokens, 9+n*3);
		}
		return tri;
	}
	
	

	public Triangle transformCam(int width, int height, double scale,
			int focal, int size, int yaw, int pitch, Matrix mat, boolean camMove, Vertex pos) {
		Triangle result = new Triangle();
		Vertex[] temp = new Vertex[3];
		double[] a1 = new double[4];
		
		for ( int n=0 ; n<3 ; n++ ) 
		{
			double offset = 10; //scale
			scale = Math.min(width/10.,height/12.5);
			temp[n] = screen2Raster(camera2Screen(world2CameraCam(coords[n], yaw, pitch,mat,pos),scale,focal,size),scale,focal,size,width,height);
			result.coords[n] = new Vertex(temp[n].x, temp[n].y,temp[n].z);
			a1 = world2CameraCam(normals[n], yaw, pitch,mat,pos); 
			result.normals[n] = new Vertex(a1[0], -a1[1],-a1[2]);
		}
		return result;
	}

public double[] world2CameraCam(Vertex from, int yaw, int pitch,Matrix mat,Vertex pos) {
	
	Triangle tri = new Triangle();
	Transformation trans = new Transformation();
	double[][] translated = new double[4][4];	double[][] rotated = new double[4][4];
	double[][] rotated2 = new double[4][4];
	double[][] temp1 = new double[4][4]; double[][] temp2 = new double[4][4]; double[][] translated2 = new double[4][4];
	double[] tr = new double[4];	
	
	
        double[][] temp = new double[4][4];
   //     temp = multiplyTwoMatric(mat.mat, )
		double[] vec = new double[4];double[] tem = new double[4];
		vec[0] = from.x;vec[1] = from.y;vec[2] = from.z;vec[3] = 1;
		tem = multiplywithVertex(vec, mat.mat); //the whole transformation matrix is now multiplied by the point
		
		translated = trans.translationMatrix(new  Vertex(-pos.x,-pos.y,-pos.z));
		translated2 = trans.translationMatrix(new Vertex(pos.x, pos.y ,pos.z));
	
		
		rotated = trans.rotateMatrix(Math.toRadians(yaw), 1); // computing the matrix of roation by yaw
		rotated2 = trans.rotateMatrix(Math.toRadians(pitch), 2); // computing the matrix of rotation by pitch

		temp1 = multiplyTwoMatric(rotated, translated2);
		temp2 = multiplyTwoMatric(temp1, rotated2);
		temp1 = multiplyTwoMatric(temp2, translated2);

        
		tem = multiplywithVertex(tem, temp1); //the whole transformation matrix is now multiplied by the point
		vec[0] = tem[0]; vec[1] = tem[1];vec[2] = tem[2];vec[3] = tem[3];
		
			
		return (tem);
			
		
	
	}

	public Vertex WeightedPoint(double[] w)
	{
		Vertex v = new Vertex(0,0,0);
		
		for(int i=0;i<3;i++)
		{
			v.x+=(coords[i].x*w[i]);
			v.y+=(coords[i].y*w[i]);
			v.z+=(coords[i].z*w[i]);
		}
		
		return v;
	}
	
	
	// transform the triangle to the screen coordinate
	public Triangle transform(int width, int height,double scale,int focal,int size, int yaw,int pitch,Matrix mat, boolean camMove, Vertex pos) 
	{
		Triangle result = new Triangle();
		Vertex[] temp = new Vertex[3];
		Vertex[] temp2 = new Vertex[3];
		double[] a1 = new double[4];
		
		for ( int n=0 ; n<3 ; n++ ) 
		{
			
			double offset = 10; //scale
			scale = Math.min(width/10.,height/12.5);
	
			temp[n] = screen2Raster(camera2Screen(world2Camera(coords[n], yaw, pitch,mat,camMove,pos),scale,focal,size),scale,focal,size,width,height);
			result.coords[n] = new Vertex(temp[n].x, temp[n].y,temp[n].z);
			//result.normals[n] = new Vertex(normals[n].x, -normals[n].y, -normals[n].z);
			a1 = world2Camera(normals[n], yaw, pitch,mat,camMove,pos); 
			//temp1[n] = screen2Raster(camera2Screen(world2Camera(normals[n], yaw, pitch),scale,focal,size),scale,focal,size,width,height);
			//result.normals[n] = new Vertex(temp1[n].x, -temp1[n].y,-temp1[n].z);
			result.normals[n] = new Vertex(a1[0], -a1[1],-a1[2]);
		}
		return result;
	}

	// convert 3D triangle to a 2D polygon
	public Polygon toPolygon() {
		Polygon poly = new Polygon();
		poly.addPoint((int)coords[0].x, (int)coords[0].y);
		poly.addPoint((int)coords[1].x, (int)coords[1].y);
		poly.addPoint((int)coords[2].x, (int)coords[2].y);
		return poly;
	}
	public double[] change2double(Vertex coords)
	{   double[] temp = new double[3];
		for(int i=0; i<3; i++)
		{
			temp[0] = coords.x;
			temp[1] = coords.y;
			temp[2] = coords.z;
		}
		return temp;
	}
public Vertex camera2Screen(double[] from,double scale, int focal, int size) {
		
		double[] temp = new double[4]; Vertex vecOut = new Vertex(0,0,0);
		
		temp[0]= from[0]; temp[1] =from[1]; temp[2] = from[2];	
		vecOut = makeHomog(temp, focal, size); //call the functin to multiply point with the projection matrix
		return new Vertex(vecOut.x, vecOut.y, vecOut.z);
		
	}
	// Transform between screen and raster coordinates
	public Vertex screen2Raster(Vertex from, double scale, int focal, int size, int width, int height) {
		return new Vertex(from.x*scale+width/2, -from.y*scale+height/2, 0);
	}
	public double[] world2Camera(Vertex from, int yaw, int pitch, Matrix mat,boolean camMove, Vertex pos) {
	
		
		Triangle tri = new Triangle();
		Transformation trans = new Transformation();
		double[][] translated = new double[4][4];	double[][] rotated = new double[4][4];
		double[][] rotated2 = new double[4][4];
		double[][] temp1 = new double[4][4]; double[][] temp2 = new double[4][4]; double[][] translated2 = new double[4][4];
		double[] tr = new double[4];
		// here I first compute the translation to origion matrix 
		//***********move to where keyboard last moved object
		tr[0] = from.x; tr[1] = from.y; tr[2] = from.z; 

		tr = tri.multiplywithVertex(tr, mat.mat);
	

		rotated = trans.rotateMatrix(Math.toRadians(yaw), 1); // computing the matrix of roation by yaw
		rotated2 = trans.rotateMatrix(Math.toRadians(pitch), 2); // computing the matrix of rotation by pitch

       temp1 = multiplyTwoMatric( rotated,rotated2);
        
		double[] vec = new double[4];double[] tem = new double[4];
	//	vec[0] = from.x;vec[1] = from.y;vec[2] = from.z;vec[3] = 1;
		tem = multiplywithVertex(tr, temp1); //the whole transformation matrix is now multiplied by the point
		
		vec[0] = tem[0]; vec[1] = tem[1];vec[2] = tem[2];vec[3] = tem[3];
		
		return (tem);
	
	}
	
	public double[] multiplywithVertex(double[] vec, double[][] matrix )
	{
		double[] temp = new double[4];
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				temp[i]+= matrix[i][j]*vec[j];
			}
		}
		return temp;
	}
public double[][] multiplyTwoMatric(double a[][], double b[][]){
		
		double[][] temp = new double[4][4];
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<4; k++)
				{
				temp[i][j] += a[i][k]*b[k][j];
				}
			}
		}
		return temp;}
	public Vertex makeHomog(double[] ver, int focal, int size)
	{
		//compute the homogenous perspective projection matrix.
		double[] vertex = new double[4];
		
  		vertex[0] = (ver[0]*focal)/(size*(focal+ver[2]));
		vertex[1] = (ver[1]*focal)/(size*(focal+ver[2]));
		vertex[2] = 0;
		vertex[3] = 1;
	
		return new Vertex(vertex[0],vertex[1],vertex[2]);
	}

}
class Matrix{
	Transformation tra = new Transformation();
	Triangle tri = new Triangle();
	double[][] mat;
	public Matrix(Vertex dx, Vertex dy, Vertex dz){
		
	mat = new double[4][4];
	mat[0][0] = dx.x ; mat[0][1] = dy.x; mat[0][2] = dz.x; mat[0][3] =0;
	mat[1][0] = dx.y;  mat[1][1] = dy.y; mat[1][2] =dz.y;  mat[1][3] =0;
	mat[2][0] = dx.z;  mat[2][1] = dx.z; mat[2][2] = dz.z; mat[2][3] =0;
	mat[3][0] = 0;     mat[3][1] = 0;    mat[3][2] = 0;    mat[3][3] =1;
	}
	public void translate(Vertex translate)
	{
		double[][] translation = new double[4][4];
		double[][] result = new double[4][4];
		translation = tra.translationMatrix(translate);
		result = tri.multiplyTwoMatric(mat, translation);
		mat = result;
	//	mat[0][3] = translate.x; mat[1][3] = translate.y;mat[2][3] = translate.z;
		//return result = tri.multiplyTwoMatric(mat,tra.translationMatrix(translate));
	}
}
class Transformation{
	public double[][] rotateMatrix(double pitch, int type)// if yaw 1 
	{
		// compute the rotation matrices along x or z
		double[][] rotation = new double[4][4];
		if(type ==1) //yaw 
		{
		rotation[0][0] = Math.cos(pitch);  rotation[0][1] = - Math.sin(pitch);rotation[0][2] = 0;rotation[0][3] =0;
		rotation[1][0] = Math.sin(pitch);  rotation[1][1] =  Math.cos(pitch); rotation[1][2] = 0;rotation[1][3] =0;
		rotation[2][0] = 0;                rotation[2][1] = 0;                rotation[2][2] = 1;rotation[2][3] =0;
		rotation[3][0] = 0;                rotation[3][1] = 0;                rotation[3][2] = 0;rotation[3][3] =1;
		}
		if(type ==2) //pitch
		{
			rotation[0][0] = 1;rotation[0][1] = 0;                   rotation[0][2] = 0;                 rotation[0][3] =0;
			rotation[1][0] = 0;rotation[1][1] =  Math.cos(pitch);    rotation[1][2] =- Math.sin(pitch);  rotation[1][3] =0;
			rotation[2][0] = 0;rotation[2][1] = Math.sin(pitch);     rotation[2][2] =  Math.cos(pitch);  rotation[2][3] =0;
			rotation[3][0] = 0;rotation[3][1] = 0;                   rotation[3][2] = 0;                 rotation[3][3] =1;
		}
		else //along Y
		{
			rotation[0][0] = Math.cos(pitch);  rotation[0][1] = 0;  rotation[0][2] = Math.sin(pitch);rotation[0][3] =0;
			rotation[1][0] = 0;  			   rotation[1][1] = 1;  rotation[1][2] = 0;			     rotation[1][3] =0;
			rotation[2][0] = -Math.sin(pitch); rotation[2][1] = 0 ; rotation[2][2] = Math.cos(pitch);rotation[2][3] =0;
			rotation[3][0] = 0;                rotation[3][1] = 0;  rotation[3][2] = 0;			     rotation[3][3] =1;
		}
		return rotation;
	}
	//to compute translation matrices 
	public double[][] translationMatrix(Vertex vec)
	{
		double[][] translation = new double[4][4];
		translation[0][0] = 1;translation[0][1] = 0;translation[0][2] = 0;translation[0][3] = vec.x;
		translation[1][0] = 0;translation[1][1] = 1;translation[1][2] = 0;translation[1][3] = vec.y;
		translation[2][0] = 0;translation[2][1] = 0;translation[2][2] = 1;translation[2][3] = vec.z;
		translation[3][0] = 0;translation[3][1] = 0;translation[3][2] = 0;translation[3][3] = 1;

		return translation;
		
	}
}