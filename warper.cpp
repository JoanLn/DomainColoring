/*
 * Program sample to warp an image using matrix-based warps
 * 
 * Command line parameters are as follows:
 *
 * warper infile.png [outfile.png] 
 *
 * Author: Joanna Lin, 11/28/22, modifed on 12/08/22

 */

#include <OpenImageIO/imageio.h>
#include <GL/glut.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <complex>
#include <stdlib.h>

using namespace std;
using std::string;
OIIO_NAMESPACE_USING

struct Pixel{ // defines a pixel structure
	unsigned char r,g,b,a;
}; 

float MAX(float a, float b){return  b < a ? a : b;}
float MIN(float a, float b){return  b > a ? a : b;}

const double PI = 3.1415926535897932384626433832795;
const double E  = 2.7182818284590452353602874713527;

//
// Global variables and constants
//
const int DEFAULTWIDTH = 600;	// default window dimensions if no image
const int DEFAULTHEIGHT = 600;

int WinWidth, WinHeight;	// window width and height
int ImWidth = DEFAULTWIDTH, ImHeight = DEFAULTHEIGHT;		// image width and height
int ImChannels;           // number of channels per image pixel
int OutWidth = DEFAULTWIDTH, OutHeight = DEFAULTHEIGHT;		// image width and height

int VpWidth, VpHeight;		// viewport width and height
int Xoffset, Yoffset;     // viewport offset from lower left corner of window

Pixel **pixmap = NULL;     // the image pixmap used for OpenGL display
Pixel **pixmapOut = NULL;  // the image pixmap used for OpenGL display
int pixformat; 	   // the pixel format used to correctly  draw the image

string outputFileName;  // string to save output file name
string inputFileName;   // string to save input file name

char func;
int imgFlag = 0;
int firstFlag = 0;

/*
Trans color system from hsv to RGB
*/

void SetHSV(double h, double s, double v, unsigned char color[3]) {
    double r, g, b;
    if(s==0)
        r = g = b = v;

    else {
        if(h==1) h = 0;
        double z = floor(h*6); int i = int(z);
        double f = double(h*6 - z);
        double p = v*(1-s);
        double q = v*(1-s*f);
        double t = v*(1-s*(1-f));

        switch(i){
        case 0: r=v; g=t; b=p; break;
        case 1: r=q; g=v; b=p; break;
        case 2: r=p; g=v; b=t; break;
        case 3: r=p; g=q; b=v; break;
        case 4: r=t; g=p; b=v; break;
        case 5: r=v; g=p; b=q; break;
        }
    }
    int c;
    c = int(256*r); if(c>255) c = 255; color[0] = c;
    c = int(256*g); if(c>255) c = 255; color[1] = c;
    c = int(256*b); if(c>255) c = 255; color[2] = c;
}

/*
Get complex function
*/
complex<double> fun(complex<double>& c , char func){
    const complex<double> i(0., 1.);
    
    switch(func){
		case 'a':
			return c;
		case 'b':
			return pow(E,c);
		case 'c':
			return (pow(c,2) -1.) *pow(c -2. -i, 2) /(pow(c,2) +2. +2. *i);
		case 'd':
			return (pow(c,2)-1.)*(pow(c,2)+1.);
		case 'e':
			return 1./c;
		case 'f':
			return c+i/c;
		case 'h':
			return (1./pow(i*c,18)-1./(i*c))/(1./(i*c)-1.);
		case 'i':
			return log(c);
		case 'j':
			return sin(c);
		case 'k':
			return asin(c);
		case 'l':
			return cos(c);
		case 'm':
			return acos(c);
		case 'n':
			return tan(c);
		case 'o':
			return atan(c);
		case 'p':
			return sinh(c);
		case 'q':
			return cosh(c);
		case 'r':
			return tanh(c);
		//case 'g':
		//	return tgammaf(c);
		default:
			return c;
	
	}
    //std::tgamma(c);

}
/*
Create domain coloring image
*/
void colorDomainImage(){
	//cout << "colorDomainImage start" << endl;
	int dimx = DEFAULTWIDTH; 
	int dimy = DEFAULTHEIGHT;
	
	pixmapOut = new Pixel*[dimy];
	if(pixmapOut != NULL)
		pixmapOut[0] = new Pixel[dimx * dimy];
	for(int i = 1; i < dimy; i++)
		pixmapOut[i] = pixmapOut[i - 1] + dimx;
	
	double reMin = -2;
	double reMax =  2;
    double imMin = -2;
    double imMax =  2;
    
    unsigned char RGBcolor[3];
    
    int row,col;
    for(row=0;row<dimy;++row){
        double im = imMax - (imMax-imMin)*row/(dimy-1);
        for(col=0;col<dimx;++col){            
            double re = reMax - (reMax-reMin)*col/(dimx-1);
            complex<double> c(re, im); 
            complex<double> v = fun(c, func);   //
            
            double a = arg(v); // get the angle of v
            
            while(a<0) a += 2*PI;
            a = a/(2*PI); // get h ([0:1])
            double m = abs(v);
            double ranges = 0;
            double rangee = 1;

            while(m>rangee){
                ranges = rangee;
                rangee *= E;
            }
            
            double k   = (m-ranges)/(rangee-ranges);
            double sat = k < 0.5 ? k *2: 1 -(k -0.5) *2;
            sat = 1 - pow(1-sat, 3); sat = 0.4 + sat*0.6;

            
            double val = k < 0.5 ? k *2: 1 -(k -0.5) *2; val = 1 - val;
            val = 1 - pow(1-val, 3); val = 0.6 + val*0.4;

            unsigned char color[3];
            SetHSV(a,sat,val,color);          
            
            pixmapOut[row][col].r = color[0];
			pixmapOut[row][col].g = color[1];
			pixmapOut[row][col].b = color[2];
			pixmapOut[row][col].a = 255;
			
        }                     
    }
    //cout << "colorDomainImage End" << endl;

}

/*
Create warping image
*/
void imgDomainImage(){
	//cout << "imgDomainImage start" << endl;
	int dimx = DEFAULTWIDTH; 
	int dimy = DEFAULTHEIGHT;
	
	pixmapOut = new Pixel*[dimy];
	if(pixmapOut != NULL)
		pixmapOut[0] = new Pixel[dimx * dimy];
	for(int i = 1; i < dimy; i++)
		pixmapOut[i] = pixmapOut[i - 1] + dimx;
	
	double reMin = -2;
	double reMax =  2;
    double imMin = -2;
    double imMax =  2;
    
    unsigned char RGBcolor[3];
    
    int row,col;
    
    float xMax=0, xMin=0, yMax=0, yMin=0;
    float xArray[dimy];
    float yArray[dimx];
    

    
    for(row=0;row<dimy;++row){
        double im = imMax - (imMax-imMin)*row/(dimy-1);
        for(col=0;col<dimx;++col){            
            double re = reMax - (reMax-reMin)*col/(dimx-1);
            complex<double> c(re, im); 
            complex<double> v = fun(c ,func);   //
            
            double a = arg(v); // get the angle of v
            
            while(a<0) a += 2*PI;
            a = a/(2*PI); // get h ([0:1])
            double m = abs(v);
            double ranges = 0;
            double rangee = 1;

            while(m>rangee){
                ranges = rangee;
                rangee *= E;
            }

            
            double k   = (m-ranges)/(rangee-ranges);
            double sat = k < 0.5 ? k *2: 1 -(k -0.5) *2;
            sat = 1 - pow(1-sat, 3); sat = 0.4 + sat*0.6;      
            
            float r = m;
            float theta = a*2*PI;

            
            }
            
            xMax = MAX(xMax,r*cos(theta));
            xMin = MIN(xMin,r*cos(theta));
            
            yMax = MAX(yMax,r*sin(theta));
            yMin = MIN(yMin,r*sin(theta));
            
       }     
  
        for(row=0;row<dimy;++row){    
			double im = imMax - (imMax-imMin)*row/(dimy-1);
			for(col=0;col<dimx;++col){    
				double re = reMax - (reMax-reMin)*col/(dimx-1);
				complex<double> c(re, im); 
				complex<double> v = fun(c, func);   //
            
				double a = arg(v); // get the angle of v
            
				while(a<0) a += 2*PI;
				a = a/(2*PI); // get h ([0:1])
				double m = abs(v);
				double ranges = 0;
				double rangee = 1;

				while(m>rangee){
					ranges = rangee;
					rangee *= E;
				}

            
				double k   = (m-ranges)/(rangee-ranges);
				double sat = k < 0.5 ? k *2: 1 -(k -0.5) *2;
				sat = 1 - pow(1-sat, 3); sat = 0.4 + sat*0.6;

				//cout << "a: " << a << ", sat: " << sat << endl;          
            
				float r = m;
				float theta = a*2*PI;
	
				
				int x = (r*cos(theta)-xMin)/(xMax-xMin)*dimy;
				int y = (r*sin(theta)-yMin)/(yMax-yMin)*dimx;
				int tempX, tempY;
				//cout << "x: " << x << ", y: " << y << endl; 

				int count = 0;
				while((x>=ImWidth or y >= ImHeight or x<0 or y<0) and count < 10000){
					count = count + 1;
					if(x>ImWidth){
						x = ImWidth - (x-ImWidth);
					}
					if(y >= ImHeight){
						y = ImHeight - (y-ImHeight);
					}
					if(x<0){
						x = -x;
					}
					if(y<0){
						y = -y;
					}
				}
			
				if(x>=ImWidth or y >= ImHeight or x<0 or y<0 ){
					pixmapOut[row][col].r = 0;
					pixmapOut[row][col].g = 0;
					pixmapOut[row][col].b = 0;
					pixmapOut[row][col].a = 0;
				}else{
					pixmapOut[row][col].r = pixmap[y][x].r;
					pixmapOut[row][col].g = pixmap[y][x].g;
					pixmapOut[row][col].b = pixmap[y][x].b;
					pixmapOut[row][col].a = pixmap[y][x].a;
				}
			}
		}          
    }	
}

/*
Print all function options user can use.
*/

void functionPrint(){
	system("clear");
	cout << "Function List:" << endl;
	cout << "(a) sin(z)" << endl;
	cout << "(b) pow(E,z)" << endl;
	cout << "(c) (pow(z,2) -1.) *pow(z -2. -i, 2) /(pow(z,2) +2. +2. *i);" << endl;
	cout << "(d) (pow(z,2)-1.)*(pow(z,2)+1.);" << endl;
	cout << "(e) 1./z" << endl;
	cout << "(f) z+i/z" << endl;
	cout << "(g) gamma(z)" << endl;
	cout << "(h) (1./pow(i*c,18)-1./(i*c))/(1./(i*c)-1.)" << endl;
	cout << "(i) log(z)" << endl;
	cout << "(j) sin(z)" << endl;
	cout << "(k) asin(z)" << endl;
	cout << "(l) cos(z)" << endl;
	cout << "(m) acos(z)" << endl;
	cout << "(n) tan(z)" << endl;
	cout << "(o) atan(z)" << endl;
	cout << "(p) sinh(z)" << endl;
	cout << "(q) cosh(z)" << endl;
	cout << "(r) tanh(z)" << endl;

}

/*
Read user input
*/
void read_input() {
	functionPrint();

	cout << "Select a function" << endl;
	cin >> func;
	
}

//
//  Routine to cleanup the memory.   
//
void destroy(){
 if (pixmap){
     delete pixmap[0];
	 delete pixmap;  
  }
  if (pixmapOut){
     delete pixmapOut[0];
	 delete pixmapOut;  
  }
}
//
//  Routine to read an image file and store in a pixmap
//  returns the size of the image in pixels if correctly read, or 0 if failure
//
int readImage(string infilename){
  
  if(firstFlag ==1){
		cout << "Enter input file name" << endl;
		cin >> infilename;
  }
  
  firstFlag = 1;
  // Create the oiio file handler for the image, and open the file for reading the image.
  // Once open, the file spec will indicate the width, height and number of channels.
  std::unique_ptr<ImageInput> infile = ImageInput::open(infilename);
  if(!infile){
    cerr << "Could not input image file " << infilename << ", error = " << geterror() << endl;
    return 0;
  }

  // Record image width, height and number of channels in global variables
  ImWidth = infile->spec().width;
  ImHeight = infile->spec().height;
  ImChannels = infile->spec().nchannels;

 
  // allocate temporary structure to read the image 
  unsigned char tmp_pixels[ImWidth * ImHeight * ImChannels];

  // read the image into the tmp_pixels from the input file, flipping it upside down using negative y-stride,
  // since OpenGL pixmaps have the bottom scanline first, and 
  // oiio expects the top scanline first in the image file.
  int scanlinesize = ImWidth * ImChannels * sizeof(unsigned char);
  if(!infile->read_image(TypeDesc::UINT8, &tmp_pixels[0] + (ImHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
    cerr << "Could not read image from " << infilename << ", error = " << geterror() << endl;
    return 0;
  }
 
 // get rid of the old pixmap and make a new one of the new size
  destroy();
  
 // allocate space for the Pixmap (contiguous approach, 2d style access)
  pixmap = new Pixel*[ImHeight];
  if(pixmap != NULL)
	pixmap[0] = new Pixel[ImWidth * ImHeight];
  for(int i = 1; i < ImHeight; i++)
	pixmap[i] = pixmap[i - 1] + ImWidth;
 
 //  assign the read pixels to the the data structure
 int index;
  for(int row = 0; row < ImHeight; ++row) {
    for(int col = 0; col < ImWidth; ++col) {
		index = (row*ImWidth+col)*ImChannels;
		
		if (ImChannels==1){ 
			pixmap[row][col].r = tmp_pixels[index];
			pixmap[row][col].g = tmp_pixels[index];
			pixmap[row][col].b = tmp_pixels[index];
			pixmap[row][col].a = 255;
		}
		else{
			pixmap[row][col].r = tmp_pixels[index];
			pixmap[row][col].g = tmp_pixels[index+1];
			pixmap[row][col].b = tmp_pixels[index+2];			
			if (ImChannels <4) // no alpha value is present so set it to 255
				pixmap[row][col].a = 255; 
			else // read the alpha value
				pixmap[row][col].a = tmp_pixels[index+3];			
		}
    }
  }
  
  cout << "End of read image" << endl;
 
  // close the image file after reading, and free up space for the oiio file handler
  infile->close();
  
  // set the pixel format to GL_RGBA and fix the # channels to 4  
  pixformat = GL_RGBA;  
  ImChannels = 4;
 
  // return image size in pixels
  return ImWidth * ImHeight;
}


//
// Routine to display a pixmapOut of warping result in the current window
//
void displayImage1(){
  // if the window is smaller than the image, scale it down, otherwise do not scale
  if(WinWidth < OutWidth  || WinHeight < OutHeight)
    glPixelZoom(float(VpWidth) / OutWidth, float(VpHeight) / OutHeight);
  else
    glPixelZoom(1.0, 1.0);
 
  // display starting at the lower lefthand corner of the viewport
  //glRasterPos2i(0, 0);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glDrawPixels(OutWidth, OutHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixmapOut[0]);
}

//
// Routine to display a pixmap of input image in the current window
//
void displayImage2(){
  // if the window is smaller than the image, scale it down, otherwise do not scale
  if(WinWidth < ImWidth  || WinHeight < ImHeight)
    glPixelZoom(float(VpWidth) / ImWidth, float(VpHeight) / ImHeight);
  else
    glPixelZoom(1.0, 1.0);
  
  // display starting at the lower lefthand corner of the viewport
  //glRasterPos2i(0, 0);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glDrawPixels(ImWidth, ImHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixmap[0]);
  	
}


//
// Routine to write the current framebuffer to an image file
//
void writeImage(string outfilename){
  // make a pixmap that is the size of the window and grab OpenGL framebuffer into it
  // alternatively, you can read the pixmap into a 1d array and export this 
   ImChannels = 4;
   
   unsigned char local_pixmap[WinWidth * WinHeight * ImChannels];
   glReadPixels(0, 0, WinWidth, WinHeight, GL_RGBA, GL_UNSIGNED_BYTE, local_pixmap);
  
  // create the oiio file handler for the image
  std::unique_ptr<ImageOutput> outfile = ImageOutput::create(outfilename);
  if(!outfile){
    cerr << "Could not create output image for " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // Open a file for writing the image. The file header will indicate an image of
  // width WinWidth, height WinHeight, and ImChannels channels per pixel.
  // All channels will be of type unsigned char
  ImageSpec spec(WinWidth, WinHeight, ImChannels, TypeDesc::UINT8);
  if(!outfile->open(outfilename, spec)){
    cerr << "Could not open " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // Write the image to the file. All channel values in the pixmap are taken to be
  // unsigned chars. While writing, flip the image upside down by using negative y stride, 
  // since OpenGL pixmaps have the bottom scanline first, and oiio writes the top scanline first in the image file.
  int scanlinesize = WinWidth * ImChannels * sizeof(unsigned char);
  if(!outfile->write_image(TypeDesc::UINT8, local_pixmap + (WinHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
    cerr << "Could not write image to " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // close the image file after the image is written and free up space for the
  // ooio file handler
  outfile->close();
}

//
//   Display Callback Routine: clear the screen and draw the current image
//
void handleDisplay1(){
  

  // specify window clear (background) color to be opaque black
  glClearColor(0, 0, 0, 1);
  // clear window to background color
  glClear(GL_COLOR_BUFFER_BIT);  
  
  // only draw the image if it is of a valid size
  if(OutWidth > 0 && OutHeight > 0)
    displayImage1();
  
  // flush the OpenGL pipeline to the viewport
  glFlush();
}

//
//   Display Callback Routine: clear the screen and draw the current image
//
void handleDisplay2(){
  // specify window clear (background) color to be opaque black
  glClearColor(0, 0, 0, 1);
  // clear window to background color
  glClear(GL_COLOR_BUFFER_BIT);  
  
  // only draw the image if it is of a valid size
  //if(ImWidth > 0 && ImHeight > 0)
    displayImage2();
  
  // flush the OpenGL pipeline to the viewport
  glFlush();
}

//
//  Keyboard Callback Routine: 'r' - read and display a new image,
//  'w' - write the current window to an image file, 'q' or ESC - quit
//
void handleKey(unsigned char key, int x, int y){
  string infilename, outfilename;
  int ok;
  
  switch(key){
    case 'r':		// 'r' - read an image from a file
    case 'R':
      if(inputFileName==NULL){
          cout << "Input image filename? ";	  // prompt user for input filename
          cin >> inputFileName;
      }
      ok = readImage(inputFileName);
      if(ok)
        glutReshapeWindow(ImWidth, ImHeight); // OpenGL window should match new image
      glutPostRedisplay();
      break;
      
    case 'w':		// 'w' - write the image to a file
    case 'W':
      if(outputFileName==NULL){
      	cout << "Output image filename? ";  // prompt user for output filename
      	cin >> outputFileName;
      }
      writeImage(outputFileName);
      break;  
 
 	case 'f':		// 
    case 'F':
		//destroy();
		read_input() ;
		if(imgFlag == 1){
			imgDomainImage();
		}else{
			colorDomainImage();
		}
		glutPostRedisplay();
		
	
	case 'q':		// q or ESC - quit
    case 'Q':
    case 27:
      destroy();
      exit(0);
      
    default:		// not a valid key -- just ignore it
      return;
  }
}


//
//  Reshape Callback Routine: If the window is too small to fit the image,
//  make a viewport of the maximum size that maintains the image proportions.
//  Otherwise, size the viewport to match the image size. In either case, the
//  viewport is centered in the window.
//
void handleReshape1(int w, int h){
  float imageaspect = (float)ImWidth / (float)ImHeight;	// aspect ratio of image
  float newaspect = (float)w / (float)h; // new aspect ratio of window
  
  // record the new window size in global variables for easy access
  WinWidth = w;
  WinHeight = h;
  
  // if the image fits in the window, viewport is the same size as the image
  if(w >= OutWidth && h >= OutHeight){
    Xoffset = (w - OutWidth) / 2;
    Yoffset = (h - OutHeight) / 2;
    VpWidth = OutWidth;
    VpHeight = OutHeight;
  }
  // if the window is wider than the image, use the full window height
  // and size the width to match the image aspect ratio
  else if(newaspect > imageaspect){
    VpHeight = h;
    VpWidth = int(imageaspect * VpHeight);
    Xoffset = int((w - VpWidth) / 2);
    Yoffset = 0;
  }
  // if the window is narrower than the image, use the full window width
  // and size the height to match the image aspect ratio
  else{
    VpWidth = w;
    VpHeight = int(VpWidth / imageaspect);
    Yoffset = int((h - VpHeight) / 2);
    Xoffset = 0;
  }
  
  // center the viewport in the window
  glViewport(Xoffset, Yoffset, VpWidth, VpHeight);
  
  // viewport coordinates are simply pixel coordinates
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, VpWidth, 0, VpHeight);
  glMatrixMode(GL_MODELVIEW);
}

void handleReshape2(int w, int h){
  float imageaspect = (float)ImWidth / (float)ImHeight;	// aspect ratio of image
  float newaspect = (float)w / (float)h; // new aspect ratio of window
  
  // record the new window size in global variables for easy access
  WinWidth = w;
  WinHeight = h;
  
  // if the image fits in the window, viewport is the same size as the image
  if(w >= ImWidth && h >= ImHeight){
    Xoffset = (w - ImWidth) / 2;
    Yoffset = (h - ImHeight) / 2;
    VpWidth = ImWidth;
    VpHeight = ImHeight;
  }
  // if the window is wider than the image, use the full window height
  // and size the width to match the image aspect ratio
  else if(newaspect > imageaspect){
    VpHeight = h;
    VpWidth = int(imageaspect * VpHeight);
    Xoffset = int((w - VpWidth) / 2);
    Yoffset = 0;
  }
  // if the window is narrower than the image, use the full window width
  // and size the height to match the image aspect ratio
  else{
    VpWidth = w;
    VpHeight = int(VpWidth / imageaspect);
    Yoffset = int((h - VpHeight) / 2);
    Xoffset = 0;
  }
  
  // center the viewport in the window
  glViewport(Xoffset, Yoffset, VpWidth, VpHeight);
  
  // viewport coordinates are simply pixel coordinates
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, VpWidth, 0, VpHeight);
  glMatrixMode(GL_MODELVIEW);
}

/*
   Main program to read an image file, then ask the user
   for transform information, transform the image and display
   it using the appropriate warp.  Optionally save the transformed
   images in  files.
*/
int main(int argc, char *argv[]){
	
   //your code to read in the input image
     // set up the default window and empty pixmap if no image or image fails to load
    WinWidth = DEFAULTWIDTH;
    WinHeight = DEFAULTHEIGHT;
    ImWidth = 0;
    ImHeight = 0;
   
	//pixmap = new Pixel*[ImHeight];
  
    // load the image if present, and size the window to match
    if(argc >= 2){
		imgFlag = 1;
		inputFileName = argv[1];
		readImage(inputFileName);
		read_input();
		imgDomainImage();
		
	}else{
		imgFlag = 0;
		read_input();
		colorDomainImage();
	}	

	if(argc >= 3){
		outputFileName = argv[2];
		writeImage(outputFileName);
	}

	// start up GLUT
	glutInit(&argc, argv);
  
	GLint window1, window2;
  
	// create the graphics window, giving width, height, and title text
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowSize(OutWidth, OutHeight);
	glutInitWindowPosition(0,0);
	window1 = glutCreateWindow("Output image");
  
	// set up the callback routines
	glutDisplayFunc(handleDisplay1); // display update callback
	glutKeyboardFunc(handleKey);	  // keyboard key press callback
	glutReshapeFunc(handleReshape1); // window resize callback
  
  	if(argc >= 2){
		// create the graphics window, giving width, height, and title text
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
		glutInitWindowSize(ImWidth, ImHeight);
		glutInitWindowPosition(620,0);
		window2 = glutCreateWindow("Image Viewer");
  
		// set up the callback routines
		glutDisplayFunc(handleDisplay2); // display update callback
		glutKeyboardFunc(handleKey);	  // keyboard key press callback
		glutReshapeFunc(handleReshape2); // window resize callback
	}
  
	// Enter GLUT's event loop
	glutMainLoop();

   return 0;
}



