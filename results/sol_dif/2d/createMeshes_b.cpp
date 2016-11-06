#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;

int main()
{
	double L=1.00;    			// Global scale factor
	double lc = .1; //= L/1;				// Factor for the mesh refinement solid interface
	double li = .1; //L/1;
	double lw = .5; //L/1;
	double width  = 6; //40*L;		// Tube width
	double height = 6; //40*L;			// Tube height
//	double R = 0.2*L; 			// Radius of the bubbles
	int    N = 2;				// Number of bubbles in the tube
	double coord[N][2]; 			// Coordenadas centrais das bolhas
	double *ri = new double[N];   //Raio de cada bolha menor
	double *rI = new double[N];   //Raio mayor
	double ri_min = .125, ri_max = .125; // L_min =1.2, L_max = 2;  //
	double rI_min = .125, rI_max = .200;
    int    FS = 101;                  //physical fsi starting tag
	double rho = 1.5;
	double pi = 3.14159265358979323846;
	double layer = 0.2*width, min, dist, epsil, epsiy, theta;
	
	// Arquivo de saida
	fstream fout;
	fout.open("Nsol1.geo",std::fstream::out);
	
	// mass radious information file
	fstream iout;
	iout.open("Nsol1.txt",std::fstream::out);

	cout << "intervalo da rand: " <<  RAND_MAX << endl;
    srand( (unsigned)time(NULL) );
	for( int i=0; i<N; i++)
	{
		coord[i][0] = (width-2*layer)*rand()/((double)RAND_MAX) + layer;
		coord[i][1] = (height-2*layer)*rand()/((double)RAND_MAX) + layer;
		
        cout << "Numero :" << i << " " << coord[i][0] << " " << coord[i][1] << endl;
        
        ri[i] = (ri_max-ri_min)*rand()/((double)RAND_MAX) + ri_min;
		rI[i] = (rI_max-rI_min)*rand()/((double)RAND_MAX) + rI_min;
		
        for( int j=0; j<i; j++)
        {
			dist = sqrt(pow(coord[i][0]-coord[j][0], 2)+pow(coord[i][1]-coord[j][1], 2));
			min = (fmax(ri[i],rI[i]) + fmax(ri[j],rI[j]))*1.2;
			int it = 0;
			while( dist < min && it < 100)
			{
				ri[i] = 0.5*ri[i];
				rI[i] = 0.5*rI[i];
				min = (fmax(ri[i],rI[i]) + fmax(ri[j],rI[j]))*1.2;
				dist = sqrt(pow(coord[i][0]-coord[j][0], 2)+pow(coord[i][1]-coord[j][1], 2));
				it++;
			}
		}

		cout << "Raio : " << ri[i] << " " << rI[i] << endl;
		if( (ri[i] <ri_min) || (rI[i] <rI_min) )
		    i--;
	}
	
	//double min = 1;
	//for( int i=0; i<N ; i++) if(min > ri[i]) min=ri[i];
	//lc = min/2;
	
	
	fout << "L = " << L << ";" << endl;    				// Global scale factor
	fout << "lc = " << lc << ";" << endl;				// Factor for the mesh refinement
	fout << "li = " << li << ";" << endl;
	fout << "lw = " << lw << ";" << endl;
	fout << "width = " << width << ";"<< endl;			// Tube width
	fout << "height = " << height << ";"<< endl;			// Tube height
	fout << "N = "<< N << ";"<< endl;					// Number of bubbles in the tube
	for( int i=0; i<N ; i++)
	{
		fout << "Point(" << i+1 << ") = {" << coord[i][0] << "," << coord[i][1]  << ",0,lc};" << endl;
	}
	
	// Definição da borda dos corpos
	for( int i=0; i<N ; i++)
	{
		fout << "Point(" << N+6*i+1 << ") = {" << coord[i][0]+rI[i] << "," << coord[i][1]       << ",0,lc};" << endl;
		fout << "Point(" << N+6*i+2 << ") = {" << coord[i][0]       << "," << coord[i][1]+ri[i] << ",0,lc};" << endl;
		fout << "Point(" << N+6*i+3 << ") = {" << coord[i][0]-rI[i] << "," << coord[i][1]       << ",0,lc};" << endl;
		fout << "Point(" << N+6*i+4 << ") = {" << coord[i][0]       << "," << coord[i][1]-ri[i] << ",0,lc};" << endl;
		epsil = 0.8*rI[i];
		epsiy = (ri[i]/rI[i])*sqrt(pow(rI[i],2)-pow(epsil,2));
		fout << "Point(" << N+6*i+5 << ") = {" << coord[i][0]-epsil << "," << coord[i][1]+epsiy << ",0,lc};" << endl;
		fout << "Point(" << N+6*i+6 << ") = {" << coord[i][0]-epsil << "," << coord[i][1]-epsiy << ",0,lc};" << endl;
		
		theta = 2*pi*rand()/((double)RAND_MAX);
		fout << "Rotate {{0,0,1},{" << coord[i][0] << "," << coord[i][1] << ",0}," << theta << "} {\n" 
		     << "\t Point{" << N+6*i+1 << "," << N+6*i+2 << "," << N+6*i+5 << 
		                "," << N+6*i+3 << "," << N+6*i+6 << "," << N+6*i+4 << "};" << "\n}";

		fout << "Ellipse(" << 6*i+1 << ") = {" << N+6*i+1 << "," <<  i+1 <<"," << N+6*i+4 << "," << N+6*i+2 << "};"<< endl;
		fout << "Ellipse(" << 6*i+2 << ") = {" << N+6*i+2 << "," <<  i+1 <<"," << N+6*i+4 << "," << N+6*i+5 << "};"<< endl;
		fout << "Ellipse(" << 6*i+3 << ") = {" << N+6*i+5 << "," <<  i+1 <<"," << N+6*i+4 << "," << N+6*i+3 << "};"<< endl;
		fout << "Ellipse(" << 6*i+4 << ") = {" << N+6*i+3 << "," <<  i+1 <<"," << N+6*i+2 << "," << N+6*i+6 << "};"<< endl;
		fout << "Ellipse(" << 6*i+5 << ") = {" << N+6*i+6 << "," <<  i+1 <<"," << N+6*i+2 << "," << N+6*i+4 << "};"<< endl;
		fout << "Ellipse(" << 6*i+6 << ") = {" << N+6*i+4 << "," <<  i+1 <<"," << N+6*i+2 << "," << N+6*i+1 << "};"<< endl;


		fout << "Line Loop(" << i+1 << ") = {" << 6*i+1 << "," << 6*i+2 << "," << 6*i+3 << "," 
		                                       << 6*i+4 << "," << 6*i+5 << "," << 6*i+6 << "};" << endl;
        fout << "Physical Line(" << FS+i << ") = {" << 6*i+1 << "," << 6*i+2 << "," 
                                                    << 6*i+5 << "," << 6*i+6 << "};" << endl;  //fsi_tags
        fout << "Physical Line(" <<FS+i+N<< ") = {" << 6*i+3 << "," << 6*i+4 << "};" << endl;  //fsi_tags                                          
        //fout << "Physical Point(" << FS+i << ") = {" << N+2*i+1 << "," << N+2*i+2 << "};" << endl; //interfacial points
		iout << rho*pi*ri[i]*rI[i] << " " << fmax(ri[i],rI[i]) << " " << pi*ri[i]*rI[i] << " " 
		     << coord[i][0] << " " << coord[i][1] << endl;
	}
	
	// Tags fisicas
//	fout << "physicalTagWall = 200;"<< endl;
//	fout << "physicalTagSurfaceBubble = 1;" << endl;
//	fout << "physicalTagOutsideBubble = 201;" << endl;
//	fout << "physicalTagInsideBubble = 100;"<< endl;    
//	fout << "Physical Line(physicalTagSurfaceBubble) = {1:2*N};" << endl;

	// Geometria do tubo
	fout << "// --------------- Defining the geometry of tube -----------------------"<< endl;
	fout << "// Points"<< endl;
	fout << "Point(7*N+1) = {0, 0, 0, lw};"<< endl;
	fout << "Point(7*N+2) = {width, 0, 0, lw};"<< endl;
	fout << "Point(7*N+3) = {width, height, 0, lw};"<< endl;
	fout << "Point(7*N+4) = {0, height, 0, lw};"<< endl;

	fout << "// Edges"<< endl;
	fout << "Line(6*N+1) = {7*N+1, 7*N+2};"<< endl;
	fout << "Line(6*N+2) = {7*N+2, 7*N+3};"<< endl;
	fout << "Line(6*N+3) = {7*N+3, 7*N+4};"<< endl;
	fout << "Line(6*N+4) = {7*N+4, 7*N+1};"<< endl;

	fout << "Line Loop(N+1) = {6*N+1:6*N+4};"<< endl;

	fout << "Physical Line(1) = {6*N+4};" << endl;
	fout << "Physical Line(2) = {6*N+1};" << endl;
	fout << "Physical Line(3) = {6*N+2};" << endl;
	fout << "Physical Line(4) = {6*N+3};" << endl;
//	fout << "Physical Point(5) = {3*N+4, 3*N+3, 3*N+2, 3*N+1};" << endl;
//	fout << "Physical Line(physicalTagWall) = {2*N+1:2*N+4};"<< endl;

    // Superficies no plano que conterão a malha
	fout << "// ------------- Defining the surface mesh -----------------------------"<< endl;
	fout << "// Fluid outside the bubbles"<< endl;
	fout << "Plane Surface(1) = {N+1,1:N};"<< endl;
	fout << "Physical Surface(10) = {1};"<< endl;

    
	fout << "// Fluid inside the bubbles"<< endl;
	fout << "For t In {1:N}"<< endl;
	fout << "Plane Surface(t+1) = {t};"<< endl;
	fout << "Physical Surface(" << FS << "+2*N-1+t) = {t+1};"<< endl;
	fout << "EndFor"<< endl;
//	fout << "Physical Surface(physicalTagInsideBubble) = {2:N+1};"<< endl;
	
	fout.close();
	iout.close();
	cout << lc << endl;

	return 0;
}
