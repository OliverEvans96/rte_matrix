// kelp_model.cpp
// Oliver Evans
// Clarkson University REU 2016
// Created: Wed 06 Jul 2016 01:11:25 PM EDT
// Last E_dited: Tue 26 Jul 2016 10:49:56 AM EDT

// Compile & run with:
// g++ -fopenmp kelp_model.cpp -o kelp_model.out && time ./kelp_model.out &

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

///////////////
// CONSTANTS //
///////////////

// Pi
const double PI = 3.141592653589793;

// Frond shape & ratio parameters
const double FS = 0.5;
const double FR = 2.0;
const double ALPHA = atan((1+FS)/(2*FR*FS));

// Vertical growth density (fronds/meter)
const double RHO = 10;

// Index of refraction of water
const double N_REFRACT = 3.0/4.0;


////////////////
// PROTOTYPES //
////////////////

// MATH FUNCTIONS //
double sign(double xx);
double min(double xx, double yy);
double max(double xx, double yy);
double factorial(int nn);
int mod(int aa,int bb);
double I0(double xx,int nn);
double I0(double xx);
inline double vonMisesPDF(double xx,double kk,double mu);

// POLAR COORDINATE FUNCTIONS //
double theta_xy(double xx,double yy);
void polar_shift(double theta,double rr,double dz,double phi_s,double theta_s,double &theta_hat,double &rr_hat);
inline double angleDiff(double theta,double phi,double theta_prime,double phi_prime);

// NUMERICAL METHODS //
double csimp(double (*ff)(double),double aa,double bb,int nn);
double csimp(double (*ff)(double,double),double aa,double bb,int nn,double pp);
double csimp(double (*ff)(double,double,double),double aa,double bb,int nn,double pp1,double pp2);
double csimp(double (*ff)(double,double*,double*,int),double aa,double bb,int nn,double* pp1,double* pp2,int pp3);

// INTERPOLATION //
double linear(double xx,double x0,double x1,double y0,double y1);
double bilinear(double* xx,double* yy,double** zz,int nn,int mm,double xval,double yval);
double vsf(double xx,double* vsf_th,double* vsf_vals,int nVSF);
double vsf1(double xx,double* vsf_th,double* vsf_vals,int nVSF);

// UTILITY FUNCTIONS //
int countLines(ifstream &inFile);
int binarySearch(double xx,double* vv,int first,int last);
inline int findBin(double aa, double bb,int nn,double xx);
inline int upperBin(double aa, double bb,int nn,double xx);
void read2column(ifstream &inFile,double* col1,double* col2,int numLines);
void print3d(double*** v, int nx, int ny,int nz, const char *name);
void write3d(double*** v, int nx, int ny,int nz, const char *name,ofstream &outFile);
void print2d(double** v, int nx, int ny, const char *name);
void write2d(double** v, int nx, int ny, const char *name,ofstream &outFile);
void write2d(double** v, int nx, int ny, const char *name,int index,ofstream &outFile);
void print1d(double* v, int nn, const char *name);
void contourf(double* theta,double* rr,double** zz,int ntheta,int nr);

// MODEL-SPECIFIC FUNCTIONS //
inline double P_theta_f(double xx,double v_w,double theta_w);
inline double L_min_shade(double theta_p,double r_p,double theta_f);
inline double theta_min_shade(double theta_p,double r_p,double LL);
inline double theta_max_shade(double theta_p,double r_p,double LL);
double N_shade_3d(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L);
double P_shade_2d(double theta_p,double r_p,double v_w,double theta_w,int nLBin,double* LBin_vals,double* P_L);
double availableLight(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double E_d0,double attenuationCoef,double a_k);
inline void frondTransform(double ss,double tt,double &xx,double &yy,double theta_f,double LL);
inline double JJ(double ss,double tt,double LL);
double calculateLightAbsorption(double theta_f,double LL,double z_f,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double E_d0,double attenuationCoef,double a_k);
inline double newLength(double lightAbsorbed,double LL,double tt);
void recalculateLengthDistribution(double tt,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double dLBin,double LBin_min,double LBin_max,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double E_d0,double attenuationCoef,double a_k);


//////////
// MAIN //
//////////

int main()
{
    // Bin convention:
    // Values stored in *Bin_vals are closed lower edges of bins
    // Values stored in P_L[ii][jj] are the proportion of the
    // population of fronds whose depth is in the interval
    // [ zBin_vals[ii] , zBin_vals[ii] )
    // whose length is in the interval
    // [ LBin_vals[ii] , LBin_vals[ii] )
    // (Upper half open intervals)

    // n*Bin is the number of bins (not number of edges)

    cout << "Start!" << endl;

    ////////////////////////
    // READ VSF FROM FILE //
    ////////////////////////

    // Open stream
    ifstream vsfFile("../data/vsf/nuc_vsf.txt");

    // Count lines
    int nVSF = countLines(vsfFile)-1;

    // Skip first line
    string junk;
    getline(vsfFile,junk);

    // Allocate pointers
    double* vsf_th = new double[nVSF];
    double* vsf_vals = new double[nVSF];

    // Read arrays
    read2column(vsfFile,vsf_th,vsf_vals,nVSF);

    // Close stream
    vsfFile.close();

    ///////////////////////////
    // SIMULATION PARAMETERS //
    ///////////////////////////

    // Time
    double tt = 0;
    double dt = 1;
    double t_max = 1;
    int nTimeSteps = t_max/dt;

    // Space
    double xmin = -5;
    double xmax = 5;
    double dx = 1;

    double ymin = -5;
    double ymax = 5;
    double dy = 1;

    double zmin = 0;
    double zmax = 5;
    double dz = 1;

    // Zenith Angle
    double phimin = 0;
    double phimax = PI;
    double dphi = PI/10;

    // Azimuth Angle
    double thetamin = -PI;
    double thetamax = PI;
    double dtheta = PI/5;

    // Depth
    double dzBin = 0.5;
    double zBin_min = 0;
    double zBin_max = 10;

    // Length
    double dLBin = .5;
    double LBin_min = -dLBin;
    double LBin_max = 10;

    // Water Current
    double theta_w = 0;
    double v_w = 5;

    // Sun
    double phi_s_orig = .48*PI;
    double theta_s = 2*PI/3;
    // Refracted angle of sun
    double phi_s = asin(N_REFRACT*sin(phi_s_orig));

    // Light
    double E_d0 = 10;
    double attenuationCoef = 0.2;
    // Absorption by kelp
    double a_k = 0.7;
    // Absoroption by water
    double a_w = 0.3;
    // Scattering by water
    double b_w = csimp(vsf1,vsf_th[0],vsf_th[nVSF-1],100,vsf_th,vsf_vals,nVSF); //Integral of beta (volume scattering function)
    cout << "b_w = " << b_w << endl;

    printf("dx = %.2f\n",dx);
    printf("dy = %.2f\n",dy);
    printf("dz = %.2f\n",dz);
    printf("dtheta = %.2f\n",dtheta);
    printf("dphi = %.2f\n",dphi);
    cout << endl;


    /////////////////////////////
    // PARTITION SPACE, LENGTH //
    /////////////////////////////

    // Partition depth (units are meters)
    int nzBin = int(floor((zBin_max-zBin_min)/dzBin)+1); // Length of zBin_vals
    double* zBin_vals = new double[nzBin];
    for(int kk=0;kk<nzBin;kk++)
        zBin_vals[kk] = zBin_min + dzBin*kk;

    // Partition length (units are meters)
    // First bin [-dLBin,0] is for dead fronds
    int nLBin = int(floor((LBin_max-LBin_min)/dLBin)+1);
    double* LBin_vals = new double[nLBin];
    for(int ii=0;ii<nLBin;ii++)
        LBin_vals[ii] = LBin_min + dLBin*ii;

    // Spatial grid
    int nx = int(floor((xmax-xmin)/dx)+1);
    double* xx = new double[nx];
    for(int ii=0;ii<nx;ii++)
        xx[ii] = xmin + ii*dx;

    int ny = int(floor((ymax-ymin)/dy)+1);
    double* yy = new double[ny];
    for(int jj=0;jj<ny;jj++)
        yy[jj] = ymin + jj*dy;

    int nz = int(floor((zmax-zmin)/dz)+1);
    double* zz = new double[nz];
    for(int kk=0;kk<nz;kk++)
        zz[kk] = zmin + kk*dz;

    // Angular grid
    int ntheta = int(floor((thetamax-thetamin)/dtheta)+1);
    double* theta = new double[ntheta];
    for(int ll=0;ll<ntheta;ll++)
        theta[ll] = thetamin + ll*dtheta;

    int nphi = int(floor((phimax-phimin)/dphi)+1);
    double* phi = new double[nphi];
    for(int ll=0;ll<nphi;ll++)
        phi[ll] = phimin + ll*dphi;


    /////////////////////////////////////////////
    // ALLOCATE RADIANCE AND IRRADIANCE ARRAYS //
    /////////////////////////////////////////////

    // Indices
    // xx: 0 <= ii < nx
    // yy: 0 <= jj < ny
    // zz: 0 <= kk < nz
    // theta: 0 <= ll < ntheta
    // phi: 0 <= mm < nphi
    // tt: 0 <= nn < nTimeSteps

    // Radiance (L)
    double***** RR = new double****[nx];

    // Downward Irradiance (E_d)
    double*** E_d = new double**[nx];

    for(int ii=0;ii<nx;ii++)
    {
        RR[ii] = new double***[ny];
        E_d[ii] = new double*[ny];
        for(int jj=0;jj<ny;jj++)
        {
            RR[ii][jj] = new double**[nz];
            E_d[ii][jj] = new double[nz];
            for(int kk=0;kk<nz;kk++)
            {
                // Initialize downward irradiance to zero
                E_d[ii][jj][kk] = 0;
                RR[ii][jj][kk] = new double*[ntheta];
                for(int ll=0;ll<ntheta;ll++)
                {
                    RR[ii][jj][kk][ll] = new double[nphi];
                    for(int mm=0;mm<nphi;mm++)
                    {
                        // Radiance Initial guess
                        // ZEROED OUT
                        RR[ii][jj][kk][ll][mm] = 0*E_d0*exp(-attenuationCoef*zz[kk]);

                        // Surface BC
                        if(kk==0)
                            RR[ii][jj][kk][ll][mm] = 1;

                    }
                }
            }
        }
    }


    ///////////////////////////////////////////////
    // INITIALIZE POPULATION LENGTH DISTRIBUTION //
    ///////////////////////////////////////////////

    // Population distribution values
    // Each depth is treated as a separate population
    // Each row should sum to 1
    double** P_L = new double*[nzBin];
    double sum;
    for(int kk=0;kk<nzBin;kk++)
    {
        P_L[kk] = new double[nLBin];

        /*
        // Initialize population to all be as small as possible
        // (But not dead)
        for(int jj = 0;jj<nLBin;jj++)
            P_L[kk][jj] = 0;

        // All in smallest positive bin
        P_L[kk][1] = 1;
        */

        // Normal distribution with varying mean
        sum = 0;
        P_L[kk][0] = 0;
        for(int jj=1;jj<nLBin;jj++)
        {
            P_L[kk][jj] = exp(-pow(LBin_vals[jj]-5,2)/10);
            //P_L[kk][jj] = 1;
            sum += P_L[kk][jj];
        }

        //Normalize
        for(int jj=0;jj<nLBin;jj++)
        {
            P_L[kk][jj] /= sum;
        }
    }

    /*
    // Only fronds in top layer
    P_L[0][0] = 0;
    P_L[0][nLBin-1] = 1;
    */

    /*
    // And a layer in the middle
    P_L[5][0] = 0;
    P_L[5][nLBin-1] = 1;
    */

    /////////////////////////
    // INITIAL PYTHON CODE //
    /////////////////////////

    // Output to python
    ofstream outFile("../python/kelp_output.py");
    outFile << "import os" << endl;
    outFile << "from numpy import *" << endl;
    outFile << "import pickle" << endl;
    outFile << "os.chdir('../python')" << endl;
    outFile << "from vispy_volume import volume" << endl;
    outFile << "P_L = zeros([" << nTimeSteps+1 << "," << nzBin << "," << nLBin << "])" << endl;
    write2d(P_L,nzBin,nLBin,"P_L",0,outFile);


    ////////////////////////////
    // Test vsf interpolation //
    ////////////////////////////
    /*
    const int nvsftest = 10000;
    double dxtest = PI/nvsftest;
    double vsf_x[nvsftest];
    double vsf_y[nvsftest];
    for(int ii=0;ii<nvsftest;ii++)
        vsf_x[ii] = (ii+1)*dxtest;
    for(int ii=0;ii<nvsftest;ii++)
    {
        vsf_y[ii] = vsf(vsf_x[ii],vsf_th,vsf_vals,nVSF);
    }

    print1d(vsf_x,nvsftest,"vsf_x");
    print1d(vsf_y,nvsftest,"vsf_y");
    print1d(vsf_th,nvsftest,"vsf_th");
    print1d(vsf_vals,nvsftest,"vsf_vals");
    */


    //////////////////////////////////
    // VISUALIZE LIGHT AVAILABILITY //
    //////////////////////////////////
    /*

    // Allocate 3d arrays
    double*** xx_3d = new double**[nx];
    double*** yy_3d = new double**[nx];
    double*** zz_3d = new double**[nx];
    double*** PP_3d = new double**[nx];
    for(int ii=0;ii<nx;ii++)
    {
        xx_3d[ii] = new double*[ny];
        yy_3d[ii] = new double*[ny];
        zz_3d[ii] = new double*[ny];
        PP_3d[ii] = new double*[ny];
        for(int jj=0;jj<ny;jj++)
        {
            xx_3d[ii][jj] = new double[nzBin];
            yy_3d[ii][jj] = new double[nzBin];
            zz_3d[ii][jj] = new double[nzBin];
            PP_3d[ii][jj] = new double[nzBin];
        }
    }

    double theta_p,r_p,z_p;


    // Calculate expected number of plants shading at each point in 3d grid
    #pragma omp parallel for private(theta_p,r_p,z_p)
    for(int ii=0;ii<nx;ii++)
    {
        for(int jj=0;jj<ny;jj++)
        {
            theta_p = theta_xy(xx[ii],yy[jj]);
            r_p = sqrt(xx[ii]*xx[ii] + yy[jj]*yy[jj]);
            for(int kk=0;kk<nzBin;kk++)
            {
                //cout << "(" << ii << "," << jj << "," << kk << ")" << endl;
                z_p = zBin_vals[kk];

                //cout << "(" << ii << "," << jj << "," << kk << ")" << endl;
                xx_3d[ii][jj][kk] = xx[ii];
                yy_3d[ii][jj][kk] = yy[jj];
                zz_3d[ii][jj][kk] = z_p;
                //PP_3d[ii][jj][kk] = N_shade_3d(theta_p,r_p,zBin_vals[kk],v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L);
                // Amount of light NOT available - use reverse colormap
                PP_3d[ii][jj][kk] = E_d0 - availableLight(theta_p,r_p,z_p,v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L,E_d0,attenuationCoef,a_k);
            }
        }
    }

    theta_p = PI/2;
    r_p = 1;

    // Write variables to file
    outFile << "xlim = [" << xmin-dx/2 << "," << xmax+dx/2 << "]" << endl;
    outFile << "ylim = [" << ymin-dy/2 << "," << ymax+dy/2 << "]" << endl;
    outFile << "zlim = [" << -zBin_max-dzBin/2 << "," << -zBin_min+dzBin/2 << "]" << endl;
    outFile << "clim = array([0," << E_d0 << "])" << endl;
    write3d(xx_3d,nx,ny,nzBin,"xx_3d",outFile);
    write3d(yy_3d,nx,ny,nzBin,"yy_3d",outFile);
    write3d(zz_3d,nx,ny,nzBin,"zz_3d",outFile);
    write3d(PP_3d,nx,ny,nzBin,"PP_3d",outFile);

    // Plot volume
    outFile << "volume(xlim,ylim,zlim,PP_3d,clim)" << endl;


    */
    ////////////////////////
    // SOLVE RADIANCE PDE //
    ////////////////////////

    // Maximum number of SOR iterations
    int maxiter = 10;
    int depthLayer;
    double rr_val,th_val;
    double P_k = 0;
    double c_ijk;
    double R_star;
    double dRdx,dRdy,dRdz;
    double anglediff;
    double R_new;
    double residual;

    // Loop through iterations
    for(int iter=0;iter<maxiter;iter++)
    {
        // Reset residual
        residual = 0;
        // z
        for(int kk=1;kk<nz;kk++)
        {
            // Determine depth layer
            depthLayer = findBin(zBin_min,zBin_max,nzBin,zz[kk]);

            // x
            for(int ii=0;ii<nx;ii++)
            {
                // y
                for(int jj=0;jj<ny;jj++)
                {
                    // Reset  downward irradiance
                    E_d[ii][jj][kk] = 0;

                    rr_val = sqrt(xx[ii]*xx[ii] + yy[jj]*yy[jj]);
                    th_val = theta_xy(xx[ii],yy[ii]);

                    // Probability of kelp
                    // ZEROED OUT
                    //P_k = P_shade_2d(th_val,rr_val,v_w,theta_w,nLBin,LBin_vals,P_L[depthLayer]);

                    // Fractional radiance loss
                    c_ijk = (1 - P_k) * (a_w + b_w) + P_k * a_k;

                    // theta
                    for(int ll=0;ll<ntheta;ll++)
                    {
                        // phi
                        for(int mm=0;mm<nphi;mm++)
                        {
                            // Spacial finite differences (CD2)
                            // Periodic BC
                            dRdx = (RR[mod(ii+1,nx)][jj][kk][ll][mm]
                                    - RR[mod(ii-1,nx)][jj][kk][ll][mm]) / (2*dx);
                            dRdy = (RR[ii][mod(jj+1,ny)][kk][ll][mm]
                                    - RR[ii][mod(jj-1,ny)][kk][ll][mm]) / (2*dy);
                            if(kk<nz-1)
                                dRdz = (RR[ii][jj][kk+1][ll][mm]
                                        - RR[ii][jj][kk-1][ll][mm]) / (2*dz);
                            // BD2 on last point
                            else
                                dRdz = (3*RR[ii][jj][kk][ll][mm]
                                        - 4*RR[ii][jj][kk-1][ll][mm] + RR[ii][jj][kk-2][ll][mm]);
                            for(int pp=0;pp<ntheta;pp++)
                            {
                                // Gain from scattering
                                R_star = 0;

                                // phi'
                                for(int qq=0;qq<nphi;qq++)
                                {
                                    // Skip qq==mm
                                    if (qq!=mm)
                                    {
                                        // Calculate angle between (theta,phi) and (theta',phi')
                                        anglediff = angleDiff(theta[mm],phi[ll],theta[pp],phi[qq]);
                                        R_star += vsf(anglediff,vsf_th,vsf_vals,nVSF)
                                            * RR[ii][jj][kk][pp][qq] * sin(phi[qq])
                                            * dphi * dtheta;
                                    }
                                    else
                                        continue;
                                }

                                // PDE
                                R_new = (dRdx*sin(phi[mm])*cos(theta[ll])
                                        + dRdy*sin(phi[mm])*sin(theta[ll])
                                        + dRdz*cos(phi[mm]) + R_star) / c_ijk;

                                // Calculate residual
                                residual += abs(RR[ii][jj][kk][ll][mm] - R_new)/(nx*ny*nz*ntheta*nphi);

                                // Update radiance value
                                RR[ii][jj][kk][ll][mm]  = R_new;

                                // Calculate downward irradiance
                                // Integrate radiance over upper hemisphere
                                if(phi[mm] >= PI)
                                {
                                    E_d[ii][jj][kk] -= RR[ii][jj][kk][ll][mm]
                                        * cos(phi[mm]) * sin(phi[mm])
                                        * dtheta * dphi;
                                }

                            }
                        }
                    }
                }
            }
        }

        cout << "iter=" << iter << ": resid = " << residual << endl;
    }
    // Print irradiance array
    // print3d(E_d,nx,ny,nz,"E_d");

    ////////////////////
    // RUN SIMULATION //
    ////////////////////

    /*
    // Loop through time
    for(int nn = 0;nn<nTimeSteps;nn++)
    {
        cout << "t[" << nn << "] = " << tt << endl;
        recalculateLengthDistribution(tt,v_w,theta_w,phi_s,theta_s,nLBin,dLBin,LBin_min,LBin_max,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L,E_d0,attenuationCoef,a_k);
        write2d(P_L,nzBin,nLBin,"P_L",nn+1,outFile);

        // Increment time
        tt += dt;
    }
    */

    //////////////
    // CLEAN UP //
    //////////////

    // Close output file
    outFile.close();

    // Run newly created python script
    system("mkdir -p ../python/kelp_pickle");
    // system("python ../python/kelp_output.py");

    // Delete visualization arrays
    /*
    for(int ii=0;ii<nx;ii++)
    {
        for(int jj=0;jj<ny;jj++)
        {
            delete [] xx_3d[ii][jj];
            delete [] yy_3d[ii][jj];
            delete [] zz_3d[ii][jj];
            delete [] PP_3d[ii][jj];
        }
        delete [] xx_3d[ii];
        delete [] yy_3d[ii];
        delete [] zz_3d[ii];
        delete [] PP_3d[ii];
    }
    delete [] xx_3d;
    delete [] yy_3d;
    delete [] zz_3d;
    delete [] PP_3d;
    */

    for(int ii=0;ii<nzBin;ii++)
        delete P_L[ii];
    delete [] P_L;
    delete [] zBin_vals;
    delete [] LBin_vals;

    // Radiance and irradiance arrays
    delete [] xx,yy,zz;
    for(int ii=0;ii<nx;ii++)
    {
        for(int jj=0;jj<ny;jj++)
        {
            for(int kk=0;kk<nz;kk++)
            {
                for(int ll=0;ll<ntheta;ll++)
                {
                    delete [] RR[ii][jj][kk][ll];
                }
                delete [] RR[ii][jj][kk];
            }
            delete [] RR[ii][jj];
            delete [] E_d[ii][jj];
        }
        delete [] RR[ii];
        delete [] E_d[ii];
    }
    delete [] RR;
    delete [] E_d;

    // VSF
    delete [] vsf_th,vsf_vals;

    cout << "Finish!" << endl;

    return 0;
}


////////////////////
// MATH FUNCTIONS //
////////////////////

// Sign function
double sign(double xx)
{
    if(xx > 0)
        return 1;
    else if(xx < 0)
        return -1;
    else
        return 0;
}

// Minimum of two numbers
double min(double xx, double yy)
{
    if(xx<yy)
        return xx;
    else
        return yy;
}

// Maximum of two numbers
double max(double xx, double yy)
{
    if(xx>yy)
        return xx;
    else
        return yy;
}


// Factorial function
double factorial(int nn)
{
    double result = 1;
    for(int ii=1;ii<=nn;ii++)
        result *= ii;
    return result;
}

// Modulus
// Always return positive remainder, like python
// e.g. mod(-1,4) = 3. By default, C++ returns -1.
int mod(int aa,int bb)
{
    if(aa>=0)
        return aa%bb;
    else
        return bb+aa%bb;
}

// Modified bessel function of the first kind of order 0
// Summation formula from Wolfram Mathworld: http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
// xx: x value
// nn: number of terms
// return: I0(xx)
double I0(double xx,int nn)
{
    double result = 0;

    for(int kk=0;kk<nn;kk++)
        result += pow((xx*xx/4),kk)/(pow(factorial(kk),2));

    return result;
}
// 20 terms seems like a reasonable default
double I0(double xx)
{
    double result = 0;

    // 20 terms seems reasonable
    for(int kk=0;kk<20;kk++)
        result += pow((xx*xx/4),kk)/(pow(factorial(kk),2));

    return result;
}

// von Mises Distribution PDF
// kk: kappa - sharpness parameter
// mu: horizontal shift
// xx: input value
// returns: PDF evaluated at xx
inline double vonMisesPDF(double xx,double kk,double mu)
{
    return exp(kk*cos(xx-mu))/(2*PI*I0(kk));
}

////////////////////////////////
// POLAR COORDINATE FUNCTIONS //
////////////////////////////////

// Calculate polar theta from x,y values
double theta_xy(double xx,double yy)
{
    double theta = sign(yy)*PI/2;
    if(xx!=0)
    {
        theta = atan(yy/xx);
        if(xx<0)
            theta += PI;
    }

    // Shift theta to the interval [-PI,PI]
    if(theta>=PI)
        theta -= 2*PI;

    return theta;
}

// Shift polar coordinate values along a cartesian line *towards* the sun
// theta: initial angular position
// rr: initial radial position
// dz: vertical distance *upward* to shift
// phi_s: zenith angle of the sun
// theta_s: azimuth angle of the sun
// theta_hat: new angular position
// rr_hat: new radial position
// * NOT YET CONSIDERING CHANGE IN PHI_S DUE TO REFRACTION *
void polar_shift(double theta,double rr,double dz,double phi_s,double theta_s,double &theta_hat,double &rr_hat)
{
    // Initial cartesian coordinates
    double xx = rr*cos(theta);
    double yy = rr*sin(theta);

    // Cartesian shifts
    double dx = dz*tan(phi_s)*cos(theta_s);
    double dy = dz*tan(phi_s)*sin(theta_s);

    // New cartesian coordinates
    double xx_hat = xx + dx;
    double yy_hat = yy + dy;

    // New polar coordinates
    theta_hat = theta_xy(xx_hat,yy_hat);
    rr_hat = sqrt(xx_hat*xx_hat + yy_hat*yy_hat);
}

// Calculate angle between two unit vectors in the directions
// (theta,phi) and (theta',phi')
// (azim,zenith)
inline double angleDiff(double theta,double phi,double theta_prime,double phi_prime)
{
    return acos(sin(phi)*sin(phi_prime) * cos(theta-theta_prime) + cos(phi)*cos(phi_prime));
}

///////////////////////
// NUMERICAL METHODS //
///////////////////////

// Composite Simpson's rule for integration (1d)
// ff: integrand
// aa: lower limit of integration
// bb: upper limit of integration
// nn: number of sub-intervals
// return: definite integral of f from a to b with n sub-intervals
double csimp(double (*ff)(double),double aa,double bb,int nn)
{
    double total = 0;
    double hh = (bb-aa)/nn;
    double x_i,x_ip1;

    // Loop through intervals
    for(int ii=0;ii<nn;ii++)
    {
        //x_i
        x_i = aa + hh*ii;
        //x{i+1}
        x_ip1 = aa + hh*(ii+1);
        total += hh/6 * (ff(x_i) + 4*ff((x_i+x_ip1)/2) + ff(x_ip1));
    }

    return total;
}
// Allow for one double parameter, pp, to be passed to ff
double csimp(double (*ff)(double,double),double aa,double bb,int nn,double pp)
{
    double total = 0;
    double hh = (bb-aa)/nn;
    double x_i,x_ip1;

    // Loop through intervals
    for(int ii=0;ii<nn;ii++)
    {
        //x_i
        x_i = aa + hh*ii;
        //x{i+1}
        x_ip1 = aa + hh*(ii+1);
        total += hh/6 * (ff(x_i,pp) + 4*ff((x_i+x_ip1)/2,pp) + ff(x_ip1,pp));
    }

    return total;
}
// Allow for two double parameters, pp1 and pp2, to be passed to ff
double csimp(double (*ff)(double,double,double),double aa,double bb,int nn,double pp1,double pp2)
{
    double total = 0;
    double hh = (bb-aa)/nn;
    double x_i,x_ip1;

    // Loop through intervals
    for(int ii=0;ii<nn;ii++)
    {
        //x_i
        x_i = aa + hh*ii;
        //x{i+1}
        x_ip1 = aa + hh*(ii+1);
        total += hh/6 * (ff(x_i,pp1,pp2) + 4*ff((x_i+x_ip1)/2,pp1,pp2) + ff(x_ip1,pp1,pp2));
    }

    return total;
}
// Allow for two double pointer parameters, pp1 and pp2 and one int parameter, pp3, to be passed to ff
double csimp(double (*ff)(double,double*,double*,int),double aa,double bb,int nn,double* pp1,double* pp2,int pp3)
{
    double total = 0;
    double hh = (bb-aa)/nn;
    double x_i,x_ip1;

    // Loop through intervals
    for(int ii=0;ii<nn;ii++)
    {
        //x_i
        x_i = aa + hh*ii;
        //x{i+1}
        x_ip1 = aa + hh*(ii+1);
        total += hh/6 * (ff(x_i,pp1,pp2,pp3) + 4*ff((x_i+x_ip1)/2,pp1,pp2,pp3) + ff(x_ip1,pp1,pp2,pp3));
    }

    return total;
}


///////////////////
// INTERPOLATION //
///////////////////

// Linear interpolation
// xx: value to interpolate
// x0: low x value
// x1: high x value
// y0: y value corresponding to x0
// y1: y value corresponding to x1
// return: interpolated y value
double linear(double xx,double x0,double x1,double y0,double y1)
{
    return y0 + (xx-x0)/(x1-x0) * (y1 - y0);
}

// Bilinear interpolation
// xx: 1d ordered array of known x values
// yy: 1d ordered array of known y values
// zz: 2d array of known z values
// nn: length of xx
// mm: length of yy
// xval: x coordinate of interpolation point
// yval: y coordinate of interpolation point
// return: interpolated z value
double bilinear(double* xx,double* yy,double** zz,int nn,int mm,double xval,double yval)
{
    // Find closest known points
    int ii = findBin(xx[0],xx[nn-1],nn,xval);
    int jj = findBin(yy[0],yy[mm-1],mm,yval);

    double result;

    // Upper left
    double A = (xval-xx[ii])*(yy[jj+1]-yval);
    // Upper right
    double B = (xx[ii+1]-xval)*(yy[jj+1]-yval);
    // Lower left
    double C = (xval-xx[ii])*(yval-yy[jj]);
    // Lower right
    double D = (xx[ii+1]-xval)*(yval-yy[jj]);

    // Not on upper edges
    if(ii < nn-1 && jj < mm-1)
        result = (A*zz[ii+1][jj] + B*zz[ii][jj] + C*zz[ii+1][jj+1] + D*zz[ii][jj+1])/(A+B+C+D);

    // Right edge
    else if (ii == nn-1 && jj < mm-1)
    {
        A = yy[jj+1]-yval;
        C = yval-yy[jj];
        result = (A*zz[nn-1][jj] + C*zz[nn-1][jj+1])/(A+C);
    }

    // Upper edge
    else if (ii < nn-1 && jj == mm-1)
    {
        B = xx[ii+1]-xval;
        D = xx[ii+1]-xval;
        result = (B*zz[ii][mm-1] + D*zz[ii+1][mm-1])/(B+D);
    }

    // Upper left corner
    else if (ii == nn-1 && jj == mm-1)
        result = zz[nn-1][mm-1];

    // Out of bounds
    else
    {
        cout << "Out of bounds!" << endl;
        result = 0;
    }
    return result;
}

// Evaluate volume scattering function by interpolating given points
// STRANGE THING!!!
// It seems that this interpolation is very close, but for some reason off by ~1e-2.
// The interpolated values seem more accurate near the given points, but bulge between them
// Similar to the shape of x(x-1) between 0 and 1
double vsf(double xx,double* vsf_th,double* vsf_vals,int nVSF)
{
    int ind = binarySearch(xx,vsf_th,0,nVSF-1);
    if(0 <= ind && ind < nVSF)
        return linear(xx,vsf_th[ind],vsf_th[ind+1],vsf_vals[ind],vsf_vals[ind+1]);
    else if(ind >= nVSF)
        return vsf_vals[nVSF-1];
    else
        return vsf_vals[0];
}
double vsf1(double xx,double* vsf_th,double* vsf_vals,int nVSF)
{
    return vsf(xx,vsf_th,vsf_vals,nVSF)*sin(xx)*2*PI;
}

///////////////////////
// UTILITY FUNCTIONS //
///////////////////////

// Count lines in a file
//Count the number of timesteps
int countLines(ifstream &inFile)
{
    string line;
    int numLines=0;
    int lineNum=0;

    //Count number of timesteps
    while(getline(inFile,line))
        numLines++;

    //Unset eof flag (if set) & return to beginning of file
    inFile.clear();
    inFile.seekg(0,ios::beg);

    return numLines;
}


// Binary search
// Search through an ordered list by splitting it in half until the closest lower value is found
// xx: value to search for
// vv: double pointer to search through
// nn: length of vv
// return: index of bin that xx belongs to
int binarySearch(double xx,double* vv,int first,int last)
{
    int middle = (first+last)/2; // Integer division
    // cout << "x=" << xx << endl;
    // cout << "f: " << first << ", m: " << middle << ", l: " << last << endl;
    // cout << "v[m]="  << vv[middle] << endl;

    if(last - first > 1)
    {
        // cout << "l-f = " << last - first << endl;
        if(xx < vv[middle])
        {
            // cout << "<" << endl;
            return binarySearch(xx,vv,first,middle);
        }
        else if(xx > vv[middle])
        {
            // cout << ">" << endl;
            return binarySearch(xx,vv,middle,last);
        }
        else
        {
            // cout << "=" << endl;
            return middle;
        }
    }
    else
    {
        // cout << "first=" << first << endl;
        return first;
    }
}

// Bin finder
// Given an interval between aa and bb divided into nn bins,
// find the number of the bin which xx belongs to.
// aa: min value
// bb: max value
// nn: number of bins
// xx: given value
// return: index of bin that xx belongs to
inline int findBin(double aa, double bb,int nn,double xx)
{
    return floor((xx-aa)/(bb-aa)*nn);
}

// Index of the next bin if xx does not fall
// exactly on a bin edge, or the current bin if it does
// Remember: bins are upper half open.
// **This is not necessarily the bin that xx belongs to**
// aa: min value
// bb: max value
// nn: number of bins
// xx: given value
// return: index of bin
inline int upperBin(double aa, double bb,int nn,double xx)
{
    return ceil((xx-aa)/(bb-aa)*nn);
}

// Read two column double data from file
// inFile: input stream to read from
// col1,col2: allocated doubles whose length equals the number of rows in the file
void read2column(ifstream &inFile,double* col1,double* col2,int numLines)
{
    string line;
    for(int ii=0;ii<numLines;ii++)
        inFile >> col1[ii] >> col2[ii];
}

// Print 3d array for input to python
void print3d(double*** v, int nx, int ny,int nz, const char *name)
{
    cout << name << " = array([";
    for(int ii=0;ii<nx;ii++)
    {
        cout << "[";
        for(int jj=0;jj<ny;jj++)
        {
            cout << "[";
            for(int kk=0;kk<nz;kk++)
            {
                cout << v[ii][jj][kk];
                if(kk<nz-1)
                    cout << ",";
            }
            cout << "]";
            if(jj<ny-1)
                cout << "," << endl;
        }
        cout << "]";
        if(ii<nx-1)
            cout << "," << endl << endl;
    }
    cout << "])" << endl;
}

// Write 3d array for input to python to file
void write3d(double*** v, int nx, int ny,int nz, const char *name,ofstream &outFile)
{
    outFile << name << " = array([";
    for(int ii=0;ii<nx;ii++)
    {
        outFile << "[";
        for(int jj=0;jj<ny;jj++)
        {
            outFile << "[";
            for(int kk=0;kk<nz;kk++)
            {
                outFile << v[ii][jj][kk];
                if(kk<nz-1)
                    outFile << ",";
            }
            outFile << "]";
            if(jj<ny-1)
                outFile << "," << endl;
        }
        outFile << "]";
        if(ii<nx-1)
            outFile << "," << endl << endl;
    }
    outFile << "])" << endl;

    // Save with pickle
    outFile << "with open('../python/kelp_pickle/" << name << ".pickle','wb') as pickle_file:" << endl;
    outFile << "    pickle.dump(" << name << ",pickle_file)" << endl;
}

// Print 2d array for input to python
void print2d(double** v, int nx, int ny, const char *name)
{
    cout << name << " = array([" << endl;
    for(int ii=0;ii<nx-1;ii++)
    {
        cout << "[";
        for(int jj=0;jj<ny-1;jj++)
            cout << v[ii][jj] << ",";
        cout << v[ii][ny-1] << "]," << endl;
    }
    cout << "[";
    for(int jj=0;jj<ny-1;jj++)
        cout << v[nx-1][jj] << ",";
    cout << v[nx-1][ny-1] << "]])" << endl;
}

// Write 2d array for input to python
void write2d(double** v, int nx, int ny, const char *name,ofstream& outFile)
{
    outFile << name << " = array([" << endl;
    for(int ii=0;ii<nx-1;ii++)
    {
        outFile << "[";
        for(int jj=0;jj<ny-1;jj++)
            outFile << v[ii][jj] << ",";
        outFile << v[ii][ny-1] << "]," << endl;
    }
    outFile << "[";
    for(int jj=0;jj<ny-1;jj++)
        outFile << v[nx-1][jj] << ",";
    outFile << v[nx-1][ny-1] << "]])" << endl;
}

// Write 2d array for input to python with index
void write2d(double** v, int nx, int ny, const char *name,int index,ofstream& outFile)
{
    outFile << name << "[" << index << "] = array([" << endl;
    for(int ii=0;ii<nx-1;ii++)
    {
        outFile << "[";
        for(int jj=0;jj<ny-1;jj++)
            outFile << v[ii][jj] << ",";
        outFile << v[ii][ny-1] << "]," << endl;
    }
    outFile << "[";
    for(int jj=0;jj<ny-1;jj++)
        outFile << v[nx-1][jj] << ",";
    outFile << v[nx-1][ny-1] << "]])" << endl;
}

// Print 1d array for input to python
void print1d(double* v, int nn, const char *name)
{
    cout << name << " = array([";
    for(int ii=0;ii<nn-1;ii++)
        cout << v[ii] << ",";
    cout << v[nn-1] << "])" << endl;
}

// Create polar filled contour plot of values in Python
void contourf(double* theta,double* rr,double** zz,int ntheta,int nr)
{
    print1d(theta,ntheta,"th");
    print1d(rr,nr,"rr");
    print2d(zz,ntheta,nr,"zz");
    cout << "rr,th = meshgrid(rr,th)" << endl;
    cout << "clf()" << endl;
    cout << "gca(projection='polar')" << endl;
    cout << "contourf(th,rr,zz,cmap='YlGn')" << endl;
    cout << "colorbar()" << endl;
}


//////////////////////////////
// MODEL-SPECIFIC FUNCTIONS //
//////////////////////////////

// PDF of frond angle distribution
// v_w: water current speed
// theta_w: angle of water current in radians
// xx: input value
// return: probability density at xx
inline double P_theta_f(double xx,double v_w,double theta_w)
{
    double result = 1/(2*PI);
    if(v_w != 0)
        result = vonMisesPDF(xx,v_w,theta_w);
    return result;
}

// Minimum LL (frond length) for shading a given point
// as a function of theta_f (frond angle)
inline double L_min_shade(double theta_p,double r_p,double theta_f)
{
    double theta_prime = theta_p - theta_f + PI/2;
    double S = sign(PI/2-theta_prime);
    return r_p*(sin(theta_prime)+2*S*FR/(1+FS)*cos(theta_prime));
}

// Minimum theta_f (frond angle ) for shading a given point
// as a function of LL (frond length)
// theta_max_shade = 2*theta_p - theta_min_shade
inline double theta_min_shade(double theta_p,double r_p,double LL)
{
    double qq = 2*FR/(1+FS);
    // Unshifted value, assuming point is at theta_p = PI/2
    double unshifted = min( atan(1/qq) + acos(LL/r_p/sqrt(1+qq*qq)) , PI/2 );
    // Rotate point to theta_p
    return unshifted + theta_p - PI/2;
}

// Expected number of fronds shading a point in 3d space
// considering only incident light (not diffuse)
// Integrate P_shade_2d over z, shifting theta_p and r_p towards the sun appropriately
double N_shade_3d(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L)
{
    // Probability of shading
    double PP = 0;

    // Shifted polar coordinates
    double theta_hat,r_hat;

    // Current depth
    double z_f;

    // Index of z bin containing point
    int k_p = findBin(zBin_min,zBin_max,nzBin,z_p);

    // Loop over depths from surface to z_p, inclusive
    for(int kk=0;kk<k_p;kk++)
    {
        z_f = zBin_vals[kk];
        polar_shift(theta_p,r_p,z_p-z_f,phi_s,theta_s,theta_hat,r_hat);
        PP += RHO/dzBin * P_shade_2d(theta_hat,r_hat,v_w,theta_w,nLBin,LBin_vals,P_L[kk]);
    }

    /*
    print1d(z_f,nz,"z_f");
    print1d(PP_z,nz,"PP_z");
    cout << "plot(z_f,PP_z)" << endl;
    */

    return PP;
}

double P_shade_2d(double theta_p,double r_p,double v_w,double theta_w,int nLBin,double* LBin_vals,double* P_L)
{
    // Recalculate information about dLBin_vals
    double dLBin = LBin_vals[1] - LBin_vals[0];
    double LBin_min = LBin_vals[0];
    double LBin_max = LBin_vals[nLBin-1] + dLBin;

    // Whether integration has to stop in theta-variable region
    bool early_cutoff = false;

    // Calculate LL limits of numerical integration
    double L_star = L_min_shade(theta_p,r_p,theta_p+ALPHA);
    int k0 = findBin(LBin_min,LBin_max,nLBin,r_p);
    // Don't let k_star exceed number of bins
    int k_star = upperBin(LBin_min,LBin_max,nLBin,L_star);

    // Check for early cutoff
    if(nLBin <= k_star)
    {
        early_cutoff = true;
        k_star = nLBin;
    }

    // theta_f limits of integration for this length & next & avg
    double theta_min_k;
    double theta_min_kp1;
    double theta_min,theta_max;

    // Number of sub-intervals to use for integrating P_theta_f
    int n_csimp = 20;

    // Probability of shading in theta-variable, theta-constant regions and total
    double P_var = 0;
    double P_const = 0;
    double P_total = 0;

    // Integrate very first bin separately
    // Don't let theta_min exceed theta_p
    // Left endpoint
    theta_min_k = theta_min_shade(theta_p,r_p,LBin_vals[k0]);
    // Right endpoint
    theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[k0+1]);
    // Use average of left and right endpoints for integration limits
    theta_min = (theta_min_k + theta_min_kp1)/2;
    theta_max = 2*theta_p - theta_min;

    // Accumulate integral
    P_var += csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w) * P_L[k0];

    // Integrate over remainder of theta-variable region
    for(int kk=k0+1;kk<k_star-1;kk++)
    {
        // Left endpoint
        theta_min_k = theta_min_kp1;

        // Right endpoint
        theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[kk+1]);

        // Use average of left and right endpoints for integration limits
        theta_min = (theta_min_k + theta_min_kp1)/2;
        theta_max = 2*theta_p - theta_min;

        // Accumulate integral
        P_var += csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w) * P_L[kk];
    }

    // Integrate over last theta-variable bin separately to efficiently check for early cutoff
    // Left endpoint
    theta_min_k = theta_min_kp1;

    // Right endpoint
    if(early_cutoff)
        theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[k_star-1]+dLBin);
    else
        theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[k_star]);

    // Use average of left and right endpoints for integration limits
    theta_min = (theta_min_k + theta_min_kp1)/2;
    theta_max = 2*theta_p - theta_min;

    // Accumulate integral
    P_var += csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w) * P_L[k_star-1];

    // Sum lengths
    theta_min = theta_p - ALPHA;
    theta_max = theta_p + ALPHA;
    for(int kk=k_star;kk<nLBin;kk++)
        P_const += P_L[kk];

    // Multiply by integral over theta_f
    P_const *= csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w);

    // Calculate total integral
    P_total = P_var + P_const;

    return P_total;
}

// Calculate available light at a given 3d point
// Not considering ambient light
double availableLight(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double E_d0,double attenuationCoef,double a_k)
{
    double nShade = N_shade_3d(theta_p,r_p,z_p,v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L);
    double absorptionFactor = pow((1-a_k),nShade);
    double attenuationFactor = exp(-attenuationCoef*z_p);

    return E_d0 * absorptionFactor * attenuationFactor;
}

// Transform a point (ss,tt) in the unit square: [-1,1]x[-1,1] to the corresponding point (xx,yy) on a frond of length L at an angle theta_f
inline void frondTransform(double ss,double tt,double &xx,double &yy,double theta_f,double LL)
{
    xx = LL*(-FR*(FS*(ss - 1)*(tt + 1) - (ss + 1)*(2*FS + tt + 1))*cos(theta_f)
            + (FS + 1)*(ss - tt)*sin(theta_f))/(4*FR*(FS + 1));
    yy = -LL*(FR*(FS*(ss - 1)*(tt + 1) - (ss + 1)*(2*FS + tt + 1))*sin(theta_f)
            + (FS + 1)*(ss - tt)*cos(theta_f))/(4*FR*(FS + 1));
}

// Jacobian of frond transform
inline double JJ(double ss,double tt,double LL)
{
    return LL*LL*(-FS*FS*ss - FS*FS*tt + 2*FS*FS + 4*FS + ss + tt + 2)/(16*FR*(FS + 1)*(FS + 1));
}

// Calculate the light absorbed by a frond with a particular angle (theta_f) and length (LL)
// at a particular depth (z_f)
// Use n=2 Gaussian quadrature product rule to integrate light field over frond area by transforming frond to unit square
double calculateLightAbsorption(double theta_f,double LL,double z_f,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double E_d0,double attenuationCoef,double a_k)
{
    // Abscissas on unit square
    static double sa[2] = {-1/sqrt(3),1/sqrt(3)};

    // Weights (not necessary for n=2)
    // static double ww[2] = {1,1};

    // Coordinates of frond abscissae
    double xfa,yfa,tfa,rfa;

    // Total light absorbed by frond
    double lightAbsorbed = 0;

    // Part of integral due to a specific abscissa
    double absc_val;

    // Loop through 2x2 abscissas
    for(int ii=0;ii<2;ii++)
    {
        for(int jj=0;jj<2;jj++)
        {
            // Calculate abscissae
            frondTransform(sa[ii],sa[jj],xfa,yfa,theta_f,LL);

            // Convert to polar coordinates
            tfa = theta_xy(xfa,yfa);
            rfa = sqrt(xfa*xfa + yfa*yfa);

            // Calculate contribution of this abscissa
            absc_val = availableLight(tfa,rfa,z_f,v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L,E_d0,attenuationCoef,a_k);

            // Multiply by weights according to quadrature product rule
            // absc_val *= ww[ii] * ww[jj];

            // Multiply by the Jacobian and add to total integral
            lightAbsorbed += absc_val * abs(JJ(sa[ii],sa[jj],LL));
        }
    }

    //cout << "LA: (" << theta_f << "," << LL << "," << z_f << ") = " << lightAbsorbed << endl;
    return lightAbsorbed;
}

// New frond length from growth as a function of light absorbed, frond length, and time
inline double newLength(double lightAbsorbed,double LL,double tt)
{
    double light_param = 1;
    double length_param = 1;
    double time_param = 1;
    double growth = lightAbsorbed * light_param / ((LL*length_param+1) * (tt*time_param+1))/10;
    return LL + growth;
}

// Recalculate length distribution by integrating P_theta_f*P_L over R for each L_i,
// where R is the region such that newLength(theta_f,LL) \in [ L_i , L_{i+1} )
// This function recalculates for all z bins, then updates
void recalculateLengthDistribution(double tt,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double dLBin,double LBin_min,double LBin_max,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double E_d0,double attenuationCoef,double a_k)
{
    // Number of points to sample from theta_f
    int n_theta_f = 5;
    double dtheta_f = 2*PI/n_theta_f;

    // Current sample coordinates
    double theta_f;
    double LL_current;
    double z_f;

    // Light absorbed by frond
    double lightAbsorbed;

    // theta_f coordinates and LL indices of points in R for a given LL bin
    // and new length distribution
    int** nPoints = new int*[nzBin];
    vector<double>** R_theta_f = new vector<double>*[nzBin];
    vector<int>** R_LL = new vector<int>*[nzBin];
    double** P_L_new = new double*[nzBin];
    for(int kk=0;kk<nzBin;kk++)
    {
        nPoints[kk] = new int[nLBin];
        R_theta_f[kk] = new vector<double>[nLBin];
        R_LL[kk] = new vector<int>[nLBin];
        P_L_new[kk] = new double[nLBin];
    }


    // Number of points which fall into each LL bin

    // New length coordinate
    double LL_new;

    // Bin of new length coordinate
    int new_bin;

    // Loop through depth bins
    #pragma omp parallel for private(theta_f,z_f,LL_current,lightAbsorbed,LL_new,new_bin)
    for(int kk=0;kk<nzBin;kk++)
    {
        // Reset point counters and P_L_new values
        for(int jj=0;jj<nLBin;jj++)
        {
            nPoints[kk][jj] = 0;
            P_L_new[kk][jj] = 0;
        }

        // Reset theta_f sample point
        theta_f = -PI;

        // Get z_f value for this index (use bin midpoint)
        z_f = zBin_vals[kk] + dzBin/2;

        // Loop through theta_f-LL space of possible living fronds (LL>0)
        for(int ii=0;ii<n_theta_f;ii++)
        {
            for(int jj=1;jj<nLBin;jj++)
            {
                // Get LL value for this index (use bin midpoint)
                LL_current = LBin_vals[jj] + dLBin/2;

                // Calculate light absorbed by a frond with this particular theta_f and LL
                lightAbsorbed = calculateLightAbsorption(theta_f,LL_current,z_f,v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L,E_d0,attenuationCoef,a_k);

                // Calculate new frond size for this type of frond
                LL_new = newLength(lightAbsorbed,LL_current,tt);

                // Determine which bin this new length falls into
                new_bin = findBin(LBin_min,LBin_max,nLBin,LL_new);
                // but don't let it exceed the number of bins!
                new_bin = min(new_bin,nLBin-1);
                // nor be less than 0
                new_bin = max(new_bin,0);

                // Save the this theta_f value in the appropriate vector
                R_theta_f[kk][new_bin].push_back(theta_f);
                // Save the index of this LL value in the appropriate vector
                R_LL[kk][new_bin].push_back(jj);

                // Increment point counter for this bin
                nPoints[kk][new_bin]++;
            }
            // Dead fronds stay dead (LL<0)
            R_theta_f[kk][0].push_back(theta_f);
            R_LL[kk][0].push_back(0);
            nPoints[kk][0]++;

            // Increment sample point
            theta_f += dtheta_f;
        }

        // Integrate 2d population distribution over appropriate region for each new length bin
        // Loop through new length bins
        for(int jj=0;jj<nLBin;jj++)
        {
            // Loop over points (use 2d midpoint rule)
            for(int nn=0;nn<nPoints[kk][jj];nn++)
            {
                if(P_L[kk][R_LL[kk][jj][nn]]<0)
                {
                    cout << "PL (" << nn << "," << jj << "," << kk << ";" << R_LL[kk][jj][nn] << ") = " << P_L[kk][R_LL[kk][jj][nn]] << endl;
                }
                P_L_new[kk][jj] += P_theta_f(R_theta_f[kk][jj][nn],v_w,theta_w)*dtheta_f*P_L[kk][R_LL[kk][jj][nn]];
            }
        }
    }

    // Update length distribution
    for(int kk=0;kk<nzBin;kk++)
        for(int jj=0;jj<nLBin;jj++)
            P_L[kk][jj] = P_L_new[kk][jj];

    // Deallocate arrays
    for(int kk=0;kk<nzBin;kk++)
    {
        delete [] nPoints[kk];
        delete [] P_L_new[kk];
        delete [] R_theta_f[kk];
        delete [] R_LL[kk];
    }

    delete [] nPoints,P_L_new,R_theta_f,R_LL;
}


/////////
// FIN //
/////////

