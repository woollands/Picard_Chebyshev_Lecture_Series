/*
   Copyright (c) TAMU LASR Lab, Austin Probe 2014
   Texas A&M University - AEROSPACE department
   File name     : EGM2008.c
   Description   : Calculates gravity according to the EGM 2008 Spherical 
		   Harmonic Model
*/

#include "./EGM2008.h"

void matmul(double* A, double* B, double* OUT, int m, int n, int q)
{
	for(int i=0;i<m; i++){
		for(int j=0;j<q;j++){
			double sum = {0.0};
			for(int j1=0;j1<n;j1++)
				sum += A[IDX2F(i+1,j1+1,m)]*B[IDX2F(j1+1,j+1,n)];
			OUT[IDX2F(i+1,j+1,m)] = sum;
        }
   	}
}

void EGM2008( double* p, double* Gxyz, int DEG)
{
	// determine radius of this thread's node
	double r             = {0.0};
	double phic          = {0.0};
	double lambda        = {0.0};
	double slambda       = {0.0};
	double clambda       = {0.0};
    double smlambda[Max_Degree+1] = {0.0};
    double cmlambda[Max_Degree+1] = {0.0};
    
    double x = p[0];
    double y = p[1];
    double z = p[2];
    
    // Compute geocentric radius
    r = pow( x*x + y*y + z*z , 0.5 );
    // Compute geocentric latitude
    phic  = asin( z / r );
    // Compute lambda
    lambda  = atan2( y, x );
    while (lambda<0)
        lambda = lambda+2*Pi;
    while (lambda>=2*Pi)
        lambda = lambda-2*Pi;
    
    
    slambda = sin(lambda);
    clambda = cos(lambda);
    smlambda[0] = 0.0;
    cmlambda[0] = 1.0;
    smlambda[1] = slambda;
    cmlambda[1] = clambda;
    
    for(int m=2;m<DEG+1;m++){
        smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
        cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
    }
    // Compute normalized associated legendre polynomials
    double P [(Max_Degree+3)*(Max_Degree+3)] = {0.0};
    double scaleFactor [(Max_Degree+3)*(Max_Degree+3)] = {0.0};
    
    loc_gravLegendre( phic, scaleFactor, P , DEG);
    
    loc_gravityPCPF( p, P, DEG, smlambda, cmlambda, r, scaleFactor, Gxyz );
    
    
}


void loc_gravLegendre( double phi, double* scaleFactor, double* P, int DEG )
{
	// loc_GRAVLEGENDRE internal function computing normalized associated
    // legendre polynomials, P, via recursion relations for spherical harmonic gravity
	int k, p;
	double cphi = {0.0};
	double sphi = {0.0};
    
    
    cphi = cos(0.5*Pi - phi);
    sphi = sin(0.5*Pi - phi);
    // Seeds for recursion formula
    P[IDX2F(1,1,Max_Degree+3)] = 1.0;            // n = 0, m = 0;
    scaleFactor[IDX2F(1,1, Max_Degree+3)] = 0.0;
    P[IDX2F(2,1, Max_Degree+3)] = sqrt(3.0)*cphi ; // n = 1, m = 0;
    scaleFactor[IDX2F(2,1,Max_Degree+3)]  = 1.0;
    P[IDX2F(2,2,Max_Degree+3)] = sqrt(3.0)*sphi; // n = 1, m = 1;
    scaleFactor[IDX2F(2,2,Max_Degree+3)] = 0.0;
    
    
    for (int nn = 2; nn <= DEG+2;nn++){
        double n = (double)nn;
        k = nn + 1;
        for(int mm=0; mm<=n;mm++) {
            double m = (double)mm;
            p = mm + 1;
            // Compute normalized associated legendre polynomials, P, via recursion relations
            // Scale Factor needed for normalization of dUdphi partial derivative
            if (n == m){
                P[IDX2F(k,k,Max_Degree+3)] = sqrt(2*n+1.0)/sqrt(2.0*n)*sphi*P[IDX2F(k-1,k-1, Max_Degree+3)];
                scaleFactor[IDX2F(k,k,Max_Degree+3)] = 0.0;
            }
            else if (m == 0){
                P[IDX2F(k,p,Max_Degree+3)] = (sqrt(2*n+1)/n)*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - (n-1)/sqrt(2*n-3)* P[IDX2F(k-2,p, Max_Degree+3)] );
                scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+1)*(n)/2);
            }
            else {
                P[IDX2F(k,p,Max_Degree+3)] = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1.0)*cphi*P[IDX2F(k-1,p,Max_Degree+3)] - sqrt(n+m-1.0)*sqrt(n-m-1)/sqrt(2*n-3)*P[IDX2F(k-2,p,Max_Degree+3)] );
                scaleFactor[IDX2F(k,p,Max_Degree+3)] = sqrt( (n+m+1)*(n-m));
            }
        }
    }
}


void loc_gravityPCPF( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor, double* Gxyz )
{
	// loc_GRAVITYPCPF internal function computing gravity in planet-centered
    // planet-fixed (PCEF) coordinates using PCPF position, desired
    // degree/order, normalized associated legendre polynomials, normalized
    // spherical harmonic coefficients, trigonometric functions of geocentric
    // latitude and longitude, planetary constants, and radius to center of
    // planet. Units are MKS.
    
	int k, j;
	double radu;
	double rRatio  = {0.0};
	double rRatio_n  = {0.0};
    
	// initialize summation of gravity in radial coordinates
	double dUdrSumN   = {1.0};
	double dUdphiSumN    = {0.0};
	double dUdlambdaSumN = {0.0};
    
	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};
    
	double dUdr    = {0.0};
	double dUdphi  = {0.0};
	double dUdlambda = {0.0};
    
    
    
    double x = p[0];
    double y = p[1];
    double z = p[2];
    radu = r;
    rRatio = REQ/radu;
    rRatio_n = rRatio;
    // summation of gravity in radial coordinates
    for (int n = 2; n <= DEG; n++) {
        k = n+1;
        rRatio_n = rRatio_n*rRatio;
        dUdrSumM      = 0.0;
        dUdphiSumM    = 0.0;
        dUdlambdaSumM = 0.0;
        for (int m = 0; m <= n; m++){
            j = m+1;
            dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
        }
        dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
        dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
        dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;
    }
    
    // gravity in spherical coordinates
    dUdr      = -MU/(radu*radu)*dUdrSumN ;
    dUdphi    =  MU/radu*dUdphiSumN ;
    dUdlambda =  MU/radu*dUdlambdaSumN ;
    
    //gravity in ECEF coordinates
    Gxyz[0] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y;
    Gxyz[1] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x;
    Gxyz[2] = (1.0/radu)*dUdr*z + ((sqrt(x*x + y*y))/(radu*radu))*dUdphi;
	
    // special case for poles
    /*
     atPole = abs(atan2(p[IDX2F(i+1,3,M)],sqrt(p[IDX2F(i+1,1,M)]^2 + p[IDX2F(i+1,2,M)]^2)))==pi/2;
     if any(atPole){
     gx(atPole) = 0;
     gy(atPole) = 0;
     gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
     }
     */ 

}
