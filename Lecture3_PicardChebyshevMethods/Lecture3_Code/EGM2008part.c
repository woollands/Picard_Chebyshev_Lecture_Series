/*
   Copyright (c) TAMU LASR Lab, Austin Probe 2014
   Texas A&M University - AEROSPACE department
   File name     : EGM2008.c
   Description   : Calculates gravity according to the EGM 2008 Spherical
		   Harmonic Model
*/

#include "./EGM2008part.h"

//Round a / b to nearest higher integer value
int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

void matmul(double* A, double* B, double* OUT, int m, int n, int q)
{
    int i;
    int j;
    int j1;
	for(i=0;i<m; i++){
		for(j=0;j<q;j++){
			double sum = {0.0};
			for(j1=0;j1<n;j1++)
				sum += A[IDX2F(i+1,j1+1,m)]*B[IDX2F(j1+1,j+1,n)];
			OUT[IDX2F(i+1,j+1,m)] = sum;
        }
   	}
}

void EGM2008( double* p, double* Gxyz, double* delV_del_r_phi_lambda, int DEG)
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
    int m;
    double P [(Max_Degree+3)*(Max_Degree+3)] = {0.0};
    double scaleFactor [(Max_Degree+3)*(Max_Degree+3)] = {0.0};

    // Compute geocentric radius
    r = pow( x*x + y*y + z*z , 0.5 );
    // Compute geocentric latitude
    phic  = asin( z / r );
    // Compute lambda
    lambda  = atan2( y, x );
    while (lambda<0){
        lambda = lambda+2*Pi;
    }
    while (lambda>=2*Pi){
        lambda = lambda-2*Pi;
    }

    slambda = sin(lambda);
    clambda = cos(lambda);
    smlambda[0] = 0.0;
    cmlambda[0] = 1.0;
    smlambda[1] = slambda;
    cmlambda[1] = clambda;


    for(m=2;m<DEG+1;m++){
        smlambda[m] = 2.0*clambda*smlambda[m-1] - smlambda[m-2];
        cmlambda[m] = 2.0*clambda*cmlambda[m-1] - cmlambda[m-2];
    }
    // Compute normalized associated legendre polynomials

    loc_gravLegendre( phic, scaleFactor, P , DEG);

    loc_gravityPCPF( p, P, DEG, smlambda, cmlambda, r, scaleFactor, Gxyz, delV_del_r_phi_lambda );


}


void loc_gravLegendre( double phi, double* scaleFactor, double* P, int DEG )
{
	// loc_GRAVLEGENDRE internal function computing normalized associated
    // legendre polynomials, P, via recursion relations for spherical harmonic gravity
	int k, p, nn, mm;
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

    for ( nn = 2; nn <= DEG+2;nn++){
        double n = (double)nn;
        k = nn + 1;
        for(mm=0; mm<=n;mm++) {
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


void loc_gravityPCPF( double* p, double* P, int DEG, double* smlambda, double* cmlambda, double r, double* scaleFactor, double* Gxyz, double* delV_del_r_phi_lambda )
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

	// initialize summation of gravity partials in radial coordinates
	double d2Udr2SumN          = {2.0};     // del^2(U)/del(r)^2
	double d2Udphi2SumN        = {0.0};     // del^2(U)/del(phi)^2
	double d2Udlambda2SumN     = {0.0};     // del^2(U)/del(lambda)^2

	double d2UdrdphiSumN       = {0.0};     // del^2(U)/(del(r)del(phi))
	double d2UdrdlambdaSumN    = {0.0};     // del^2(U)/(del(r)del(lambda))
	double d2UdphidlambdaSumN  = {0.0};     // del^2(U)/(del(phi)del(lambda))

	double d2UdphidrSumN       = {0.0};     // del^2(U)/(del(phi)del(r))
	double d2UdlambdadrSumN    = {0.0};     // del^2(U)/(del(lambda)del(r))
	double d2UdlambdadphiSumN  = {0.0};     // del^2(U)/(del(lambda)del(phi))

	double dUdrSumM   = {0.0};
	double dUdphiSumM = {0.0};
	double dUdlambdaSumM = {0.0};

	double d2Udr2SumM       = {0.0};
	double d2Udphi2SumM     = {0.0};
	double d2Udlambda2SumM  = {0.0};

	double d2UdrdphiSumM      = {0.0};
	double d2UdrdlambdaSumM   = {0.0};
	double d2UdphidlambdaSumM = {0.0};

	double d2UdphidrSumM      = {0.0};
	double d2UdlambdadrSumM   = {0.0};
	double d2UdlambdadphiSumM = {0.0};

	double dUdr           = {0.0};
	double dUdphi         = {0.0};
	double dUdlambda      = {0.0};

	double d2Udr2         = {0.0};
	double d2Udphi2       = {0.0};
	double d2Udlambda2    = {0.0};

	double d2Udrdphi      = {0.0};
	double d2Udrdlambda   = {0.0};
	double d2Udphidlambda = {0.0};

	double d2Udphidr      = {0.0};
	double d2Udlambdadr   = {0.0};
	double d2Udlambdadphi = {0.0};

    double x = p[0];
    double y = p[1];
    double z = p[2];
    int n;
    int m;
    
    // Initialize partials of Associated Legendre Functions
    double dPdphi_nm = {0.0};
    double dPdphi_nmp1 = {0.0};
    double d2Pdphi2_nm = {0.0};
    
    double d2Udrdphi_diff, d2Udrdlambda_diff, d2Udphidlambda_diff;

	// Commonly used terms
	double tanPhi, secPhi;
	tanPhi = z/(sqrt(x*x + y*y));
	secPhi = r/(sqrt(x*x + y*y));
//	tanPhi = p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2));
//	secPhi = r./(sqrt(p(:,1).^2 + p(:,2).^2));
    radu = r;
    rRatio = REQ/radu;
    rRatio_n = rRatio;



    // summation of gravity in radial coordinates
    for (n = 2; n <= DEG; n++) {
        k = n+1;
        rRatio_n = rRatio_n*rRatio;

        dUdrSumM      = 0.0;
        dUdphiSumM    = 0.0;
        dUdlambdaSumM = 0.0;

		d2Udr2SumM         = 0.0;
		d2Udphi2SumM       = 0.0;
		d2Udlambda2SumM    = 0.0;

		d2UdrdphiSumM      = 0.0;
		d2UdrdlambdaSumM   = 0.0;

		d2UdphidlambdaSumM = 0.0;
		d2UdphidrSumM      = 0.0;
		d2UdlambdadrSumM   = 0.0;
		d2UdlambdadphiSumM = 0.0;

        for (m = 0; m <= n; m++){
            j = m+1;

	    // First and second partials of ALFs
	    dPdphi_nm     = P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)] - tanPhi*m*P[IDX2F(k,j,Max_Degree+3)];
	    dPdphi_nmp1   = P[IDX2F(k,j+2,Max_Degree+3)]*scaleFactor[IDX2F(k,j+1,Max_Degree+3)] - tanPhi*(m+1)*(P[IDX2F(k,j+1,Max_Degree+3)]);
	    d2Pdphi2_nm   = dPdphi_nmp1*scaleFactor[IDX2F(k,j,Max_Degree+3)] - m*(secPhi*secPhi)*P[IDX2F(k,j,Max_Degree+3)] - m*tanPhi*dPdphi_nm;
            //printf("\n");
            //printf("dPdphi_nm \n");
            //printf("   %3.20f \n",dPdphi_nm);
            //printf("\n");
            //printf("dPdphi_nmp1 \n");
            //printf("   %3.20f \n",dPdphi_nmp1);
            //printf("\n");
            //printf("d2Pdphi2_nm \n");
            //printf("   %3.20f \n",d2Pdphi2_nm);
	    //dPdphi_nm     = ((reshape(P(k,j+1,:),size(r)).*reshape(scaleFactor(k,j,:),size(r))) - tanPhi.*m.*reshape(P(k,j,:),size(r)));
	    //dPdphi_nmp1   = ((reshape(P(k,j+2,:),size(r)).*reshape(scaleFactor(k,j+1,:),size(r))) - tanPhi.*(m+1).*reshape(P(k,j+1,:),size(r)));
	    //d2Pdphi2_nm   = dPdphi_nmp1.*reshape(scaleFactor(k,j,:),size(r)) - m.*(secPhi.^2).*reshape(P(k,j,:),size(r)) - m.*tanPhi.*dPdphi_nm;

	    // First partials
            dUdrSumM      = dUdrSumM + P[IDX2F(k,j,Max_Degree+3)] *(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("dUdrSumM \n");
            //printf("   %3.20f \n",dUdrSumM);
            dUdphiSumM    = dUdphiSumM + ((P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)]) - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("dUdphiSumM \n");
            //printf("   %3.20f \n",dUdphiSumM);
            dUdlambdaSumM = dUdlambdaSumM + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("dUdlambdaSumM \n");
            //printf("   %3.20f \n",dUdlambdaSumM);
	    // First partials - from Matlab
	    //dUdrSumM      = dUdrSumM + reshape(P(k,j,:),size(r)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
	    //dUdphiSumM    = dUdphiSumM + dPdphi_nm.*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
	    //dUdphiSumM    = dUdphiSumM + ( (reshape(P(k,j+1,:),size(r)).*reshape(scaleFactor(k,j,:),size(r))) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m.*reshape(P(k,j,:),size(r))).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
	    //dUdlambdaSumM = dUdlambdaSumM + m*reshape(P(k,j,:), size(r)).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j));

	    // Second partials
	    d2Udr2SumM         = d2Udr2SumM         + P[IDX2F(k,j,Max_Degree+3)]*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]); 
            //printf("\n");
            //printf("d2Udr2SumM \n");
            //printf("   %3.20f \n",d2Udr2SumM);
	    d2Udphi2SumM       = d2Udphi2SumM       + (d2Pdphi2_nm)*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("d2Udphi2SumM \n");
            //printf("   %3.20f \n",d2Udphi2SumM);
	    d2Udlambda2SumM    = d2Udlambda2SumM    + m*m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*smlambda[j-1] + C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1]);
            //printf("\n");
            //printf("d2Udlambda2SumM");
            //printf("   %3.20f \n",d2Udlambda2SumM);

    	    d2UdrdphiSumM      = d2UdrdphiSumM      + dPdphi_nm*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("d2UdrdphiSumM \n");
            //printf("   %3.20f \n",d2UdrdphiSumM);
	    d2UdrdlambdaSumM   = d2UdrdlambdaSumM   + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("d2UdrdlambdaSumM \n");
            //printf("   %3.20f \n",d2UdrdlambdaSumM);
	    d2UdphidlambdaSumM = d2UdphidlambdaSumM + m*( P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)] - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(S[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] - C[IDX2F(k,j,Max_Degree)]*smlambda[j-1]); 
            //printf("\n");
            //printf("d2UdphidlambdaSumM \n");
            //printf("   %3.20f \n",d2UdphidlambdaSumM);

	    d2UdphidrSumM      = d2UdphidrSumM      + dPdphi_nm*(C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1] + S[IDX2F(k,j,Max_Degree)]*smlambda[j-1]);
            //printf("\n");
            //printf("d2UdphidrSumM \n");
            //printf("   %3.20f \n",d2UdphidrSumM);

	    d2UdlambdadrSumM    = d2UdlambdadrSumM  + m*P[IDX2F(k,j,Max_Degree+3)]*(S[IDX2F(k,j,Max_Degree)]*smlambda[j-1] + C[IDX2F(k,j,Max_Degree)]*cmlambda[j-1]);
//             printf("\n");
//             printf("d2UdlambdadrSumM \n");
//             printf("   %3.20f \n",d2UdlambdadrSumM);
	    d2UdlambdadphiSumM = d2UdlambdadphiSumM + m*( P[IDX2F(k,j+1,Max_Degree+3)]*scaleFactor[IDX2F(k,j,Max_Degree+3)] - z/(sqrt(x*x + y*y))*m*P[IDX2F(k,j,Max_Degree+3)])*(S[IDX2F(k,j,Max_Degree)]*smlambda[j-1] - C[IDX2F(k,j,Max_Degree)]);

	    // Second partials - from Matlab
	    //d2Udr2SumM         = d2Udr2SumM + reshape(P(k,j,:),size(r)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
	    //d2Udphi2SumM = d2Udphi2SumM + (d2Pdphi2_nm).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
	    //d2Udlambda2SumM    = d2Udlambda2SumM + m*m.*reshape(P(k,j,:), size(r)).*(S(k,j).*smlambda(:,j) + C(k,j).*cmlambda(:,j));

	    //d2UdrdphiSumM      = d2UdrdphiSumM      + dPdphi_nm.*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
	    //d2UdrdlambdaSumM   = d2UdrdlambdaSumM + m*reshape(P(k,j,:), size(r)).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j));
	    //d2UdphidlambdaSumM = d2UdphidlambdaSumM + m*( (reshape(P(k,j+1,:),size(r)).*reshape(scaleFactor(k,j,:),size(r))) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m.*reshape(P(k,j,:),size(r))).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j)); 

	    //d2UdphidrSumM      = d2UdphidrSumM + dPdphi_nm.*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j));
	    //d2UdlambdadrSumM    = d2UdlambdadrSumM + m*reshape(P(k,j,:), size(r)).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j));
	    //d2UdlambdadphiSumM = d2UdlambdadphiSumM + m*( (reshape(P(k,j+1,:),size(r)).*reshape(scaleFactor(k,j,:),size(r))) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m.*reshape(P(k,j,:),size(r))).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j)); 

        }

	// First partials
        dUdrSumN      = dUdrSumN      + dUdrSumM*rRatio_n*k;
        dUdphiSumN    = dUdphiSumN    + dUdphiSumM*rRatio_n;
        dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM*rRatio_n;

	// Second partials
	d2Udr2SumN         = d2Udr2SumN + d2Udr2SumM*rRatio_n*k*(k+1);
	d2Udphi2SumN       = d2Udphi2SumN    + d2Udphi2SumM*rRatio_n;
	d2Udlambda2SumN    = d2Udlambda2SumN + d2Udlambda2SumM*rRatio_n;

	d2UdrdphiSumN      = d2UdrdphiSumN + d2UdrdphiSumM*rRatio_n*k;
	d2UdrdlambdaSumN   = d2UdrdlambdaSumN + d2UdrdlambdaSumM*rRatio_n*k;
	d2UdphidlambdaSumN = d2UdphidlambdaSumN + d2UdphidlambdaSumM*rRatio_n;

	d2UdphidrSumN      = d2UdphidrSumN + d2UdphidrSumM*rRatio_n*k;
	d2UdlambdadrSumN   = d2UdlambdadrSumN + d2UdlambdadrSumM*rRatio_n*k;
	d2UdlambdadphiSumN = d2UdlambdadphiSumN + d2UdlambdadphiSumM*rRatio_n;
    }

//     // gravity in spherical coordinates
//     dUdr      = -MU/(radu*radu)*dUdrSumN ;
//     dUdphi    =  MU/radu*dUdphiSumN ;
//     dUdlambda =  MU/radu*dUdlambdaSumN ;
// 
//     //gravity in ECEF coordinates
//     Gxyz[0] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y;
//     Gxyz[1] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x;
//     Gxyz[2] = (1.0/radu)*dUdr*z + ((sqrt(x*x + y*y))/(radu*radu))*dUdphi;

    // Gravity partials in spherical coordinates (second partial of potential)
    d2Udr2         =  MU/(radu*radu*radu)*d2Udr2SumN;
    d2Udphi2       =  MU/radu*d2Udphi2SumN;
    d2Udlambda2    = -MU/radu*d2Udlambda2SumN;

    d2Udrdphi      = -MU/(radu*radu)*d2UdrdphiSumN;
    d2Udrdlambda   = -MU/(radu*radu)*d2UdrdlambdaSumN;
    d2Udphidlambda =  MU/radu*d2UdphidlambdaSumN;

    d2Udphidr      = -MU/(radu*radu)*d2UdphidrSumN;
    d2Udlambdadr   = -MU/(radu*radu)*d2UdlambdadrSumN;
    d2Udlambdadphi =  MU/radu*d2UdlambdadphiSumN;

    // Gravity partials - from Matlab
    //d2Udr2         =  GM./(r.*r.*r).*d2Udr2SumN;
    //d2Udphi2       =  GM./r.*d2Udphi2SumN;
    //d2Udlambda2    = -GM./r.*d2Udlambda2SumN;

    //d2Udrdphi      = -GM./(r.*r).*d2UdrdphiSumN;
    //dUdrdlambda    = -GM./(r.*r).*d2UdrdlambdaSumN;
    //d2Udphidlambda =  GM./r.*d2UdphidlambdaSumN;

    //d2Udphidr      = -GM./(r.*r).*d2UdphidrSumN;
    //d2Udlambdadr   = -GM./(r.*r).*d2UdlambdadrSumN;
    //d2Udlambdadphi =  GM./r.*d2UdlambdadphiSumN;

    //Directly compare gravity partials that should be equal
    d2Udrdphi_diff = d2Udrdphi - d2Udphidr;
    d2Udrdlambda_diff = d2Udrdlambda - d2Udlambdadr;
    d2Udphidlambda_diff = d2Udphidlambda - d2Udlambdadphi;

    // Output below is used to compute matrix G in Perturbed_Accel.m for the state transition matrix
    //delV_del_r_phi_lambda = [dUdr dUdphi dUdlambda d2Udr2 d2Udphi2 d2Udlambda2 d2Udrdphi dUdrdlambda d2Udphidlambda];
    delV_del_r_phi_lambda[0] = dUdr;
    delV_del_r_phi_lambda[1] = dUdphi;
    delV_del_r_phi_lambda[2] = dUdlambda;
    delV_del_r_phi_lambda[3] = d2Udr2;
    delV_del_r_phi_lambda[4] = d2Udphi2;
    delV_del_r_phi_lambda[5] = d2Udlambda2;
    delV_del_r_phi_lambda[6] = d2Udrdphi;
    delV_del_r_phi_lambda[7] = d2Udrdlambda;
    delV_del_r_phi_lambda[8] = d2Udphidlambda;

    // special case for poles
    /*
     atPole = abs(atan2(p[IDX2F(i+1,3,M)],sqrt(p[IDX2F(i+1,1,M)]^2 + p[IDX2F(i+1,2,M)]^2)))==pi/2;
     if any(atPole){
     gx(atPole) = 0;
     gy(atPole) = 0;
     gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
     }
     */
    // gravity in spherical coordinates
    dUdr      = -MU/(radu*radu)*dUdrSumN ;
    dUdphi    =  MU/radu*dUdphiSumN ;
    dUdlambda =  MU/radu*dUdlambdaSumN ;

    //gravity in ECEF coordinates
    Gxyz[0] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*x - (dUdlambda/(x*x + y*y))*y;
    Gxyz[1] = ((1.0/radu)*dUdr - (z/(radu*radu*sqrt(x*x + y*y)))*dUdphi)*y + (dUdlambda/(x*x + y*y))*x;
    Gxyz[2] = (1.0/radu)*dUdr*z + ((sqrt(x*x + y*y))/(radu*radu))*dUdphi;
}
