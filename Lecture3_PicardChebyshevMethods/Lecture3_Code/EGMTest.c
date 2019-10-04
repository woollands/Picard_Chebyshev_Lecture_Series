/*
   Copyright (c) TAMU LASR Lab, Austin Probe 2014
   Texas A&M University - AEROSPACE department
   File name     : EGMTest.cpp
   Description   : Main script to test EGM 2008 Gravity output
*/

#include "./EGM2008.c"

// main function
int main() {


    double r    = 25512548.0;
    double lat = 0.277878;
    double lon = 4.86994;
    double pos [3] = {0.0};

    pos[0] = 3849366380.185;
    pos[1] = -24230014423.921;
    pos[2] = 6998491653.534;
    //pos [0] = r * cos(lat) * cos(lon);
    //pos [1] = r * cos(lat) * sin(lon);
    //pos [2] = r * sin(lat);


    printf("x \n");
    printf("   %f \n",pos[0]);
    printf("y \n");
    printf("   %f \n",pos[1]);
    printf("z \n");
    printf("   %f \n",pos[2]);
    //int count = 0;

    clock_t EGMstart = clock();            /* get initial time in sec's */

    // Initialize
    double Gxyz [3]                 = {0.0};
    double delV_del_r_phi_lambda[9] = {0.0};

    //for(int i = 0;i<10000;i++)
    //{

        // Calculate the EGM2008 acceleration in the Earth fixed coordinates
        EGM2008(pos, Gxyz, delV_del_r_phi_lambda, 10);
        printf("\n");
        printf("dUdr \n");
        printf("   %3.20f, \n",delV_del_r_phi_lambda[0]);
        printf("\n");
        printf("dUdphi \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[1]);
        printf("\n");
        printf("dUdlambda \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[2]);
        printf("\n");
        printf("d2Udr2 \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[3]);
        printf("\n");
        printf("d2Udphi2 \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[4]);
        printf("\n");
        printf("d2Udlambda2 \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[5]);
        printf("\n");
        printf("d2Udrdphi \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[6]);
        printf("\n");
        printf("d2Udrdlambda \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[7]);
        printf("\n");
        printf("d2Udphidlambda \n");
        printf("   %3.20f \n",delV_del_r_phi_lambda[8]);
        printf(" \n");
    //    count++;
    //}

    double timeEGM = ((double)(clock() - EGMstart)) / (10*(double)CLOCKS_PER_SEC);

    printf("\n");
    printf( "Result: %1.16f\t%1.16f\t%1.16f\n", Gxyz[0],Gxyz[1],Gxyz[2] );
    printf("\n");
    printf( "Norm: %f\n", sqrt(Gxyz[0]*Gxyz[0]+Gxyz[1]*Gxyz[1]+Gxyz[2]*Gxyz[2]) );
    printf( "==========================================================================\n\n" );
    printf( "Computation time for EGM (ms): %lf\n\n", timeEGM );

    return 0;
}
