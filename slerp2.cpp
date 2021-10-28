// slerp2.cpp: using the recurrence formula to replace the trig functions in the loop

#include "Vector.h"
#include <iostream>
#include <cstdlib>

int main( int argc, char* argv[] ) {

   va::Vector i( 1., 0., 0. ), j( 0., 1., 0. ), k( 0., 0., 1. );
   
   va::Vector u1 = i;
   va::Vector u2 = j;
   if ( argc > 1 ) { // specify initial and final vectors on the command line
   
      u1 = va::Vector( atof( argv[1] ), atof( argv[2] ), atof( argv[3] ) );
      u2 = va::Vector( atof( argv[4] ), atof( argv[5] ), atof( argv[6] ) );
   }
   
   double u12 = u1 * u2;
   va::Vector u0 = ( u2 - u12 * u1 ) / sqrt( ( 1. - u12 ) * ( 1. +  u12 ) );
   
   const int N = 1000;
   const double TH = acos( u12 );
   double th = TH / double( N-1 );
   va::Vector u_2, u_1, u;
   
   u_2 = u1;                                // n = 0 will become u at n - 2;
   u_1 = cos( th ) * u1 + sin( th ) * u0;   // n = 1 will become u at n - 1
   
   std::cout << u_2 << std::endl;
   std::cout << u_1 << std::endl;
   
   const double C = 2. * cos( th );
   
   for ( int n = 2; n < N; n++ ) {
   
      u   = C * u_1 - u_2;
      u_2 = u_1;
      u_1 = u;
      std::cout << u << std::endl;
   }

   return EXIT_SUCCESS;
}
