// Rotation.h: Rotation class definition for the algebra of 3D rotations
// Ref: Kuipers, J. B., Quaternions and Rotation Sequences, Princeton, 1999.
//      Altmann, S. L., Rotations, Quaternions, and Double Groups, 1986.
//      Doran & Lasenby, Geometric Algebra for Physicists, Cambridge, 2003.
// R. Saucier, March 2005 (Last revised June 2016)

#ifndef ROTATION_H
#define ROTATION_H

#include "Vector.h"
#include "Random.h"
#include <iostream>

namespace va {   // vector algebra namespace

const Vector DEFAULT_UNIT_VECTOR( 0., 0., 1. );   // arbitrarily choose k
const double DEFAULT_ROTATION_ANGLE( 0. );        // arbitrarily choose 0
const double TWO_PI( 2. * M_PI );

enum ORDER {   // order of rotation sequence about body axes

   // six distinct principal axes factorizations
   ZYX,   // first about z-axis, second about y-axis and third about x-axis
   XYZ,   // first about x-axis, second about y-axis and third about z-axis
   YXZ,   // first about y-axis, second about x-axis and third about z-axis
   ZXY,   // first about z-axis, second about x-axis and third about y-axis
   XZY,   // first about x-axis, second about z-axis and third about y-axis
   YZX,   // first about y-axis, second about z-axis and third about x-axis
   
   // six repeated principal axes factorizations
   ZYZ,   // first about z-axis, second about y-axis and third about z-axis
   ZXZ,   // first about z-axis, second about x-axis and third about z-axis
   YZY,   // first about y-axis, second about z-axis and third about y-axis
   YXY,   // first about y-axis, second about x-axis and third about y-axis
   XYX,   // first about x-axis, second about y-axis and third about x-axis
   XZX    // first about x-axis, second about z-axis and third about x-axis
};

struct quaternion {   // q = w + v, where w is scalar part and v is vector part

   quaternion( void ) {
   }
   
   quaternion( double scalar, Vector vector ) : w( scalar ), v( vector ) {
   }

  ~quaternion( void ) {
   }
   
   // overloaded multiplication of two quaternions
   friend quaternion operator*( const quaternion& q1, const quaternion& q2 ) {
   
      return quaternion(
         ( q1.w * q2.w ) - ( q1.v * q2.v ),                     // scalar part
         ( q1.w * q2.v ) + ( q2.w * q1.v ) + ( q1.v ^ q2.v )    // vector part
      );
   }
   
   double w;   // scalar part
   Vector v;   // vector part (actually a bivector disguised as a vector)
};

struct sequence {   // rotation sequence about three principal body axes
   
   sequence( void ) {
   }
   
   sequence( double phi_1, double phi_2, double phi_3 ) : first( phi_1 ), second( phi_2 ), third( phi_3 ) {
   }
   
  ~sequence( void ) {
   }
   
   // factor quaternion into Euler rotation sequence about three distinct principal axes
   sequence factor( const quaternion& p, ORDER order ) {
      
      double p0 = p.w;
      double p1 = x( p.v );
      double p2 = y( p.v );
      double p3 = z( p.v );
      
      double phi_1, phi_2, phi_3;
      if ( order == ZYX ) {   // distinct principal axes zyx
   
         double A = p0 * p1 + p2 * p3;
         double B = ( p2 - p0 ) * ( p2 + p0 );
         double D = ( p1 - p3 ) * ( p1 + p3 );
         
         phi_3 = atan( - 2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c1 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p1 * c1;
         //double q1 = p1 * c0 - p0 * c1;
         double q2 = p2 * c0 - p3 * c1;
         double q3 = p3 * c0 + p2 * c1;
         
         phi_1 = 2. * atan( q3 / q0 );
         phi_2 = 2. * atan( q2 / q0 );
      }
      else if ( order == XYZ ) {   // distinct principal axes xyz
      
         double A = p1 * p2 - p0 * p3;
         double B = ( p1 - p3 ) * ( p1 + p3 );
         double D = ( p0 - p2 ) * ( p0 + p2 );
         
         phi_3 = atan( - 2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c3 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p3 * c3;
         double q1 = p1 * c0 - p2 * c3;
         double q2 = p2 * c0 + p1 * c3;
         //double q3 = p3 * c0 - p0 * c3;
         
         phi_1 = 2. * atan( q1 / q0 );
         phi_2 = 2. * atan( q2 / q0 );
      }
      else if ( order == YXZ ) {   // distinct principal axes yxz
      
         double A = p1 * p2 + p0 * p3;
         double B = p1 * p1 + p3 * p3;
         double D = -( p0 * p0 + p2 * p2 );
         
         phi_3 = atan( - 2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c3 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p3 * c3;
         double q1 = p1 * c0 - p2 * c3;
         double q2 = p2 * c0 + p1 * c3;
         //double q3 = p3 * c0 - p0 * c3;
         
         phi_1 = 2. * atan( q2 / q0 );
         phi_2 = 2. * atan( q1 / q0 );
      }
      else if ( order == ZXY ) {   // distinct principal axes zxy
      
         double A = p1 * p3 - p0 * p2;
         double B = ( p0 - p1 ) * ( p0 + p1 );
         double D = ( p3 - p2 ) * ( p3 + p2 );
         
         phi_3 = atan( - 2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c2 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p2 * c2;
         double q1 = p1 * c0 + p3 * c2;
         //double q2 = p2 * c0 - p0 * c2;
         double q3 = p3 * c0 - p1 * c2;
         
         phi_1 = 2. * atan( q3 / q0 );
         phi_2 = 2. * atan( q1 / q0 );
      }
      else if ( order == XZY ) {   // distinct principal axes xzy
      
         double A = p1 * p3 + p0 * p2;
         double B = -( p0 * p0 + p1 * p1 );
         double D = p2 * p2 + p3 * p3;
         
         phi_3 = atan( - 2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c2 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p2 * c2;
         double q1 = p1 * c0 + p3 * c2;
         //double q2 = p2 * c0 - p0 * c2;
         double q3 = p3 * c0 - p1 * c2;
         
         phi_1 = 2. * atan( q1 / q0 );
         phi_2 = 2. * atan( q3 / q0 );
      }
      else if ( order == YZX ) {   // distinct principal axes yzx
      
         double A = p2 * p3 - p0 * p1;
         double B = p0 * p0 + p2 * p2;
         double D = -( p1 * p1 + p3 * p3 );
         
         phi_3 = atan( - 2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c1 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p1 * c1;
         //double q1 = p1 * c0 - p0 * c1;
         double q2 = p2 * c0 - p3 * c1;
         double q3 = p3 * c0 + p2 * c1;
         
         phi_1 = 2. * atan( q2 / q0 );
         phi_2 = 2. * atan( q3 / q0 );
      }
      else if ( order == ZYZ ) {   // repeated principal axes zyz
   
         double A = p0 * p1 + p2 * p3;
         double B = -2. * p0 * p2;
         double D = 2. * p1 * p3;
         
         phi_3 = atan( -2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c3 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p3 * c3;
         //double q1 = p1 * c0 - p2 * c3;
         double q2 = p2 * c0 + p1 * c3;
         double q3 = p3 * c0 - p0 * c3;
         
         phi_1 = 2. * atan( q3 / q0 );
         phi_2 = 2. * atan( q2 / q0 );
      }
      else if ( order == ZXZ ) {   // repeated principal axes zxz
   
         double A = p0 * p2 - p1 * p3;
         double B = 2. * p0 * p1;
         double D = 2. * p2 * p3;
         
         phi_3 = atan( -2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c3 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p3 * c3;
         double q1 = p1 * c0 - p2 * c3;
         //double q2 = p2 * c0 + p1 * c3;
         double q3 = p3 * c0 - p0 * c3;
         
         phi_1 = 2. * atan( q3 / q0 );
         phi_2 = 2. * atan( q1 / q0 );
      }
      else if ( order == YZY ) {   // repeated principal axes yzy
   
         double A = p0 * p1 - p2 * p3;
         double B = p0 * p3 + p1 * p2;
         
         phi_3 = atan( -A / B );
         double c0 = cos( 0.5 * phi_3 );
         double c2 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p2 * c2;
         //double q1 = p1 * c0 + p3 * c2;
         double q2 = p2 * c0 - p0 * c2;
         double q3 = p3 * c0 - p1 * c2;
         
         phi_1 = 2. * atan( q2 / q0 );
         phi_2 = 2. * atan( q3 / q0 );
      }
      else if ( order == YXY ) {   // repeated principal axes yxy
   
         double A = p0 * p3 + p1 * p2;
         double B = -2. * p0 * p1;
         double D = 2. * p2 * p3;
         
         phi_3 = atan( -2. * A / ( B + D ) );
         double c0 = cos( 0.5 * phi_3 );
         double c2 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p2 * c2;
         double q1 = p1 * c0 + p3 * c2;
         double q2 = p2 * c0 - p0 * c2;
         //double q3 = p3 * c0 - p1 * c2;
         
         phi_1 = 2. * atan( q2 / q0 );
         phi_2 = 2. * atan( q1 / q0 );
      }
      else if ( order == XYX ) {   // repeated principal axes xyx
   
         double A = p0 * p3 - p1 * p2;
         double B = p0 * p2 + p1 * p3;
         
         phi_3 = atan( -A / B );
         double c0 = cos( 0.5 * phi_3 );
         double c1 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p1 * c1;
         double q1 = p1 * c0 - p0 * c1;
         double q2 = p2 * c0 - p3 * c1;
         //double q3 = p3 * c0 + p2 * c1;
         
         phi_1 = 2. * atan( q1 / q0 );
         phi_2 = 2. * atan( q2 / q0 );
      }
      else if ( order == XZX ) {   // repeated principal axes xzx
   
         double A = p0 * p2 + p1 * p3;
         double B = -p0 * p3 + p1 * p2;
         
         phi_3 = atan( -A / B );
         double c0 = cos( 0.5 * phi_3 );
         double c1 = sin( 0.5 * phi_3 );
   
         double q0 = p0 * c0 + p1 * c1;
         double q1 = p1 * c0 - p0 * c1;
         //double q2 = p2 * c0 - p3 * c1;
         double q3 = p3 * c0 + p2 * c1;
         
         phi_1 = 2. * atan( q1 / q0 );
         phi_2 = 2. * atan( q3 / q0 );
      }
      else {
         std::cerr << "ERROR in Rotation: invalid sequence order: " << order << std::endl;
         exit( EXIT_FAILURE );
      }
      return sequence( phi_1, phi_2, phi_3 );
   }

   double first,    // about first body axis
          second,   // about second body axis
          third;    // about third body axis
};   // end struct sequence

struct matrix {   // all matrices here are rotations in three-space

   // overloaded multiplication of two matrices
   // (defined as a convenience to the user; not used in Rotation class)
   friend matrix operator*( const matrix& A, const matrix& B ) {
   
      matrix C;
      C.a11 = A.a11 * B.a11 + A.a12 * B.a21 + A.a13 * B.a31;
      C.a12 = A.a11 * B.a12 + A.a12 * B.a22 + A.a13 * B.a32;
      C.a13 = A.a11 * B.a13 + A.a12 * B.a23 + A.a13 * B.a33;
      
      C.a21 = A.a21 * B.a11 + A.a22 * B.a21 + A.a23 * B.a31;
      C.a22 = A.a21 * B.a12 + A.a22 * B.a22 + A.a23 * B.a32;
      C.a23 = A.a21 * B.a13 + A.a22 * B.a23 + A.a23 * B.a33;
      
      C.a31 = A.a31 * B.a11 + A.a32 * B.a21 + A.a33 * B.a31;
      C.a32 = A.a31 * B.a12 + A.a32 * B.a22 + A.a33 * B.a32;
      C.a33 = A.a31 * B.a13 + A.a32 * B.a23 + A.a33 * B.a33;
      
      return C;
   }
   
   // transpose of a matrix
   // (defined as a convenience to the user; not used in Rotation class)
   friend matrix transpose( const matrix& A ) {
   
      matrix B;
      B.a11 = A.a11;
      B.a12 = A.a21;
      B.a13 = A.a31;
      
      B.a21 = A.a12;
      B.a22 = A.a22;
      B.a23 = A.a32;
      
      B.a31 = A.a13;
      B.a32 = A.a23;
      B.a33 = A.a33;
      
      return B;
   }
   
   // inverse of a matrix
   // (defined as a convenience to the user; not used in Rotation class)
   friend matrix inverse( const matrix& A ) {

      double det = A.a11 * ( A.a22 * A.a33 - A.a23 * A.a32 ) +
                   A.a12 * ( A.a23 * A.a31 - A.a21 * A.a33 ) +
                   A.a13 * ( A.a21 * A.a32 - A.a22 * A.a31 );
      assert( det != 0. );
      matrix B;
      B.a11 = +( A.a22 * A.a33 - A.a23 * A.a32 ) / det;
      B.a12 = -( A.a12 * A.a33 - A.a13 * A.a32 ) / det;
      B.a13 = +( A.a12 * A.a23 - A.a13 * A.a22 ) / det;
   
      B.a21 = -( A.a21 * A.a33 - A.a23 * A.a31 ) / det;
      B.a22 = +( A.a11 * A.a33 - A.a13 * A.a31 ) / det;
      B.a23 = -( A.a11 * A.a23 - A.a13 * A.a21 ) / det;
   
      B.a31 = +( A.a21 * A.a32 - A.a22 * A.a31 ) / det;
      B.a32 = -( A.a11 * A.a32 - A.a12 * A.a31 ) / det;
      B.a33 = +( A.a11 * A.a22 - A.a12 * A.a21 ) / det;
      
      return B;
   }
   
   // returns the matrix with no more than 15 decimal digit accuracy
   friend matrix set_precision( const matrix& A ) {
      
      matrix B( A );
      if ( fabs(B.a11) < TOL ) B.a11 = 0.0L;
      if ( fabs(B.a12) < TOL ) B.a12 = 0.0L;
      if ( fabs(B.a13) < TOL ) B.a13 = 0.0L;
      
      if ( fabs(B.a21) < TOL ) B.a21 = 0.0L;
      if ( fabs(B.a22) < TOL ) B.a22 = 0.0L;
      if ( fabs(B.a23) < TOL ) B.a23 = 0.0L;
      
      if ( fabs(B.a31) < TOL ) B.a31 = 0.0L;
      if ( fabs(B.a32) < TOL ) B.a32 = 0.0L;
      if ( fabs(B.a33) < TOL ) B.a33 = 0.0L;
      
      return B;
   }
   
   // convenient matrix properties, but not essential to the Rotation class
   
   friend double tr( const matrix& A ) {   // trace of a matrix
   
      return A.a11 + A.a22 + A.a33;
   }

   friend double det( const matrix& A ) {   // determinant of a matrix
   
      return A.a11 * ( A.a22 * A.a33 - A.a23 * A.a32 ) +
             A.a12 * ( A.a23 * A.a31 - A.a21 * A.a33 ) +
             A.a13 * ( A.a21 * A.a32 - A.a22 * A.a31 );
   }
   
   friend Vector eigenvector( const matrix& A ) {   // eigenvector of a matrix
   
      return unit( ( A.a32 - A.a23 ) * Vector( 1., 0., 0. ) +
                   ( A.a13 - A.a31 ) * Vector( 0., 1., 0. ) +
                   ( A.a21 - A.a12 ) * Vector( 0., 0., 1. ) );
   }
   
   friend double angle( const matrix& A ) {   // angle of rotation
   
      return acos( 0.5 * ( A.a11 + A.a22 + A.a33 - 1. ) );
   }

   Vector eigenvector( void ) {   // axis of rotation
   
      return unit( ( a32 - a23 ) * Vector( 1., 0., 0. ) +
                   ( a13 - a31 ) * Vector( 0., 1., 0. ) +
                   ( a21 - a12 ) * Vector( 0., 0., 1. ) );
   }
   
   double tr( void ) {   // trace of the matrix
   
      return a11 + a22 + a33;
   }
   
   double angle( void ) {   // angle of rotation
   
      return acos( 0.5 * ( a11 + a22 + a33 - 1. ) );
   }
   
   double a11, a12, a13,   // 1st row
          a21, a22, a23,   // 2nd row
          a31, a32, a33;   // 3rd row
};   // end struct matrix

class Rotation {

   // friends list

   // overloaded multiplication of two successive rotations (using quaternions)
   // notice the order is important; first the right is applied and then the left
   friend Rotation operator*( const Rotation& R1, const Rotation& R2 ) {

      return Rotation( to_quaternion( R1 ) * to_quaternion( R2 ) );
   }

   // rotation of a vector
   // overloaded multiplication of a vector by a rotation (using quaternions)
   friend Vector operator*( const Rotation& R, const Vector& a ) {
   
      quaternion q( to_quaternion( R ) );
      double w( q.w );
      Vector v( q.v );
   
      Vector b = 2. * ( v ^ a );
      return a + ( w * b ) + ( v ^ b );
   }
   
   // spherical linear interpolation on the unit sphere from u1 to u2
   friend Vector slerp( const Vector& u1, const Vector& u2, double theta, double t ) {
   
      assert( theta != 0 );
      assert( 0. <= t && t <= 1. );
      return ( sin( ( 1. - t ) * theta ) * u1 + sin( t * theta ) * u2 ) / sin( theta );
   }
   
   // spherical linear interpolation on the unit sphere from u1 to u2
   friend Vector slerp( const Vector& u1, const Vector& u2, double t ) {
   
      double theta = angle( u1, u2 );
      if ( theta == 0. ) return u1;
      assert( 0. <= t && t <= 1. );
      return ( sin( ( 1. - t ) * theta ) * u1 + sin( t * theta ) * u2 ) / sin( theta );
   }

   // access functions
   
   // return the unit axial vector
   friend Vector vec( const Rotation& R ) {   // return axial unit eigenvector
   
      return R._vec;
   }
   
   // return the rotation angle
   friend double ang( const Rotation& R ) {   // return rotation angle (rad)
   
      return R._ang;
   }
   
   // inverse rotation
   friend Rotation inverse( Rotation R ) {
   
      return Rotation( R._vec, -R._ang );
   }
   
   // conversion to quaternion
   friend quaternion to_quaternion( const Rotation& R ) {
   
      double a = 0.5 * R._ang;
      Vector u = R._vec;
      
      return quaternion( cos( a ), u * sin( a ) );
   }
   
   // conversion to rotation matrix
   friend matrix to_matrix( const Rotation& R ) {
   
      quaternion q( to_quaternion( R ) );
      double w = q.w;
      Vector v = q.v;
      double v1 = v[ X ], v2 = v[ Y ], v3 = v[ Z ];
      
      matrix A;
      A.a11 = 2. * ( w * w - 0.5 + v1 * v1 );   // 1st row, 1st col
      A.a12 = 2. * ( v1 * v2 - w * v3 );        // 1st row, 2nd col
      A.a13 = 2. * ( v1 * v3 + w * v2 );        // 1st row, 3rd col

      A.a21 = 2. * ( v1 * v2 + w * v3 );        // 2nd row, 1st col
      A.a22 = 2. * ( w * w - 0.5 + v2 * v2 );   // 2nd row, 2nd col
      A.a23 = 2. * ( v2 * v3 - w * v1 );        // 2nd row, 3rd col
      
      A.a31 = 2. * ( v1 * v3 - w * v2 );        // 3rd row, 1st col
      A.a32 = 2. * ( v2 * v3 + w * v1 );        // 3rd row, 2nd col
      A.a33 = 2. * ( w * w - 0.5 + v3 * v3 );   // 3rd row, 3rd col
      
      return A;
   }

   // factor rotation into a rotation sequence
   friend sequence factor( const Rotation& R, ORDER order ) {
   
      sequence s;
      return s.factor( to_quaternion( R ), order );
   }

   // factor matrix representation of rotation into a rotation sequence
   friend sequence factor( const matrix& A, ORDER order ) {
   
      sequence s;
      return s.factor( to_quaternion( Rotation( A ) ), order );
   }

   // overloaded stream operators

   // input a rotation
   friend std::istream& operator>>( std::istream& is, Rotation& R ) {
   
      std::cout << "Specify axis of rotation by entering an axial vector (need not be a unit vector)" << std::endl;
      is >> R._vec;
      std::cout << "Enter the angle of rotation (deg): ";
      is >> R._ang;

      R._vec = unit( R._vec );   // store the unit vector representing the axis
      R._ang = R._ang * D2R;     // store the rotation angle in radians
      return is;
   }
   
   // output a rotation
   friend std::ostream& operator<<( std::ostream& os, const Rotation& R ) {
   
      return os << R._vec << "\t" << R._ang * R2D;
   }

   // output a quaternion
   friend std::ostream& operator<<( std::ostream& os, const quaternion& q ) {
   
      return os << q.w << "\t" << q.v;
   }
   
   // output a matrix
   friend std::ostream& operator<<( std::ostream& os, const matrix& A ) {
  
      matrix B = set_precision( A );   // no more than 15 decimal digits of accuracy
      return os << B.a11 << "\t" << B.a12 << "\t" << B.a13 << std::endl
                << B.a21 << "\t" << B.a22 << "\t" << B.a23 << std::endl
                << B.a31 << "\t" << B.a32 << "\t" << B.a33;
   }

public:

   // constructor from three angles (rad), phi_1, phi_2, phi_3 (in that order, left to right)
   // about three distinct principal body axes
   Rotation( double phi_1, double phi_2, double phi_3, ORDER order ) {
   
      double ang_1 = 0.5 * phi_1, c1 = cos( ang_1 ), s1 = sin( ang_1 );
      double ang_2 = 0.5 * phi_2, c2 = cos( ang_2 ), s2 = sin( ang_2 );
      double ang_3 = 0.5 * phi_3, c3 = cos( ang_3 ), s3 = sin( ang_3 );
      double w;
      Vector v;
             
      if ( order == ZYX ) {   // 1st about z-axis, 2nd about y-axis, 3rd about x-axis (Aerospace sequence)
      
         w = c1 * c2 * c3 + s1 * s2 * s3;
         v = Vector(  c1 * c2 * s3 - s1 * s2 * c3,
                      c1 * s2 * c3 + s1 * c2 * s3,
                     -c1 * s2 * s3 + s1 * c2 * c3 );
      }
      else if ( order == XYZ ) {   // 1st about x-axis, 2nd about y-axis, 3rd about z-axis (FATEPEN sequence)
      
         w = c1 * c2 * c3 - s1 * s2 * s3;
         v = Vector(  c1 * s2 * s3 + s1 * c2 * c3,
                      c1 * s2 * c3 - s1 * c2 * s3,
                      c1 * c2 * s3 + s1 * s2 * c3 );
      }
      else if ( order == YXZ ) {   // 1st about y-axis, 2nd about x-axis, 3rd about z-axis
         
         w = c1 * c2 * c3 + s1 * s2 * s3;
         v = Vector(  c1 * s2 * c3 + s1 * c2 * s3,
                     -c1 * s2 * s3 + s1 * c2 * c3,
                      c1 * c2 * s3 - s1 * s2 * c3 );
      }
      else if ( order == ZXY ) {   // 1st about z-axis, 2nd about x-axis, 3rd about y-axis
            
         w = c1 * c2 * c3 - s1 * s2 * s3;
         v = Vector(  c1 * s2 * c3 - s1 * c2 * s3,
                      c1 * c2 * s3 + s1 * s2 * c3,
                      c1 * s2 * s3 + s1 * c2 * c3 );
      }
      else if ( order == XZY ) {   // 1st about x-axis, 2nd about z-axis, 3rd about y-axis
            
         w = c1 * c2 * c3 + s1 * s2 * s3;
         v = Vector( -c1 * s2 * s3 + s1 * c2 * c3,
                      c1 * c2 * s3 - s1 * s2 * c3,
                      c1 * s2 * c3 + s1 * c2 * s3 );
      }
      else if ( order == YZX ) {   // 1st about y-axis, 2nd about z-axis, 3rd about x-axis
            
         w = c1 * c2 * c3 - s1 * s2 * s3;
         v = Vector(  c1 * c2 * s3 + s1 * s2 * c3,
                      c1 * s2 * s3 + s1 * c2 * c3,
                      c1 * s2 * c3 - s1 * c2 * s3 );
      }
      else if ( order == ZYZ ) {   // Euler sequence, 1st about z-axis, 2nd about y-axis, 3rd about z-axis
      
         w = c1 * c2 * c3 - s1 * c2 * s3;
         v = Vector(  c1 * s2 * s3 - s1 * s2 * c3,
                      c1 * s2 * c3 + s1 * s2 * s3,
                      c1 * c2 * s3 + s1 * c2 * c3 );
      }
      else if ( order == ZXZ ) {   // Euler sequence, 1st about z-axis, 2nd about x-axis, 3rd about z-axis
      
         w = c1 * c2 * c3 - s1 * c2 * s3;
         v = Vector(  c1 * s2 * c3 + s1 * s2 * s3,
                     -c1 * s2 * s3 + s1 * s2 * c3,
                      c1 * c2 * s3 + s1 * c2 * c3 );
      }
      else if ( order == YZY ) {   // Euler sequence, 1st about y-axis, 2nd about z-axis, 3rd about y-axis
      
         w = c1 * c2 * c3 - s1 * c2 * s3;
         v = Vector( -c1 * s2 * s3 + s1 * s2 * c3,
                      c1 * c2 * s3 + s1 * c2 * c3,
                      c1 * s2 * c3 + s1 * s2 * s3 );
      }
      else if ( order == YXY ) {   // Euler sequence, 1st about y-axis, 2nd about x-axis, 3rd about y-axis
      
         w = c1 * c2 * c3 - s1 * c2 * s3;
         v = Vector(  c1 * s2 * c3 + s1 * s2 * s3,
                      c1 * c2 * s3 + s1 * c2 * c3,
                      c1 * s2 * s3 - s1 * s2 * c3 );
      }
      else if ( order == XYX ) {   // Euler sequence, 1st about x-axis, 2nd about y-axis, 3rd about x-axis
      
         w = c1 * c2 * c3 - s1 * c2 * s3;
         v = Vector(  c1 * c2 * s3 + s1 * c2 * c3,
                      c1 * s2 * c3 + s1 * s2 * s3,
                     -c1 * s2 * s3 + s1 * s2 * c3 );
      }
      else if ( order == XZX ) {   // Euler sequence, 1st about x-axis, 2nd about z-axis, 3rd about x-axis
      
         w = c1 * c2 * c3 - s1 * c2 * s3;
         v = Vector(  c1 * c2 * s3 + s1 * c2 * c3,
                      c1 * s2 * s3 - s1 * s2 * c3,
                      c1 * s2 * c3 + s1 * s2 * s3 );
      }
      else {
      
         std::cerr << "ERROR in Rotation: invalid order: " << order << std::endl;
         exit( EXIT_FAILURE );
      }
      if ( w >= 1. || v == 0. ) {
         _ang = DEFAULT_ROTATION_ANGLE;
         _vec = DEFAULT_UNIT_VECTOR;
      }
      else {
         _ang = 2. * acos( w );
         _vec = v / sqrt( 1. - w * w );
      }
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }
   
   // constructor from a rotation sequence
   Rotation( const sequence& s, ORDER order ) {
   
      Rotation R( s.first, s.second, s.third, order );
      
      _vec = vec( R );     // set the axial vector
      _ang = ang( R );     // set the rotation angle
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }

   // constructor from an axial vector and rotation angle (rad)
   Rotation( const Vector& v, double a ) : _vec( v ), _ang( a ) {
   
      _vec = _vec.unit();   // store the unit vector representing the axis
      _set_angle();    // angle in the range [-M_PI, M_PI]
   }
   
   // constructor using sphericalCoord (of axial vector) and rotation angle (rad)
   //Rotation( rng::sphericalCoord s, double ang ) {
	Rotation( std::pair<double,double> s, double ang ) {
   
      //_vec = Vector( 1., s.theta, s.phi, POLAR );   // unit vector
		_vec = Vector( 1., s.first, s.second, POLAR );   // unit vector
      _ang = ang;
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }
   
   // constructor from the cross product of two vectors
   // generate the rotation that, when applied to vector a, will result in vector b
   Rotation( const Vector& a, const Vector& b ) {
   
      _vec = unit( a ^ b );   // unit vector
      double s = a.unit() * b.unit();
      if ( s >= 1. )
         _ang = 0.;
      else if ( s <= -1. )
         _ang= M_PI;
      else
         _ang = acos( s );
      
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }
   
   // constructor from unit quaternion
   Rotation( const quaternion& q ) {
   
      double w = q.w;
      Vector v = q.v;
      
      if ( w >= 1. || v == 0. ) {
         _ang = DEFAULT_ROTATION_ANGLE;
         _vec = DEFAULT_UNIT_VECTOR;
      }
      else {
         double n = sqrt( w * w + v * v );   // need to insure it's a unit quaternion
         w /= n;
         v /= n;
         _ang = 2. * acos( w );
         _vec = v / sqrt( 1. - w * w );
      }
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }
   
   // constructor from rotation matrix
   Rotation( const matrix& A ) {
   
      _vec = ( A.a32 - A.a23 ) * Vector( 1., 0., 0. ) +
             ( A.a13 - A.a31 ) * Vector( 0., 1., 0. ) +
             ( A.a21 - A.a12 ) * Vector( 0., 0., 1. );
      if ( _vec == 0. ) {   // then it must be the identity matrix
         _vec = DEFAULT_UNIT_VECTOR;
         _ang = DEFAULT_ROTATION_ANGLE;
      }
      else {
         _vec = _vec.unit();
         _ang = acos( 0.5 * ( A.a11 + A.a22 + A.a33 - 1. ) );
      }
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }

   // constructor from two sets of three vectors, where the pair must be related by a pure rotation
   // returns the rotation that will take a1 to b1, a2 to b2, and a3 to b3
   // Ref: Micheals, R. J. and Boult, T. E., "Increasing Robustness in Self-Localization and Pose Estimation," online paper.

   Rotation( const Vector& a1, const Vector& a2, const Vector& a3,      // initial vectors
             const Vector& b1, const Vector& b2, const Vector& b3 ) {   // rotated vectors
             
      assert( det( a1, a2, a3 ) != 0. && det( b1, b2, b3 ) != 0. );   // all vectors must be nonzero
      assert( fabs( det( a1, a2, a3 ) - det( b1, b2, b3 ) ) < 0.001 );   // if it doesn't preserve volume, it's not a pure rotation

      if ( det( a1, a2, a3 ) == 1 ) {   // these are basis vectors so use simpler method to construct the rotation
      
         matrix A;
         A.a11 = b1 * a1; A.a12 = b2 * a1; A.a13 = b3 * a1;
         A.a21 = b1 * a2; A.a22 = b2 * a2; A.a23 = b3 * a2;
         A.a31 = b1 * a3; A.a32 = b2 * a3; A.a33 = b3 * a3;
      
         _vec = ( A.a32 - A.a23 ) * Vector( 1., 0., 0. ) +
                ( A.a13 - A.a31 ) * Vector( 0., 1., 0. ) +
                ( A.a21 - A.a12 ) * Vector( 0., 0., 1. );
         if ( _vec == 0. ) {   // then it must be the identity matrix
            _vec = DEFAULT_UNIT_VECTOR;
            _ang = DEFAULT_ROTATION_ANGLE;
         }
         else {
            _vec = _vec.unit();
            _ang = acos( 0.5 * ( A.a11 + A.a22 + A.a33 - 1. ) );
         }
         _set_angle();   // angle in the range [-M_PI, M_PI]
      }
      else {   // use R. J. Micheals' closed-form solution to the absolute orientation problem
      
         Vector c1, c2, c3;
         double aaa, baa, aba, aab, caa, aca, aac;
         double q02, q0, q0q1, q1, q0q2, q2, q0q3, q3;
   
         aaa = det( a1, a2, a3 );
         baa = det( b1, a2, a3 );
         aba = det( a1, b2, a3 );
         aab = det( a1, a2, b3 );
   
         q02 = fabs( ( aaa + baa + aba + aab ) / ( 4. * aaa ) );
         q0 = sqrt( q02 );

         c1 = Vector( 0., b1[Z], -b1[Y] );
         c2 = Vector( 0., b2[Z], -b2[Y] );
         c3 = Vector( 0., b3[Z], -b3[Y] );
   
         caa = det( c1, a2, a3 );
         aca = det( a1, c2, a3 );
         aac = det( a1, a2, c3 );
   
         q0q1 = ( caa + aca + aac ) / ( 4. * aaa );
         q1 = q0q1 / q0;
   
         c1 = Vector( -b1[Z], 0., b1[X] );
         c2 = Vector( -b2[Z], 0., b2[X] );
         c3 = Vector( -b3[Z], 0., b3[X] );
   
         caa = det( c1, a2, a3 );
         aca = det( a1, c2, a3 );
         aac = det( a1, a2, c3 );
   
         q0q2 = ( caa + aca + aac ) / ( 4. * aaa );
         q2 = q0q2 / q0;
   
         c1 = Vector( b1[Y], -b1[X], 0. );
         c2 = Vector( b2[Y], -b2[X], 0. );
         c3 = Vector( b3[Y], -b3[X], 0. );
   
         caa = det( c1, a2, a3 );
         aca = det( a1, c2, a3 );
         aac = det( a1, a2, c3 );
   
         q0q3 = ( caa + aca + aac ) / ( 4. * aaa );
         q3 = q0q3 / q0;
         
         // no need to normalize since constructed to be unit quaternions
         double w( q0 );
         Vector v( q1, q2, q3 );

         if ( w >= 1. || v == 0. ) {
            _ang = DEFAULT_ROTATION_ANGLE;
            _vec = DEFAULT_UNIT_VECTOR;
         }
         else {
            _ang = 2. * acos( w );
            _vec = v / sqrt( 1. - w * w );
         }
         _set_angle();   // angle in the range [-M_PI, M_PI]
      }
   }

   // constructor for a uniform random rotation, uniformly distributed over the unit sphere,
   // by fast generation of random quaternions, uniformly-distributed over the 4D unit sphere
   // Ref: Shoemake, K., "Uniform Random Rotations," Graphic Gems III, September, 1991.
   Rotation( rng::Random& rng ) {   // random rotation in canonical form

      double s   = rng.uniform( 0., 1. );
      double s1  = sqrt( 1. - s );
      double th1 = rng.uniform( 0., TWO_PI );
      double x   = s1 * sin( th1 );
      double y   = s1 * cos( th1 );
      double s2  = sqrt( s );
      double th2 = rng.uniform( 0., TWO_PI );
      double z   = s2 * sin( th2 );
      double w   = s2 * cos( th2 );
      Vector v( x, y, z );
      
      if ( w >= 1. || v == 0. ) {
         _ang = DEFAULT_ROTATION_ANGLE;
         _vec = DEFAULT_UNIT_VECTOR;
      }
      else {
         _ang = 2. * acos( w );
         _vec = v / sqrt( 1. - w * w );
      }
      _set_angle();   // angle in the range [-M_PI, M_PI]
   }

   // default constructor
   Rotation( void ) {
   
      _vec = DEFAULT_UNIT_VECTOR;
      _ang = DEFAULT_ROTATION_ANGLE;
   }

   // default destructor
  ~Rotation( void ) {
   }
   
   // copy constructor
   Rotation( const Rotation& r ) : _vec( r._vec ), _ang( r._ang ) {

      _set_angle();   // angle in the range [-M_PI, M_PI]
   }
   
   // overloaded assignment operator
   Rotation& operator=( const Rotation& R ) {
   
      if ( this != &R ) {
         _vec = R._vec;
         _ang = R._ang;
         _set_angle();   // angle in the range [-M_PI, M_PI]
      }
      return *this;
   }
   
   // conversion operator to return the eigenvector vec
   operator Vector( void ) const {
   
      return _vec;
   }

   // conversion operator to return the angle of rotation about the eigenvector
   operator double( void ) const {
   
      return _ang;
   }
   
   // overloaded arithmetic operators
   
   // inverse rotation
   Rotation operator-( void ) {
   
      return Rotation( -_vec, _ang );
   }
   
   // triple scalar product, same as a * ( b ^ c )
   inline double det( const Vector& a, const Vector& b, const Vector& c ) {

      return a[X] * ( b[Y] * c[Z] - b[Z] * c[Y] ) +
             a[Y] * ( b[Z] * c[X] - b[X] * c[Z] ) +
             a[Z] * ( b[X] * c[Y] - b[Y] * c[X] );
   }
   
private:

   inline void _set_angle( void ) {   // always choose the smaller of the two angles
   
      if ( _ang > +TWO_PI ) _ang -= TWO_PI;
      if ( _ang < -TWO_PI ) _ang += TWO_PI;
      if ( _ang > M_PI ) {
         _ang = TWO_PI - _ang;
         _vec = -_vec;
      }
      else if ( _ang < -M_PI ) {
         _ang = TWO_PI + _ang;
         _vec = -_vec;
      }
   }

   Vector _vec;   // unit eigenvector representing the axis of rotation
   double _ang;   // angle of rotation (rad) falls in the range [-M_PI, M_PI]
};

// declaration of friends
   quaternion operator*( const quaternion& q1, const quaternion& q2 );           // product of two quaternions
   matrix operator*( const matrix& A, const matrix& B );                         // product of two matrices, first B, then A
   matrix transpose( const matrix& A );                                          // transpose of a matrix
   matrix inverse( const matrix& A );                                            // inverse of a matrix
   matrix set_precision( const matrix& A );                                      // returns a matrix with no more than 15 decimal digit accuracy
   double tr( const matrix& A );                                                 // trace of a matrix
   double det( const matrix& A );                                                // determinant of a matrix
   Vector eigenvector( const matrix& A );                                        // eigenvector of a rotation matrix
   double angle( const matrix& A );                                              // angle of rotation (rad)
   Rotation operator*( const Rotation& R1, const Rotation& R2 );                 // successive rotations, first right, then left
   Vector operator*( const Rotation& R, const Vector& a );                       // rotation of a vector
   Vector slerp( const Vector& u1, const Vector& u2, double t );                 // spherical linear interpolation on the unit sphere
   Vector slerp( const Vector& u1, const Vector& u2, double theta, double t );   // slerp, given the angle between the vectors
   Vector vec( const Rotation& R );                                              // return axial unit eigenvector
   double ang( const Rotation& R );                                              // return rotation angle (rad)
   Rotation inverse( Rotation R );                                               // inverse rotation
   quaternion to_quaternion( const Rotation& R );                                // convert rotation to a quaternion
   matrix to_matrix( const Rotation& R );                                        // convert a rotation to a rotation matrix
   sequence factor( const Rotation& R, ORDER order );                            // factor a rotation into a rotation sequence
   sequence factor( const matrix& A, ORDER order );                              // factor a rotation matrix into a rotation sequence
   std::istream& operator>>( std::istream& is, Rotation& R );                    // input rotation
   std::ostream& operator<<( std::ostream& os, const Rotation& R );              // output a rotation
   std::ostream& operator<<( std::ostream& os, const quaternion& q );            // output a quaternion
   std::ostream& operator<<( std::ostream& os, const matrix& A );                // output a rotation matrix

} // vector algebra namespace
#endif
