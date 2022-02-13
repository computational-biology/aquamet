//
// Created by parthajit on 13/07/21.
//

#ifndef BIOGEOM_GEOM3D_H
#define BIOGEOM_GEOM3D_H
#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#define PI (3.1415926)

double torad(double deg);


double todeg(double rad);

typedef struct{
      double x;
      double y;
      double z;
}Point3d;

typedef Point3d Vector3d;

Vector3d vec3d_create(double x, double y, double z);
Vector3d vec3d_cross(Vector3d u, Vector3d v);
double vec3d_dot(Vector3d u, Vector3d v);
double vec3d_norm(Vector3d v, int order);
Vector3d vec3d_add(Vector3d u, Vector3d v);
Vector3d vec3d_sub(Vector3d u, Vector3d v);
Vector3d vec3d_neg(Vector3d v);
Vector3d vec3d_scal_mult(Vector3d v, double factor);
Vector3d vec3d_unit(Vector3d v);

void vec3d_dir_cosine(double* l, double* m, double* n, const Vector3d v);

Vector3d vec3d_polar_rotation(Vector3d polar_axis, Vector3d reference_normal, double tau);


double dist(const Point3d p, const Point3d q);

double distsqr(const Point3d p, const Point3d q);

double angle3d(const Point3d p1, const Point3d mid, const Point3d p2); // returne angle in radian.

//double point_trans(Point3d* p, Point3d center);
typedef struct{
      Vector3d unit_normal;
      Point3d point;
}Plane;
Plane plane_create(Point3d a, Point3d b, Point3d c);

double plane_perp_dist(Plane plane, Point3d p);

double dihedral_angle(Plane p1, Plane p2);   // angle in radian.
double torsion_angle(Point3d a, Point3d b, Point3d c, Point3d d); // angle in radian.



typedef struct{
      Point3d p;
      Vector3d direction;
}Line;


typedef struct{
      Point3d c;
      double r;
}Sphere;



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sphere_fibsurf
 *  Description: This function generates the equidistant points on a unit sphere 
 *  		 centered at (0,0,0). The input 'surface' is an array of size n
 *  		 that will be populated with the surface Point3ds. n is the number 
 *  		 of points to be generated for the surface of the sphere. 
 * =====================================================================================
 */
void fib_unit_sphere(Point3d* surface, int n);


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sphere_jointsurf_pts
 *  Description:  This function computes the joint sueface area of n shperes.
 *        Input:  Array of pointer to surface points array. and number of spheres.
 *       Outout:  Joint surface area. 
 * =====================================================================================
 */

int sphere_jointsurf_pts(Sphere* sphs, int n, Point3d* unitsurf[], int* npts);

#endif//BIOGEOM_GEOM3D_H
