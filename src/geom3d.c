//
// Created by parthajit on 13/07/21.
//

#include "geom3d.h"

double torad(double deg)
{
    return ((PI * deg) / (double)180.0);

}

double todeg(double rad)
{
    return ((rad * 180) / (double)PI);

}


Vector3d vec3d_cross(Vector3d u, Vector3d v) {
      Vector3d result;
      result.x = (u.y * v.z - u.z * v.y);
      result.y = (u.z * v.x - u.x * v.z);
      result.z = (u.x * v.y - u.y * v.x);
      return result;
}
double vec3d_dot(Vector3d u, Vector3d v) {
      return (u.x * v.x + u.y * v.y + u.z * v.z);
}


double vec3d_norm(Vector3d v, int order) {
      if(order == 2){// as l2 norm is most frequent, this is the first.
	    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
      }else if(order == 1){
	    return (v.x + v.y + v.z);
      }else{
	    double inv = 1.0 /(double)order;
	    double result = pow(v.x, order) + pow(v.y, order) + pow(v.z, order);
	    result = pow(result, inv);
	    return result;
      }
}


Vector3d vec3d_add(Vector3d u, Vector3d v) {
      Vector3d result;
      result.x = u.x + v.x;
      result.y = u.y + v.y;
      result.z = u.z + v.z;
      return result;
}


Vector3d vec3d_sub(Vector3d u, Vector3d v) {
      Vector3d result;
      result.x = u.x - v.x;
      result.y = u.y - v.y;
      result.z = u.z - v.z;
      return result;
}


Vector3d vec3d_neg(Vector3d v) {
      Vector3d t;
      t.x = -v.x;
      t.y = -v.y;
      t.z = -v.z;
      return t;
}

Vector3d vec3d_scal_mult(Vector3d v, double factor){
      v.x = v.x * factor;
      v.y = v.y * factor;
      v.z = v.z * factor;
      return v;
}
Vector3d vec3d_unit(Vector3d v){
      double l2norm = vec3d_norm(v, 2);
      return vec3d_scal_mult(v, 1.0/l2norm);
}

Vector3d vec3d_create(double x, double y, double z) {
      Vector3d v;
      v.x = x;
      v.y = y;
      v.z = z;
      return v;
}

void vec3d_dir_cosine(double* l, double* m, double* n, const Vector3d v){
     double l2norm = vec3d_norm(v, 2);
     *l = v.x / l2norm;
     *m = v.y / l2norm;
     *n = v.z / l2norm;
}

static Vector3d vec3d_matrix_product(Vector3d v, double matrix[][3]){
      Vector3d t;
      t.x = matrix[0][0] * v.x + matrix[0][1] * v.y + matrix[0][2] * v.z;
      t.y = matrix[1][0] * v.x + matrix[1][1] * v.y + matrix[1][2] * v.z;
      t.z = matrix[2][0] * v.x + matrix[2][1] * v.y + matrix[2][2] * v.z;
      return t;
}
Vector3d vec3d_polar_rotation(Vector3d polar_axis, Vector3d reference_normal, double tau){ 
      double l;
      double m;
      double n;
      vec3d_dir_cosine(&l, &m, &n, polar_axis);
      double cos_tau = cos(tau);
      double sin_tau = sin(tau);
      double one_costau = 1.0 - cos_tau;
      double a1 = l * sin_tau;
      double a2 = m * sin_tau;
      double a3 = n * sin_tau;
      double b1 = m * n * one_costau;
      double b2 = n * l * one_costau;
      double b3 = l * m * one_costau;

      double mat[3][3];

      mat[0][0] = cos_tau + l * l * one_costau;
      mat[0][1] = b3 - a3;
      mat[0][2] = b2 + a2;
      mat[1][0] = b3 + a3;
      mat[1][1] = cos_tau + m * m * one_costau;
      mat[1][2] = b1 - a1;
      mat[2][0] = b2 - a2;
      mat[2][1] = b1 + a1;
      mat[2][2] = cos_tau + n * n * one_costau;
      Vector3d t = vec3d_matrix_product(reference_normal, mat);
      return t;
}
double dist(const Point3d p, const Point3d q){
      return sqrt ( (p.x-q.x) * (p.x-q.x) + (p.y-q.y) * (p.y-q.y) + (p.z-q.z) * (p.z-q.z) );
}

double distsqr(const Point3d p, const Point3d q){
      return ( (p.x-q.x) * (p.x-q.x) + (p.y-q.y) * (p.y-q.y) + (p.z-q.z) * (p.z-q.z) );
}

double angle3d(const Point3d p1, const Point3d mid, const Point3d p2){
	const Point3d* a = &mid;
	const Point3d* b = &p1;
	const Point3d* c = &p2;
	double x = (a->x - b->x);
    double y = (a->y - b->y);
    double z = (a->z - b->z);
    double ab = sqrt(x*x + y*y + z*z);

    x = (b->x - c->x);
    y = (b->y - c->y);
    z = (b->z - c->z);
    double bc = sqrt(x*x + y*y + z*z);

    x = (c->x - a->x);
    y = (c->y - a->y);
    z = (c->z - a->z);
    double ac = sqrt(x*x + y*y + z*z);

    double angl =  acos((ab*ab + ac*ac - bc*bc)/(2.0* ab * ac));
    return angl;
}
//double point_trans(Point3d* p, Point3d center){
//      p->x = p->x - center.x
//}
Plane plane_create(Point3d a, Point3d b, Point3d c){
      Plane p;
      p.unit_normal = vec3d_unit(vec3d_cross(vec3d_sub(b, a), vec3d_sub(c, a)));
      p.point = a;
      return p;
}

double plane_perp_dist(Plane plane, Point3d p){
      double d = vec3d_dot(plane.unit_normal, plane.point);
      double sd = vec3d_dot(plane.unit_normal, p);
      double normal_vec_len = 1.0;
      return (sd - d)/normal_vec_len;
}

double dihedral_angle(Plane p1, Plane p2){
      double N = vec3d_dot(p1.unit_normal, p2.unit_normal);
      return acos(N);
}
double torsion_angle(Point3d a, Point3d b, Point3d c, Point3d d){
      Vector3d bc = vec3d_create(c.x - b.x, c.y - b.y, c.z - b.z);

      Plane p1 = plane_create(a, b, c);
      Plane p2 = plane_create(b, c, d);
      double angle = dihedral_angle(p1, p2);

      double sign = vec3d_dot(vec3d_cross(p1.unit_normal, p2.unit_normal), bc);
      if(sign < 0.0){
	    angle = -angle; 
      }
      return angle;
}

void fib_unit_sphere(Point3d* surface, int n){
      double golden = (1.0 + sqrt(5.0)) / 2.0;
      double twopi_by_golden = (2.0 * PI) / golden;
      for(int i=0; i<n; ++i){
	    double theta = twopi_by_golden * i;
	    double phi = acos(1.0 - 2.0 * ((double)i + 0.5)/(double) n);
	    double sinphi = sin(phi);
	    surface[i].x = cos(theta) * sinphi;
	    surface[i].y = sin(theta) * sinphi;
	    surface[i].z = cos(phi);
      }
}


int sphere_jointsurf_pts(Sphere* sphs, int n, Point3d* unitsurf[], int* npts){
      Point3d surfpt;
      int surf = 0;
      int contact = 0;
      for(int i=0; i<n; ++i){
	    surf = surf + npts[i];
	    Point3d* ithsurf = unitsurf[i];
	    Point3d  c = sphs[i].c;
	    double r = sphs[i].r;

	    for(int k=0; k<npts[i]; ++k){
		  surfpt = ithsurf[k];
		  surfpt.x = surfpt.x * r + c.x;
		  surfpt.y = surfpt.y * r + c.y;
		  surfpt.z = surfpt.z * r + c.z;
		  for(int j=0; j<n; ++j){
			if(j != i && distsqr(surfpt, sphs[j].c) < sphs[j].r * sphs[j].r){
			      contact ++;
			      break;
			}
		  }
	    }
      }
      return surf - contact;
}
