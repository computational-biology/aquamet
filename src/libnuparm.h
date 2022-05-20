/*
 * =====================================================================================
 *
 *       Filename:  libnuparm.h
 *
 *    Description:  Nuparm library.
 *
 *        Version:  1.0
 *        Created:  Thursday 10 February 2022 11:10:17  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */




#ifndef  __libnuparm_H__
#define  __libnuparm_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#define FALSE (0)
#define TRUE (1)
#define PRM_LINE_MAX (512)
struct bporient{
      int    res1;
      int    res2;
      double open; // rotation
      double buckle;
      double propel;
      double stagger; // translation
      double shear;
      double stretch;
};

struct bpstep{
      int res1;
      int res2;
      double tilt; // rotation
      double roll;
      double twist;
      double shift; // translation
      double slide;
      double rise;
};



struct nuparm{
      struct bporient* orient;
      struct bpstep*   step;
      
      
      
      int size;

      struct{
	    int is_orient;
	    int is_step; 
      }flag;
};

void nuparm_init(struct nuparm* self, int size, const struct nuparm* opt);

void nuparm_read(struct nuparm* self, char* prmfile);
void nuparm_free(struct nuparm* self);





#endif   /* ----- #ifndef __libnuparm_H__  ----- */
