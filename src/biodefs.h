/*
 * =====================================================================================
 *
 *       Filename:  biodefs.h
 *
 *    Description:  Definition of rna atoms, residues etc. 
 *
 *        Version:  1.0
 *        Created:  11/08/21 08:11:39 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#ifndef  __biodefs_H__
#define  __biodefs_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geom3d.h"




#define STD_HB_DIST 3.8
#define STD_MG_BOND_DIST 2.7

#define TRUE (1)
#define FALSE (0)

enum polymer_type{
      NUC_TYPE=0,
      PRO_TYPE,
      SOLVENT_TYPE,
      METAL_TYPE,
      OTHER_TYPE,
      ALL_TYPE
};

struct atom{
      // this data is reorganized to 
      // minimize cache miss
      Point3d center;
      char resname[5];
      char chain[5];
      int resid;
      char ins[5];
      char loc[5];
      char symbol[3];
      int id;
      double occu;
      double bfact;
      int model;
      char type; // A for atom H for hetatm.
      char altloc;
};




int is_guavar(char* res);
int is_adevar(char* res);
int is_cytvar(char* res);
int is_uravar(char* res);

int is_std_nucleic(char* res);
int is_modi_nucleic(char* res);

int is_std_amino(char* res);




#endif   /* ----- #ifndef __biodefs_H__  ----- */
