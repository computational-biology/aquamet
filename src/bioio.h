/*
 * =====================================================================================
 *
 *       Filename:  bioio.h
 *
 *    Description: Biological file reading and writing. 
 *
 *        Version:  1.0
 *        Created:  11/08/21 08:31:09 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#ifndef  __bioio_H__
#define  __bioio_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#include "biodefs.h"
#include "geom3d.h"


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  scanpdb
 *  Description:  
 * =====================================================================================
 */

void scanpdb(const char* pdbfile, int (*pf)(char*), const char* chain, 
	    const char* modelno, struct atom** atomary, int* size, enum polymer_type polytype, char occurule);



void print_pdb_line(FILE* fp, const struct atom* atom);
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  scancif
 *  Description:  This function  
 * =====================================================================================
 */

void scancif(const char* ciffile, int (*pf)(char*), const char* chain, const char* modelno, struct atom** atomary, int* size, enum polymer_type polytype, char* label_or_auth, char occurule);

void printpdb(char* file_name, struct atom* atom_array, int size);
void fname_split(char *path, char *basename, char *ext, char *filename) ;

void fname_join(char *filename, const char *path, const char *basename, const char *ext) ;
#endif   /* ----- #ifndef __bioio_H__  ----- */
