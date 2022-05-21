/*
 * =====================================================================================
 *
 *       Filename:  fetch.h
 *
 *    Description:  This library fetches the strusture files from RCSB.
 *
 *        Version:  1.0
 *        Created:  Wednesday 04 May 2022 06:12:43  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */






#ifndef  __fetch_H__
#define  __fetch_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<unistd.h>

#define FETCH_BUFF_SIZE (128)
int check_internet_conn(void);
void fetch_rcsb(char* accn, char* ext);

#endif   /* ----- #ifndef __fetch_H__  ----- */
