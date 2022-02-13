/*
 * =====================================================================================
 *
 *       Filename:  biodefs.c
 *
 *    Description:   
 *
 *        Version:  1.0
 *        Created:  14/08/21 04:25:35 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "biodefs.h"

char pro_amino[20][2][10]={{"ALA","HYPHOB"}, {"ARG","BASIC"}, {"ASN","POLNUTRL"}, {"ASP", "ACIDIC"}, {"CYS",
                                  "POLNUTRL"},{"GLN", "POLNUTRL"}, {"GLU", "ACIDIC"}, {"GLY","SPL"}, {"HIS","BASIC"}, {"ILE", "HYPHOB"},
        {"LEU", "HYPHOB"}, {"LYS", "BASIC"}, {"MET", "HYPHOB"}, {"PHE","HYPHOB"}, {"PRO","SPL"},
        {"SER", "POLNUTRL"}, {"THR", "POLNUTRL"}, {"TRP","HYPHOB"}, {"TYR", "HYPHOB"}, {"VAL", "HYPHOB"}};

const int pro_amino_size = 20;



int is_std_amino(char* res){
    for (int i = 0; i < pro_amino_size; ++i) {
        if(strcmp(res, pro_amino[i][0]) == 0 ) return 1;
    }
    return 0;
}


char gua_var_names[][4]={
"G",
"  G",
" G ",
"G  ",
"DG ",
" DG",
"DG",
"GUA",
"DG5",
"DG3",
"2MG",
"OMG",
"G7M",
"GNP",
"7MG",
"1MG",
"5CG",
"+G ",
" +G",
"+G",
"GDP",
"M2G",
"  I",
" I ",
"I  ",
"I",
"GMP",
"GTP",
"DG ",
" DG",
"DG",
"PGP",
" YG",
"YG ",
"YG",
"YYG",
"2PR",
"XUG",
" RG",
"RG ",
"RG",
"RG5",
"RG3",
"6OG",
"@@@"
};
int is_guavar(char* res){
      int i=0;
      while(gua_var_names[i][0] != '@'){
	    if(strcmp(gua_var_names[i], res) == 0) return 1;
	    ++i;
      }
      return 0;

}


char ade_var_names[][4] = {
"A",
"  A",
" A ",
"A  ",
"ADE",
"DA5",
"DA3",
"1MA",
"MIA",
"+A ",
" +A",
"+A",
"AMP",
"AMO",
"12A",
"AET",
"PSD",
"AVC",
"APC",
"GOM",
"MAD",
"A23",
"ATP",
"2MA",
"A2M",
"T6A",
"RIA",
"6MZ",
"6IA",
" DA",
"DA ",
"DA",
" RA",
"RA ",
"RA",
"RA5",
"RA3",
"ADP",
"5AA",
"PR5",
"2AD",
"3DA",
"ANZ",
"AVC",
"TSB",
"QSI",
"VAA",
"@@@"
};
int is_adevar(char* res){
      int i=0;
      while(ade_var_names[i][0] != '@'){
	    if(strcmp(ade_var_names[i], res) == 0) return 1;
	    ++i;
      }
      return 0;

}

char cyt_var_names[][4] = {
"C",
"  C",
" C ",
"C  ",
"CYT",
"DC5",
"DC3",
"RC5",
"RC3",
"5MC",
"+C ",
" +C",
"+C",
"OMC",
"S4C",
"CB2",
"5IC",
"CCC",
"1SC",
" DC",
"DC ",
"DC",
" RC",
"RC ",
"RC",
"CBV",
"DCZ",
"CSL",
"CBR",
"C38",
"BLS",
"5CM",
"@@@"
};


int is_cytvar(char* res){
      int i=0;
      while(cyt_var_names[i][0] != '@'){
	    if(strcmp(cyt_var_names[i], res) == 0) return 1;
	    ++i;
      }
      return 0;

}


char ura_var_names[][4] = {
"T",
"  T",
" T ",
"T  ",
" DT",
"DT ",
"DT",
"THY",
"DT5",
"DT3",
"DU5",
"DU3",
"BRU",
"  U",
" U ",
"U  ",
"U",
"URA",
"H2U",
"5MU",
"2MU",
"4SU",
"FMU",
"CMO",
"OMU",
"70U",
" +U",
"+U ",
"+U",
"DHU",
"UR3",
" RT",
"RT ",
"RT",
" RU",
"RU ",
"RU",
"RU5",
"RU3",
"5BU",
"S4U",
"MTU",
"MNU",
"UMS",
" IU",
"IU ",
"IU",
"UD5",
"PYO",
"SUR",
"SSU",
"UCL",
"5IU",
" DU",
"DU ",
"DU",
"PSU",
"@@@"
};


int is_uravar(char* res){
      int i=0;
      while(ura_var_names[i][0] != '@'){
	    if(strcmp(ura_var_names[i], res) == 0) {
	    return 1;
	    }
	    ++i;
      }
      return 0;
}

int is_std_nucleic(char* res){
    if(strcmp(res, "G") == 0) return 1;
    if(strcmp(res, "A") == 0) return 1;
    if(strcmp(res, "C") == 0) return 1;
    if(strcmp(res, "T") == 0) return 1;
    if(strcmp(res, "U") == 0) return 1;
    if(strcmp(res, "DG") == 0) return 1;
    if(strcmp(res, "DA") == 0) return 1;
    if(strcmp(res, "DC") == 0) return 1;
    if(strcmp(res, "DT") == 0) return 1;
    return 0;
}

int is_modi_nucleic(char* res){
      if(is_guavar(res) == 1) return 1;
      if(is_adevar(res) == 1) return 1;
      if(is_cytvar(res) == 1) return 1;
      if(is_uravar(res) == 1) return 1;
      return 0;
}

