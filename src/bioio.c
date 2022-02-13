/*
 * =====================================================================================
 *
 *       Filename:  bioio.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  14/08/21 08:05:09 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "bioio.h"
struct cif_attrloc{
      int grouploc;  // group_PDB
      int type_symbloc;

      int idloc;   // id
      int locloc;  // atom_id
      int altloc;  // alt_id
      int resnameloc;  // comp_id
      int chainloc; // asym_id
      int residloc; // seq_id
      int modelloc; // pdbx_PDB_model_num
      int occuloc; //occupancy
      int insloc;
      int bfactloc;
      int xloc;
      int yloc;
      int zloc;

};

static void IUPAC_atom_loc(char* iupac, const char* loc, const char* symb){
      int loclen = strlen(loc);
      if(loclen == 4){
	    strcpy(iupac, loc);
	    return;
      }
      int symblen = strlen(symb);
      char tmploc[5];
      strcpy(tmploc, loc);
      char* symbindx = strstr(tmploc, symb);
      int symbpos = (symbindx - tmploc)/sizeof(char);
      if(symbpos == 0){
	    if(symblen == 1){
		  iupac[0] = ' ';
		  strcpy(iupac+1, loc);
	    }else{
		  strcpy(iupac, loc);
	    }
      }else{
	    strcpy(iupac, loc);
      }
      int iupaclen = strlen(iupac);
      for(int i=3; i>=iupaclen; --i){
	    iupac[i] = ' ';
      }
      iupac[4] = '\0';
}



void print_pdb_line(FILE* fp, const struct atom* atom){
      char tag[8];
      char iupacloc[6];
      char ins;


      if(atom->type == 'A'){
	    strcpy(tag, "ATOM");
      }else{
	    strcpy(tag, "HETATM");
      }
      IUPAC_atom_loc(iupacloc, atom->loc, atom->symbol);
      if(atom->ins[0] == '?' || atom->ins[0] == '.'){
	    ins = ' ';
      }else{
	    ins = atom->ins[0];
      }


      fprintf(fp, "%-6s%5d %4s%c%3s%2s%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s\n", 
		  tag, atom->id, iupacloc, atom->altloc,
		  atom->resname, atom->chain, atom->resid, ins, 
		  atom->center.x, atom->center.y, atom->center.z,
		  atom->occu, atom->bfact, atom->symbol);


      //      fprintf(stdout, "%6s |%5d|%4s|%c|%3s|%2s|%4d|%c|   %8.3lf|%8.3lf|%8.3lf|%6.2lf|%6.2lf          %2s\n", 
      //		  tag, atom->id, iupacloc, atom->altloc,
      //		  atom->resname, atom->chain, atom->resid, ins, 
      //		  atom->center.x, atom->center.y, atom->center.z,
      //		  atom->occu, atom->bfact, atom->symbol);

}


void printpdb(char* file_name, struct atom* atom_array, int size){

      FILE	*fp;										/* output-file pointer */

      fp	= fopen( file_name, "w" );
      if ( fp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      for(int i=0; i<size; ++i){
	    print_pdb_line(fp, atom_array+i);
      }
      fprintf(fp, "END\n");

      if( fclose(fp) == EOF ) {			/* close output file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }

}


static void insert_atom(struct atom** atomarray, int* index, int* max_size, struct atom atom){
      struct atom* atomarr = *atomarray;
      if(*index == *max_size){
	    *max_size *= 2;
	    atomarr = (struct atom*) realloc ( atomarr, *max_size * sizeof(struct atom) );
	    *atomarray = atomarr;

	    if ( atomarr == NULL ) {
		  fprintf ( stderr, "\ndynamic memory reallocation failed\n" );
		  exit (EXIT_FAILURE);
	    }
      }
      atomarr[*index] = atom;
      (*index) ++;
}


static int guess_size_pdb(FILE* fp, char* line, enum polymer_type polytype){
      int max_size = 0;
      int buffer = 2000;
      int proatm = 0;
      int nucatm = 0;
      int hohatm = 0;
      int metatm = 0;
      while(fgets(line, 200, fp) != NULL){
	    if(strncmp(line+13, "PROTEIN ATOMS", 13) == 0){
		  proatm = atol(line+40);
	    }else if(strncmp(line+13, "NUCLEIC ACID ATOMS", 18) == 0){
		  nucatm = atol(line+40);
	    }else if(strncmp(line+13, "HETEROGEN ATOMS", 15) == 0){
		  metatm = atol(line+40);
	    }else if(strncmp(line+13, "SOLVENT ATOMS", 13) == 0){
		  hohatm = atol(line+40);
	    }else if(strncmp(line, "ATOM", 4) == 0){
		  break;
	    }else if(strncmp(line, "HETATM", 6) == 0){
		  break;
	    }else if(strncmp(line, "MODEL", 5) == 0){
		  break;
	    }else{
		  continue;
	    }
      }
      if(polytype == NUC_TYPE){
	    max_size = nucatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == PRO_TYPE){
	    max_size = proatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == SOLVENT_TYPE){
	    max_size = nucatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == METAL_TYPE){
	    max_size = metatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == ALL_TYPE){
	    max_size = nucatm + proatm + hohatm + metatm;  // buffer is extra space for some extra atoms.
      }else{
	    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Wrong molecule type selected.\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }


      if(max_size == 0){
	    //if no such record present in the pdb file. 
	    max_size = 10000;
      }else{
	    max_size += buffer;
      }
      return max_size;

}

static int locate_model_pdb(FILE* fp, char* line, const char* modelno){
      int model = -1;
      int found = 0;
      if(modelno != NULL){
	    model = atoi(modelno);
	    if(model<=0){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s. (File: %s, Line %d)... model no cannot be zero or -ve.\n", __func__, __FILE__, __LINE__);
		  exit(EXIT_FAILURE);
	    }
      }
      int currmodel = -1;
      do{
	    if(strncmp(line, "MODEL", 5) != 0) continue;
	    sscanf(line+6, "%8d", &currmodel);
	    if(currmodel == model){
		  found = 2;
		  break;
	    }
      }while(fgets(line, 200, fp) != NULL);
      if(found != 2){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Requested model not found.\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }
      fgets(line, 200, fp);
      return model;
}

static struct atom parse_pdb_line_to_atom(const char* line, int model){

      char token[20];
      struct atom atom;
      atom.type = line[0];
      strncpy(token, line+6,5);
      token[5] = '\0';
      atom.id = atoi(token);

      strncpy(token, line+12,4);
      token[4] = '\0';
      sscanf(token, "%s", atom.loc);

      atom.altloc = line[16];

      strncpy(token, line+17,3);
      token[3] = '\0';
      sscanf(token, "%s", atom.resname);


      strncpy(token, line+20,2);
      token[2] = '\0';
      sscanf(token, "%s", atom.chain);

      strncpy(token, line+22,4);
      token[4] = '\0';
      atom.resid = atoi(token);

      atom.ins[0] = line[26];
      atom.ins[1] = '\0';

      strncpy(token, line+30,8);
      token[8] = '\0';
      atom.center.x = atof(token);

      strncpy(token, line+38,8);
      token[8] = '\0';
      atom.center.y = atof(token);

      strncpy(token, line+46,8);
      token[8] = '\0';
      atom.center.z = atof(token);

      strncpy(token, line+54,6);
      token[6] = '\0';
      atom.occu = atof(token);

      strncpy(token, line+60,6);
      token[6] = '\0';
      atom.bfact = atof(token);

      strncpy(token, line+76,2);
      token[2] = '\0';
      sscanf(token, "%s", atom.symbol);
      atom.model = model;
if(strncmp(line, "HETATM52148 MG", 14) == 0){
    printf("id=%d, and resname =%s|\n", atom.id, atom.resname);
    }
      return atom;
}

static int pdb_multi_occupancy(struct atom* atom, struct atom* buffer, struct atom* current, char rule){
      int select = 1;
      int reject = 0;
      if(  (buffer->altloc == ' ' || buffer->altloc == '.') ){
	    *current = *buffer;
	    *buffer  = *atom;
	    if(rule == 'S' && current->occu <= 0.499) { 
		  return reject; 
	    }else{
		  return select;
	    }
      }else if(rule == 'A'){
	    // all occupancy
	    *current = *buffer;
	    *buffer  = *atom;
	    return select;
      }

      if(  strcmp(atom->loc, buffer->loc) == 0 && 
		  atom->resid == buffer->resid && 
		  strcmp(atom->chain, buffer->chain) ==0 && 
		  strcmp(atom->ins, buffer->ins) == 0){

	    if(atom->occu > buffer->occu){
		  *current = *buffer;
		  *buffer  = *atom;
	    }
	    return reject;
      }else{
	    *current = *buffer;
	    *buffer  = *atom;
	    if(rule == 'B'){
		  // best occupancy
		  return select;
	    }else{
		  // standard occupancy
		  if(current->occu <= 0.499) return reject;
		  return select;
	    }
      }
}


void scanpdb(const char* pdbfile, int (*pf)(char*), const char* chain, 
	    const char* modelno, struct atom** atomary, int* size, enum polymer_type polytype, char occurule){

      FILE* fp	= fopen( pdbfile, "r" );
      if ( fp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			pdbfile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      char line[200];
      int max_size = guess_size_pdb(fp, line, polytype);
      int model = -1;
      if(modelno != NULL){
	    model = locate_model_pdb(fp, line, modelno); 
      }
      struct atom* atomarr = NULL;
      int atmindx = 0;
      struct atom atom;
      struct atom buffer = parse_pdb_line_to_atom(line, model);
      while(!  ((modelno == NULL || buffer.model == model) && ((chain == NULL)||(strcmp(chain, buffer.chain) == 0)) && pf(buffer.resname) == 1)){
	    if(fgets(line, 1024, fp) == NULL){
		  goto final;
	    }
	    buffer = parse_pdb_line_to_atom(line, model);
      }
      struct atom current;
      atomarr = (struct atom*) malloc (max_size * sizeof(struct atom));

      while(fgets(line, 200, fp) != NULL && strncmp(line, "END", 3) != 0 && strncmp(line, "ENDMDL", 6) != 0){
	    if(strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0){
		  continue;
	    }

	    atom = parse_pdb_line_to_atom(line, model);
	    if(  (pf(atom.resname) == 1)     &&    ((chain == NULL) || (strcmp(chain,atom.chain) == 0))){
		  int status = pdb_multi_occupancy(&atom, &buffer, &current, occurule);
		  if(status == 1){
			//print_pdb_line(stdout, &current);
			insert_atom(&atomarr, &atmindx, &max_size, current);
		  }
	    }
      }

      if((pf(buffer.resname) == 1) &&((chain == NULL) || (strcmp(chain,buffer.chain)== 0))){
	    struct atom dummy;
	    dummy.occu = ' ';
	    int status = pdb_multi_occupancy(&dummy, &buffer, &current, occurule);
	    if(status == 1){
		  //print_pdb_line(stdout, &current);
		  insert_atom(&atomarr, &atmindx, &max_size, current);
	    }
      }
final:
      *size = atmindx;
      *atomary = atomarr;
      if( fclose(fp) == EOF ) {			/* close input file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			pdbfile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
}


static void init_cif_attrib_loc(struct cif_attrloc* attribloc){
      attribloc->grouploc=-1;  // group_PDB
      attribloc->type_symbloc = -1;
      attribloc->idloc = -1;   // id
      attribloc->locloc = -1;  // label_atom_id
      attribloc->altloc = -1;  // label_alt_id
      attribloc->resnameloc = -1;  // label_comp_id
      attribloc->chainloc = -1; // auth_asym_id
      attribloc->residloc = -1; // auth_seq_id
      attribloc->modelloc = -1; // pdbx_PDB_model_num
      attribloc->occuloc  = -1; //occupancy
      attribloc->insloc   = -1;
      attribloc->bfactloc = -1;
      attribloc->xloc = -1;
      attribloc->yloc = -1;
      attribloc->zloc = -1;
}
static void set_cif_attrib_loc(struct cif_attrloc* attribloc, char* line, char* type, FILE* fp){
      int count = 0;
      attribloc->grouploc = count;

      count++;

      if(strcmp(type, "auth") == 0){
	    while(fgets(line, 1024, fp) != NULL   &&   strncmp(line, "_atom_site.", 11) == 0){
		  if(strncmp(line+11,"id", 2) == 0 && line[13] != '_' ){
			attribloc->idloc = count;
		  }else if(strncmp(line+11,"type_symbol", 11) == 0 && line[22] != '_' ){
			attribloc->type_symbloc = count;
		  }else if(strncmp(line+11,"auth_atom_id", 12) == 0 && line[23] != '_' ){
			attribloc->locloc = count;
		  }else if(strncmp(line+11,"auth_comp_id", 12) == 0 && line[23] != '_' ){
			attribloc->resnameloc = count;
		  }else if(strncmp(line+11,"auth_asym_id", 12) == 0 && line[23] != '_' ){
			attribloc->chainloc = count;
		  }else if(strncmp(line+11,"auth_seq_id", 11) == 0 && line[22] != '_' ){
			attribloc->residloc = count;
		  }else if(strncmp(line+11,"Cartn_x", 7) == 0 && line[18] != '_' ){
			attribloc->xloc = count;
		  }else if(strncmp(line+11,"Cartn_y", 7) == 0 && line[18] != '_' ){
			attribloc->yloc = count;
		  }else if(strncmp(line+11,"Cartn_z", 7) == 0 && line[18] != '_' ){
			attribloc->zloc = count;
		  }else if(strncmp(line+11,"occupancy", 9) == 0 && line[20] != '_' ){
			attribloc->occuloc = count;
		  }else if(strncmp(line+11,"pdbx_PDB_model_num", 18) == 0 && line[29] != '_' ){
			attribloc->modelloc = count;
		  }else if(strncmp(line+11,"B_iso_or_equiv", 14) == 0 && line[25] != '_' ){
			attribloc->bfactloc = count;
		  }else if(strncmp(line+11,"pdbx_PDB_ins_code", 17) == 0 && line[28] != '_' ){
			attribloc->insloc = count;
		  }else if(strncmp(line+11,"label_alt_id", 12) == 0 && line[23] != '_' ){
			attribloc->altloc = count;
		  }else{
			;
		  }
		  count ++;
	    }

      }else if(strcmp(type, "label") == 0){
	    while(fgets(line, 1024, fp) != NULL   &&   strncmp(line, "_atom_site.", 11) == 0){
		  if(strncmp(line+11,"id", 2) == 0 && line[13] != '_' ){
			attribloc->idloc = count;
		  }else if(strncmp(line+11,"type_symbol", 11) == 0 && line[22] != '_' ){
			attribloc->type_symbloc = count;
		  }else if(strncmp(line+11,"label_atom_id", 13) == 0 && line[24] != '_' ){
			attribloc->locloc = count;
		  }else if(strncmp(line+11,"label_comp_id", 13) == 0 && line[24] != '_' ){
			attribloc->resnameloc = count;
		  }else if(strncmp(line+11,"label_asym_id", 13) == 0 && line[24] != '_' ){
			attribloc->chainloc = count;
		  }else if(strncmp(line+11,"label_seq_id", 12) == 0 && line[23] != '_' ){
			attribloc->residloc = count;
		  }else if(strncmp(line+11,"Cartn_x", 7) == 0 && line[18] != '_' ){
			attribloc->xloc = count;
		  }else if(strncmp(line+11,"Cartn_y", 7) == 0 && line[18] != '_' ){
			attribloc->yloc = count;
		  }else if(strncmp(line+11,"Cartn_z", 7) == 0 && line[18] != '_' ){
			attribloc->zloc = count;
		  }else if(strncmp(line+11,"occupancy", 9) == 0 && line[20] != '_' ){
			attribloc->occuloc = count;
		  }else if(strncmp(line+11,"pdbx_PDB_model_num", 18) == 0 && line[29] != '_' ){
			attribloc->modelloc = count;
		  }else if(strncmp(line+11,"B_iso_or_equiv", 14) == 0 && line[25] != '_' ){
			attribloc->bfactloc = count;
		  }else if(strncmp(line+11,"pdbx_PDB_ins_code", 17) == 0 && line[28] != '_' ){
			attribloc->insloc = count;
		  }else if(strncmp(line+11,"label_alt_id", 12) == 0 && line[23] != '_' ){
			attribloc->altloc = count;
		  }else{
			;
		  }
		  count ++;
	    }

      }else{
	    /* Exception Handling */ 
	    fprintf(stderr, "Error...  Please supply either \"auth\" or \"label\" for cif reading.\n");
	    exit(EXIT_FAILURE);

      }

}
static int guess_size_cif(FILE* fp, char* line, enum polymer_type polytype){
      int max_size = 0;
      int buffer = 2000;
      int proatm = 0;
      int nucatm = 0;
      int hohatm = 0;
      int ligatm = 0;
      int total  = 0;
      while(fgets(line, 1024, fp) != NULL    &&    strncmp(line, "_atom_site.", 11) != 0){
	    if(strncmp(line, "_refine_hist.pdbx_number_atoms_protein", 38) == 0){
		  proatm = atoi(line+38);
	    }else if(strncmp(line, "_refine_hist.pdbx_number_atoms_nucleic_acid", 43) == 0){
		  nucatm = atoi(line+43);
	    }else if(strncmp(line, "_refine_hist.pdbx_number_atoms_ligand", 37) == 0){
		  ligatm = atoi(line+37);
	    }else if(strncmp(line, "_refine_hist.number_atoms_solvent", 33) == 0){
		  hohatm = atoi(line+33);
	    }else if(strncmp(line, "_refine_hist.number_atoms_total", 31) == 0){
		  total = atoi(line+31);
	    }
      }
      if(polytype == NUC_TYPE){
	    max_size = nucatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == PRO_TYPE){
	    max_size = proatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == SOLVENT_TYPE){
	    max_size = nucatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == METAL_TYPE){
	    max_size = ligatm;  // buffer is extra space for some extra atoms.
      }else if(polytype == ALL_TYPE){
	    max_size = total;  // buffer is extra space for some extra atoms.
      }else{
	    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Wrong molecule type selected.\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }


      if(max_size == 0){
	    //if no such record present in the cif file. 
	    max_size = 10000;
      }else{
	    max_size += buffer;
      }
      return max_size;
}
static struct atom parse_cif_line_to_atom(char* line, struct cif_attrloc* attrloc){
      int indx = 0;
      char tokary[32][32];
      char sep[10] = "\t \"\n";
      char* token;
      token = strtok(line, sep);
      struct atom atom;


      while(token != NULL){
	    strcpy(tokary[indx], token);
	    indx++;
	    token = strtok(NULL, sep);
      }

      atom.id = atol(tokary[attrloc->idloc]);
      if(strlen(tokary[attrloc->type_symbloc]) > 2){
	    printf("atom is=%d\n", atom.id);
      }
      strcpy(atom.symbol, tokary[attrloc->type_symbloc]);
      atom.resid = atol(tokary[attrloc->residloc]);
      if(strlen(tokary[attrloc->locloc]) > 4){
	    printf("atom is=%d\n", atom.id);
      }
      strcpy(atom.loc, tokary[attrloc->locloc]);
      atom.altloc = tokary[attrloc->altloc][0]== '.'  ? ' '  :   tokary[attrloc->altloc][0];
      if(strlen(tokary[attrloc->resnameloc]) > 3){
	    printf("atom is=%d\n", atom.id);
      }
      strcpy(atom.resname, tokary[attrloc->resnameloc]);
      if(strlen(tokary[attrloc->chainloc]) > 2){
	    printf("atom is=%d\n", atom.id);
      }
      strcpy(atom.chain, tokary[attrloc->chainloc]);
      if(strlen(tokary[attrloc->insloc]) > 4){
	    printf("atom is=%d\n", atom.id);
      }
      strcpy(atom.ins, tokary[attrloc->insloc]);
      atom.model = atoi(tokary[attrloc->modelloc]);
      if(strncmp(tokary[attrloc->grouploc], "HETATM", 6) == 0)
	    atom.type = 'H';
      else
	    atom.type = 'A';

      atom.center.x = atof(tokary[attrloc->xloc]);
      atom.center.y = atof(tokary[attrloc->yloc]);
      atom.center.z = atof(tokary[attrloc->zloc]);
      atom.occu = atof(tokary[attrloc->occuloc]);
      atom.bfact = atof(tokary[attrloc->bfactloc]);
      return atom;

}
void scancif(const char* ciffile, int (*pf)(char*), const char* chain, const char* modelno, struct atom** atomary, int* size, enum polymer_type polytype, char* label_or_auth, char occurule){

      FILE	*fp;										/* input-file pointer */
      fp	= fopen( ciffile, "r" );
      if ( fp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			ciffile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }

      int model;
      if(modelno == NULL){
	    model = -1;
      }else{
	    model = atoi(modelno);
      }

      char line[1024];
      char type[10];
      if(label_or_auth == NULL){
	    strcpy(type , "auth");
      }else{
	    strcpy(type, label_or_auth);
      }
      struct cif_attrloc attribloc;
      int max_size = guess_size_cif(fp, line, polytype);
      set_cif_attrib_loc(&attribloc, line, type, fp);
      int atmindx = 0;


      struct atom atom;
      struct atom* atomarr= NULL;
      //fprintf(stdout, "line=%s\n", line);
      struct atom buffer = parse_cif_line_to_atom(line, &attribloc);
      while(!  ((modelno == NULL || buffer.model == model) && ((chain == NULL)||(strcmp(chain, buffer.chain) == 0)) && pf(buffer.resname) == 1)){
	    if(fgets(line, 1024, fp) == NULL || line[0] == '#' ){
		  goto final;
	    }
	    buffer = parse_cif_line_to_atom(line, &attribloc);
      }
      struct atom current;
      atomarr = (struct atom*) malloc (max_size * sizeof(struct atom));
      while(fgets(line, 1024, fp) != NULL && line[0] != '#'){
	    if(strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0) continue;
	    atom = parse_cif_line_to_atom(line, &attribloc);
	    if(atom.model != buffer.model) break;
	    if(  (pf(atom.resname) == 1)     &&    ((chain == NULL) || (strcmp(chain,atom.chain) == 0))){
		  int status = pdb_multi_occupancy(&atom, &buffer, &current, occurule);
		  if(status == 1){
			//print_pdb_line(stdout, &current);
			insert_atom(&atomarr, &atmindx, &max_size, current);
		  }
	    }
      }
      if((pf(buffer.resname) == 1) && ((chain == NULL) || (strcmp(chain,buffer.chain)== 0))){
	    struct atom dummy;
	    dummy.occu = ' ';
	    int status = pdb_multi_occupancy(&dummy, &buffer, &current, occurule);
	    if(status == 1){
		  //print_pdb_line(stdout, &current);
		  insert_atom(&atomarr, &atmindx, &max_size, current);
	    }
      }
final:
      *size = atmindx;
      *atomary = atomarr; 
      if( fclose(fp) == EOF ) {			/* close input file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			ciffile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      return ;
}
void fname_split(char *path, char *basename, char *ext, char *filename) {
      char* extn = strrchr(filename,'.');
      if(extn == NULL || extn == filename){
	    fprintf(stderr, "Error...  in function %s(). Improper file extension encountered. (Line: %d, File:%s\n", __func__ , __LINE__, __FILE__);
	    exit(EXIT_FAILURE);
      }
      strcpy(ext, extn);
      char* sep = strrchr(filename, '/');

      if(sep == NULL){
	    path[0] ='\0' ;

	    strncpy(basename, filename, (size_t)(extn - filename));
	    basename[(size_t)(extn - filename)] = '\0';
      }else{
	    strncpy(path, filename, (size_t)(sep - filename)+1);
	    path[(size_t)(sep - filename) + 1] = '\0';
	    strncpy(basename, sep + 1, (size_t)(extn-sep) - 1);
	    basename[(size_t)(extn - sep) - 1] = '\0';
      }
}
void fname_join(char *filename, const char *path, const char *basename, const char *ext) {
      strcpy(filename, path);
      strcat(filename, basename);
      strcat(filename, ext);
}
//void scancif1(const char* ciffile, int (*pf)(char*), const char* chain, const char* modelno, struct atom** atomary, int* size, enum polymer_type polytype){
//
//      FILE* fp = fopen(ciffile, "r");
//      if(fp == NULL){    /* Exception Handling */ 
//	    fprintf(stderr, "Error in function %s(). (File: %s, Line %d)... cannot open cif file.\n", __func__, __FILE__, __LINE__);
//	    exit(EXIT_FAILURE);
//      }
//      char line[1024];
////      int grouploc=-1;  // group_PDB
////      int type_symbloc = -1;
////
////      int idloc = -1;   // id
////      int locloc = -1;  // label_atom_id
////      int altloc = -1;  // label_alt_id
////      int resnameloc = -1;  // label_comp_id
////      int chainloc = -1; // auth_asym_id
////      int residloc = -1; // auth_seq_id
////      int modelloc = -1; // pdbx_PDB_model_num
////      int occuloc  = -1; //occupancy
////      int insloc   = -1;
////      int bfactloc = -1;
////      int xloc = -1;
////      int yloc = -1;
////      int zloc = -1;
////      int count = -1;
//      int max_size = guess_size_cif(fp, polytype);
//      int atomindex = 0;
//
//      struct atom* atomarr = (struct atom*) malloc (max_size * sizeof(struct atom));
////      while(fgets(line, 1024, fp) != NULL){
////	    if(strncmp(line, "_atom_site.", 11) != 0 && count == -1) continue;
////	    if(strncmp(line, "_atom_site.", 11) != 0 && count != -1) break;
//////	    count ++;
//////	    if(strncmp(line,"_atom_site.group_PDB", 20) == 0 && line[20] != '_' ){
//////		  grouploc = count;
//////	    }else if(strncmp(line,"_atom_site.id", 13) == 0 && line[13] != '_' ){
//////		  idloc = count;
//////	    }else if(strncmp(line,"_atom_site.type_symbol", 22) == 0 && line[22] != '_' ){
//////		  type_symbloc = count;
//////	    }else if(strncmp(line,"_atom_site.label_atom_id", 24) == 0 && line[24] != '_' ){
//////		  locloc = count;
//////	    }else if(strncmp(line,"_atom_site.label_comp_id", 24) == 0 && line[24] != '_' ){
//////		  resnameloc = count;
//////	    }else if(strncmp(line,"_atom_site.auth_asym_id", 23) == 0 && line[23] != '_' ){
//////		  chainloc = count;
//////	    }else if(strncmp(line,"_atom_site.auth_seq_id", 22) == 0 && line[22] != '_' ){
//////		  residloc = count;
//////	    }else if(strncmp(line,"_atom_site.Cartn_x", 18) == 0 && line[18] != '_' ){
//////		  xloc = count;
//////	    }else if(strncmp(line,"_atom_site.Cartn_y", 18) == 0 && line[18] != '_' ){
//////		  yloc = count;
//////	    }else if(strncmp(line,"_atom_site.Cartn_z", 18) == 0 && line[18] != '_' ){
//////		  zloc = count;
//////	    }else if(strncmp(line,"_atom_site.occupancy", 20) == 0 && line[20] != '_' ){
//////		  occuloc = count;
//////	    }else if(strncmp(line,"_atom_site.pdbx_PDB_model_num", 29) == 0 && line[29] != '_' ){
//////		  modelloc = count;
//////	    }else if(strncmp(line,"_atom_site.B_iso_or_equiv", 25) == 0 && line[25] != '_' ){
//////		  bfactloc = count;
//////	    }else if(strncmp(line,"_atom_site.pdbx_PDB_ins_code", 28) == 0 && line[28] != '_' ){
//////		  insloc = count;
//////	    }else if(strncmp(line,"_atom_site.label_alt_id", 23) == 0 && line[23] != '_' ){
//////		  altloc = count;
//////	    }else{
//////		  ;
//////	    }
//////
////      }
//      struct atom atom;
//      struct atom buffer = parse_cif_line_to_atom(line, attrloc);
//      struct atom current;
//      while(line[0] != '#'){
//	    int indx=0;
//	    token = strtok(line, sep);
//	    while(token != NULL){
//		  strcpy(tokary[indx], token);
//		  indx++;
//		  token = strtok(NULL, sep);
//	    }
//
//	    if(tokary[type_symbloc][0] == 'H'){ // We dont consider H
//		  if(fgets(line, 1024, fp) == NULL) break;
//		  continue;
//	    }
//	    if(pf(tokary[resnameloc]) < 0 ){
//		  if(fgets(line, 1024, fp) == NULL) break;
//		  continue;
//	    }
//
//	    long resid = atol(tokary[residloc]);
//	    char chn[10];
//	    strcpy(chn, tokary[chainloc]);
//	    if(prevresid != resid || strcmp(prevchain, chn) != 0){
//		  resindx ++;
//		  prevresid = resid;
//		  strcpy(prevchain, chn);
//	    }
//
//	    struct atom atom;
//	    char atomhetatom[20];
//	    strcpy(atomhetatom, tokary[grouploc]);
//
//	    atom.id = atol(tokary[idloc]);
//	    strcpy(atom.symbol, tokary[type_symbloc]);
//	    atom.resid = atol(tokary[residloc]);
//	    strcpy(atom.loc, tokary[locloc]);
//	    atom.altloc = tokary[altloc][0]== '.'? ' ': tokary[altloc][0];
//	    strcpy(atom.resname, tokary[resnameloc]);
//	    strcpy(atom.chain, tokary[chainloc]);
//	    strcpy(atom.ins, tokary[insloc]);
//	    if(strncmp(tokary[grouploc], "HETATM", 5) == 0)
//		  atom.type = 'H';
//	    else
//		  atom.type = 'A';
//
//	    atom.center.x = atof(tokary[xloc]);
//	    atom.center.y = atof(tokary[yloc]);
//	    atom.center.z = atof(tokary[zloc]);
//	    atom.occu = atof(tokary[occuloc]);
//	    atom.bfact = atof(tokary[bfactloc]);
//	    if(atomindex == max_arr_size){
//		  max_arr_size *= 2;
//		  atomarr = (struct atom*) realloc(atomarr, max_arr_size * sizeof(struct atom));
//	    }
//	    atomarr[atomindex] = atom;
//	    atomindex ++;
//	    print_pdb_line(stdout, &atom);
//	    if(fgets(line, 1024, fp) == NULL) break;
//	    //printf("%s",line);
//      }
//      *size = atomindex;
//      *atomary = atomarr; 
//      fclose(fp);
//      //if(this->size>0)
//      //printf("Moltest1=%s %lf\n", this->residue[0].atom[0].resname, residue[0].atom[0].x);
//      //else{
//      //   printf("No Mol found\n");
//      //}
//
//}
