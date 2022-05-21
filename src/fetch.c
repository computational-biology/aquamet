/*
 * =====================================================================================
 *
 *       Filename:  fetch.c
 *
 *    Description:  Fetching from RCSB.
 *
 *        Version:  1.0
 *        Created:  Thursday 19 May 2022 05:22:57  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "fetch.h"


int check_internet_conn(void) {
    char cmd[] = "curl -Is http://www.rcsb.com|head -1";    
    FILE* fp = NULL;
    char buf[FETCH_BUFF_SIZE];

    if ((fp = popen(cmd, "r")) == NULL) {
        fprintf(stderr, "Error opening pipe.\n");
        exit(EXIT_FAILURE);
    }
    
    
    if(fgets(buf, FETCH_BUFF_SIZE, fp) != NULL) {
        // Do whatever you want here...
        //fprintf(stdout, "Status: %s", buf);
        fprintf(stdout, "Connection established...\n");
    }else{
        fprintf(stderr, "Connection couldn't be established.\nPlease check internet connection!\n");
        exit(EXIT_FAILURE);
    }
    
    

    if (pclose(fp)) {
        printf("Command not found or exited with error status\n");
        exit(EXIT_FAILURE);
    }

    return 0;
}
void fetch_rcsb(char* accn, char* ext)
{
      FILE* fp = stdout;
      fprintf(fp, "Preparing file dowinload......\n");
      fprintf(fp, "Checking internet connection......\n");
      if(check_internet_conn() == 0){
          fprintf(fp, "Connection Detected...\n");
          fprintf(fp, "Download started... accn:%s%s (This may take several minutes. Please wait).\n", accn,ext);
          char cmd[512];
          sprintf(cmd,"curl -# -s -f https://files.rcsb.org/download/%s%s -o %s%s", accn, ext, accn, ext);
          /*char url[512];
          
          sprintf(url,"https://files.rcsb.org/download/%s.cif", accn);
          int pid = fork();
          if(pid == 0){
              execlp("curl", "-s", "-f", "https://files.rcsb.org/download/1ehz.cif", "-o", accn, NULL);
          }else{
              wait(&pid);
          }*/
          system(cmd);

          fprintf(fp, "Download finished...\n");
          fprintf(fp, "Starting Computations...\n");
          
      }
      
      
}
