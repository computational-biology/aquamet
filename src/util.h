//
// Created by parthajit on 11/7/20.
//

#ifndef CPPMET_UTIL_H
#define CPPMET_UTIL_H
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>

#include "biodefs.h"

int find_amino(char* amino);

int find_nucleic(char* nucleic);
int is_HOH(char* res);

char get_res_loc(struct atom* atom, char restype);

void file_name_split(char* file_path, char* file_name, char* file_ext, char* src_file);

void file_name_join(char* joined_file_name, const char* path, const char* file_name, const char* ext);
void now(char output[]);
void today(char output[]);
#endif //CPPMET_UTIL_H
