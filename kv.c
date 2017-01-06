/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kv.h"


// take data in stream f and put in a ptr to an array of ptr with size count
// adapted from assignment 49
void readfile(FILE * f, char *** str_array, size_t *count) 
{
  char * str = NULL;
  size_t size = 0;
  while( getline(&str, &size, f) >= 0 ) 
  {
      *str_array = realloc(*str_array, (*count+1)*sizeof(**str_array));
      if(*str != '\n') //skip empty line; TODO: write a better treatment!
      { 
         (*str_array)[*count] = str;
         (*count)++;
      }
      else // str points to '\n', which needs to be freed before proceeding!
         free(str);
      str = NULL;
  }
  free(str); // the last step in while loop gives str an empty string "", so feel free to free it.
}


kvpair_t * create_kvpair(char * str)
{
   kvpair_t * kvpair = malloc(sizeof(*kvpair));
   kvpair->key = NULL;
   kvpair->value = NULL;
   
   char * equal_sign1 = strchr(str, '=');
   if(!equal_sign1) 
   {
      fprintf(stderr, "Equal sign is not found in the intput. Abort!\n");
      exit(EXIT_FAILURE);
   }

   char * line_end = strchr(str, '\n'); //str is read using getline, so the string is ended with "\n\0"
   if(!line_end)
   {
      fprintf(stderr, "Line delimiter is not found in the input. Abort!\n");
      exit(EXIT_FAILURE);
   }
   line_end--;

   //get key
   kvpair->key = realloc(kvpair->key, (equal_sign1-str+1)*sizeof(char) );
   strncpy(kvpair->key, str, equal_sign1-str);
   kvpair->key[equal_sign1-str] = '\0';
   
   //get value
   kvpair->value = realloc(kvpair->value, (line_end-equal_sign1+1)*sizeof(char) );
   strncpy(kvpair->value, equal_sign1+1, line_end-equal_sign1);
   kvpair->value[line_end-equal_sign1] = '\0';

   return kvpair;
}


kvarray_t * readKVs(const char * fname) {
  kvarray_t * kvarray = malloc(sizeof(*kvarray));
  kvarray->kvpair = NULL;
  kvarray->kvpair_len = 0;
  
  //open file
  FILE * kv_input_file = fopen(fname, "r");
  if(!kv_input_file)
  {
     perror("Could not open key-value input file.\n");
     exit(EXIT_FAILURE);
  }

  //read raw strings from file
  char ** raw_data = NULL; 
  size_t size = 0;
  readfile(kv_input_file, &raw_data, &size);

  //create kvpair_t object
  for(int i=0; i<size; i++)
  {
     kvarray->kvpair = realloc(kvarray->kvpair, (kvarray->kvpair_len+1)*sizeof(*(kvarray->kvpair)) );
     kvarray->kvpair[i] = create_kvpair(raw_data[i]);
     kvarray->kvpair_len++;
  }

  //quick check
  if(kvarray->kvpair_len != size) 
  {
     fprintf(stderr, "Something is wrong. Abort!\n");
     exit(EXIT_FAILURE);
  }   

  //close file
  if( fclose(kv_input_file) )
  {
     perror("Error when closing the key-value input file.\n");
     exit(EXIT_FAILURE);
  }

  //free raw data
  for(int i=0; i<size; i++)
     free(raw_data[i]);
  free(raw_data);

  return kvarray;
}


void freeKVs(kvarray_t * pairs) {
  for(int i=0; i<pairs->kvpair_len; i++)
  {
     free(pairs->kvpair[i]->key);
     free(pairs->kvpair[i]->value);
     free(pairs->kvpair[i]);
  }
  free(pairs->kvpair);
  free(pairs);
}


void printKVs(kvarray_t * pairs) {
  for(int i=0; i<pairs->kvpair_len; i++)
     printf("key = '%s' value = '%s'\n", pairs->kvpair[i]->key, pairs->kvpair[i]->value);
}


char * lookupValue(kvarray_t * pairs, const char * key) {
  //const char * str1=NULL;
  //const char * str2=NULL;
  int key_index = 0;
  int found_key = 0;

  for(int i=0; i<pairs->kvpair_len; i++)
  {
     if(strcmp(pairs->kvpair[i]->key, key)==0)
     //str1 = strstr(pairs->kvpair[i]->key, key);
     //str2 = strstr(pairs->kvpair[i]->key+1, key);
     //if( str1 && !str2 )
     {
        key_index = i;
        found_key = 1;
        break;
     }
  }
  
  if(!found_key)
  {
      printf("%s: input parameter %s is not found.\n", __func__, key);
      return NULL;
//      exit(EXIT_FAILURE);
  }
  else
      return pairs->kvpair[key_index]->value;
}
