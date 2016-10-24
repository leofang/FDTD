#ifndef __KV_H__
#define __KV_H__

// This struct is borrowed from the assignment (52_kvs) of ECE551
// that I took (credit: Prof. Andrew Hilton). I wrote the entire 
// body of the given prototypes in kv.c. This struct reads the 
// input parameters given as key-value pair (such as "nx=500") 
// and then stores the information.

struct _kvpair_t { 
  char * key;
  char * value;

};
typedef struct _kvpair_t kvpair_t;

struct _kvarray_t { 
  kvpair_t ** kvpair;
  size_t kvpair_len;

};
typedef struct _kvarray_t kvarray_t;


kvarray_t * readKVs(const char * fname);

void freeKVs(kvarray_t * pairs);

void printKVs(kvarray_t * pairs);

char * lookupValue(kvarray_t * pairs, const char * key);

#endif
