#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef HASHMAP_H
#define HASHMAP_H


typedef struct Pair {
    char* key;
    uint32_t value;
    struct Pair* next;
} Pair;

typedef struct Hashmap {
    int size;
    int capacity;
    Pair** buckets;
} Hashmap;

bool CreateHashMap(Hashmap* map, int size);
int HashFunction(int* result, char* key, const int size);
bool InsertPair(Hashmap* map, char* key, uint32_t* value);
bool GetValue(Hashmap* map, char* key, uint32_t* value);
bool DestroyHashMap(Hashmap* map);
void PrintHashMap(Hashmap* map);


#endif