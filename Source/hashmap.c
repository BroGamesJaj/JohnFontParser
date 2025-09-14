
#include "hashmap.h"

bool CreateHashMap(Hashmap* map,const int size) {
    map->size = size;
    map->buckets = (Pair**)calloc(size, sizeof(Pair*));
    
    if(!map->buckets) return true;
    return false;
}

int HashFunction(int* result, char* key, const int size){
    if(!key) return true;
    if(size == 0) return true;
    int hash = 0; // unneeded, but more understandable
    while (*key) {
        hash = (hash * 17) + *key++; // 17 is just a rando prime
    }
    *result = hash % size;
    return false;
}

bool InsertPair(Hashmap* map, char* key, uint32_t* value){
    if(!map || !map->buckets) return true;
    if(!key) return true;
    if(!value) return true;
    int index;
    if(HashFunction(&index, key, map->size)) return true;
    Pair* n_pair = calloc(1, sizeof(Pair));
    if(!n_pair) return true;
    n_pair->key = strdup(key);
    n_pair->value = *value;
    n_pair->next = map->buckets[index];
    map->buckets[index] = n_pair;
    return false;
}

bool GetValue(Hashmap* map, char* key, uint32_t* value){
    if(!map || !map->buckets) return true;
    if(!key) return true;
    int index;
    if(HashFunction(&index, key, map->size)) return true;
    if(!map->buckets[index]) return true;
    Pair* pair;
    pair = map->buckets[index];
    while(pair){
        if(strcmp(pair->key, key) == 0) {
            *value = pair->value;
            return false;
        }
        pair = pair->next;
    }
    return true;
}

bool DestroyHashMap(Hashmap* map){
    if(!map || !map->buckets) return true;

    for (int i = 0; i < map->size; i++) {
        Pair* curr = map->buckets[i];
        while (curr) {
            Pair* tmp = curr;
            curr = curr->next;
            free(tmp->key);
            free(tmp);
        }
        map->buckets[i] = NULL;
    }

    free(map->buckets);
    map->buckets = NULL;
    map->size = 0;

    return false;
}

void PrintHashMap(Hashmap* map) {
    if (!map || !map->buckets) {
        printf("HashMap is empty or uninitialized.\n");
        return;
    }

    for (int i = 0; i < map->size; i++) {
        Pair* curr = map->buckets[i];
        if (curr) {
            printf("Bucket %d:\n", i);
            while (curr) {
                printf("  Key: %s, Value: %u\n", curr->key, curr->value);
                curr = curr->next;
            }
        }
    }
}