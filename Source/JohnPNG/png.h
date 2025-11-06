

#ifndef _PNG_H_
#define _PNG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef enum ColourType {
    Type_Greyscale = 0,
    Type_RGB = 2,
    Type_Indexedcolour = 3,
    Type_GreyscaleAlpha = 4,
    Type_RGBA = 6
} ColourType;

typedef struct JohnPNG {
    uint32_t width;
    uint32_t height;
    uint32_t capacity;
    ColourType type;
    uint32_t crc; // CRC32 chunk validation
} JohnPNG;

#ifdef __cplusplus
}
#endif
#endif
