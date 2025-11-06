#ifndef BMPMAKER_H
#define BMPMAKER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <memory.h>
#include <math.h>

int MakeBMP(FILE* bmp, uint8_t* pixel_data, const int32_t width, const int32_t height, const uint16_t bits_per_pixel, const uint32_t compression_method, const int32_t horizontal_resolution, const int32_t vertical_resolution);

#endif
