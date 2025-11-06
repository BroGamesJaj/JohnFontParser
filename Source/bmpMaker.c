#include "bmpMaker.h"
#pragma pack()

struct header {
    uint32_t bitmap_size;
    uint32_t reserved;
    uint32_t img_offset;
};

struct bitmapInfoHeader {
    uint32_t info_header_size;
    int32_t bitmap_width;
    int32_t bitmap_height;
    uint16_t num_of_color_planes;
    uint16_t bits_per_pixel;
    uint32_t compression_method;
    uint32_t image_size;
    int32_t horizontal_resolution;
    int32_t vertical_resolution;
    uint32_t num_of_colors_pallette;
    uint32_t num_of_important_colors;
};

struct bitmap {
    struct header header;
    struct bitmapInfoHeader bitmapInfoHeader;
};

int TestEndian(){
    uint32_t n = 1;

    if(*(char *)&n == 1 && *(char*)&n+'0' == '1'){
        return 1; //Little endian
    } else {
        return 0; //Big endian
    }
}

uint16_t ChangeEndianUInt16_t(uint16_t num){
    num = (num >> 8) | (num << 8);
    return num;
}

uint32_t ChangeEndianUInt32_t(uint32_t num){
    num = ((num >> 24) & 0x000000FF) |
            ((num >> 8)  & 0x0000FF00) |
            ((num << 8)  & 0x00FF0000) |
            ((num << 24) & 0xFF000000);
    return num;
}

int32_t ChangeEndianInt32_t(int32_t num){
    num =  (int32_t)(uint32_t)((num >> 24) & 0x000000FF) |
            ((num >> 8)  & 0x0000FF00) |
            ((num << 8)  & 0x00FF0000) |
            ((num << 24) & 0xFF000000);
    return num;
}

int MakeBMP(FILE* bmp, uint8_t* pixel_data, int32_t width, const int32_t height, const uint16_t bits_per_pixel, const uint32_t compression_method,
    const int32_t horizontal_resolution, const int32_t vertical_resolution){

    if(!TestEndian()) return 0;

    struct bitmap *bitmap = (struct bitmap*)calloc(1,sizeof(struct bitmap));
    
    uint8_t signature[2];
    signature[0] = 'B';
    signature[1] = 'M';
    fwrite(signature, 1, sizeof(signature), bmp);

    int32_t t = (int32_t)((width*bits_per_pixel+31)/32*4*abs(height) + sizeof(struct bitmap) + sizeof(signature));
    bitmap->header.bitmap_size = t;
    bitmap->header.img_offset = sizeof(struct bitmap) + sizeof(signature);

    bitmap->bitmapInfoHeader.info_header_size = sizeof(struct bitmapInfoHeader);
    bitmap->bitmapInfoHeader.bitmap_width = width;
    bitmap->bitmapInfoHeader.bitmap_height = height;
    bitmap->bitmapInfoHeader.num_of_color_planes = 1;
    bitmap->bitmapInfoHeader.bits_per_pixel = bits_per_pixel;
    bitmap->bitmapInfoHeader.compression_method = compression_method;
    bitmap->bitmapInfoHeader.image_size = (int32_t)ceil((double)(width*bits_per_pixel)/32)*4*abs(height);
    bitmap->bitmapInfoHeader.horizontal_resolution = horizontal_resolution;
    bitmap->bitmapInfoHeader.vertical_resolution = vertical_resolution;
    bitmap->bitmapInfoHeader.num_of_colors_pallette = 0;

    fwrite(bitmap, 1, sizeof(struct bitmap), bmp);

    int32_t used_bytes = ceil((double)width*bits_per_pixel/8);
    int32_t padded_row_bytes = ceil((double)(width*bits_per_pixel)/32)*4;
    int32_t pad = padded_row_bytes - used_bytes;

    uint32_t padded_pixel_data_size = padded_row_bytes * abs(height);

    uint8_t* padded_pixel_data = malloc(padded_pixel_data_size);

    if(padded_pixel_data == NULL) return 0;

    for(int i = 0; i < height; i++){
        memmove(padded_pixel_data + (padded_row_bytes * i),pixel_data + (used_bytes * i), used_bytes);
        memset(padded_pixel_data + (padded_row_bytes * i)+used_bytes, 0, pad);
    }
   
    fwrite(padded_pixel_data, 1, padded_pixel_data_size, bmp);

    free(padded_pixel_data);
    free(bitmap);
    return 1;
}

