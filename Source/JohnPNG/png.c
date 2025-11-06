#include "png.h"

uint32_t crcTable[256];
char crcTableComputed = 0;

void MakeCRCTable(void)
{
  uint32_t c;
  int n, k;

  for (n = 0; n < 256; n++) {
    c = (uint32_t) n;
    for (k = 0; k < 8; k++) {
      if (c & 1)
        c = 0xedb88320L ^ (c >> 1);
      else
        c = c >> 1;
    }
    crcTable[n] = c;
  }
  crcTableComputed = 1;
}

uint32_t UpdateCRC(uint32_t crc, uint8_t *buf,
                         int len)
{
  unsigned long c = crc;
  int n;

  if (!crcTableComputed)
    MakeCRCTable();
  for (n = 0; n < len; n++) {
    c = crcTable[(c ^ buf[n]) & 0xff] ^ (c >> 8);
  }
  return c;
}

int Adduint32ToBuffer(uint8_t* buffer, int* pos, uint32_t data){
    buffer[(*pos)++] = (data >> 24) & 0xFF;
    buffer[(*pos)++] = (data >> 16) & 0xFF;
    buffer[(*pos)++] = (data >> 8) & 0xFF;
    buffer[(*pos)++] = data & 0xFF;
    return 0;
}


int CreateIHDR(FILE* outputFile, uint32_t width, uint32_t height, uint8_t bitDepth, uint8_t colorType, uint8_t compressionMethod, uint8_t filterMethod, uint8_t interlaceMethod){
    uint8_t buffer[22];
    int pos = 0;
    uint32_t dataFieldLength = 0;
    fwrite(&dataFieldLength,sizeof(dataFieldLength),1,outputFile);
    uint8_t IHDRMagicNumbers[4] = {0x49, 0x48, 0x44, 0x52};
    memcpy(buffer,IHDRMagicNumbers, sizeof(IHDRMagicNumbers));
    pos += 4;

    Adduint32ToBuffer(buffer,&pos,width);
    Adduint32ToBuffer(buffer,&pos,height);

    buffer[pos++] = bitDepth;
    buffer[pos++] = colorType;
    buffer[pos++] = compressionMethod;
    buffer[pos++] = filterMethod;
    buffer[pos++] = interlaceMethod;

    uint32_t crc = UpdateCRC(0xFFFFFFFF,buffer,pos);
    Adduint32ToBuffer(buffer,&pos,crc);
    fwrite(buffer, pos, 1, outputFile);
    return 0; 
}

int CreateIDAT(FILE* outputFile){
    
    uint8_t* buffer;
    int pos = 0;

    uint8_t IDATMagicNumbers[4] = {0x49, 0x44, 0x41, 0x54};
    memcpy(buffer, IDATMagicNumbers, 4);
    pos += 4;


}

int Test(){
    FILE* tf = fopen("red.png","wb");

    uint8_t pngMagicNumbers[8] = {0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A };
    fwrite(pngMagicNumbers,sizeof(pngMagicNumbers),1,tf);

    CreateIHDR(tf, 1, 1, 8, 2, 0, 0, 0);
    uint8_t* buffer = malloc(20);
    int pos = 0;
    //IDAT

    

    fclose(tf);
    free(buffer);
    return 0;
}
int main(int argc, char** argv){
    if(argc > 1 && !strcmp(argv[1], "test"))
        Test();

    return 0;
}
