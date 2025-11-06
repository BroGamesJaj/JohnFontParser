1. PNG magic number <br>
89 50 4E 47 0D 0A 1A 0A

2. Chunks

basic specs: <br>
- Length: uint32, only the chunk data. <br>
- chunk type: 4 bytes -> binary value. <br>
    - 4 bits (bit 5 in every byte (ie: 32)) are chunk properties.
    - ancillary bit: 0 = critical, 1 = ancillary.
    - private bit: 0 = public, 1 = private.
    - reserved bit: 0 = currently, 1 = not good.
    - safe-to-copy bit: 0 = unsafe to copy, 1 = safe to copy.
- chunk data: size depends on chunk type, can be 0. <br>
- CRC: calculated on the chunk excluding length.

chunk layout: 
- IHDR
- PLTE
- IDAT - one or multiple
- IEND

critical chunks:

- IHDR
    - chunk type: 49 48 44 52
    - chunk data:
        - Width 	4 bytes - uint
        - Height 	4 bytes - uint
        - Bit depth 	1 byte
        - Color type 	1 byte
        - Compression method 	1 byte
        - Filter method 	1 byte
        - Interlace method 	1 byte

- PLTE
    - required for color type 3, may appear for 2,6 and it shall not appear for 0,4.
    - chunk type: 50 4C 54 45
    - chunk data:
        - Red   1 byte
        - Green 1 byte
        - Blue  1 byte

- IDAT
    - multiple IDAT need to appear consecutively.
    - chunk type: 49 44 41 54
    - chunk data:
        - Image data

- IEND
    - chunk type: 49 45 4E 44 
