#include <stdio.h>
#include <math.h>
#include "hashmap.h"

// if 1 all 16bit or larger integers get swapped to Big-Endian
#define CONVERT_TO_B 1 
#define FWord int16_t
#define F2DOT14 uint16_t

// 0 if anything other than apple or windows ei.: unicode platform

#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WINDOWS__)
static int platform = 3;
#else 
static int platform = 0;
#endif


typedef struct Bounds {
    FWord xMin, yMin, xMax, yMax;
} Bounds;

typedef struct Point {
    int16_t x;
    int16_t y;
    bool isOnCurve;
} Point;

typedef struct Glyph Glyph;

typedef struct Glyph {
    uint16_t glyphIndex;
    int16_t numOfContours;
    uint16_t* endPtsOfContours;
    uint16_t instructionLength; // not yet needed
    //uint8_t* instructions; not yet needed, some small font size hint bullshit.
    uint8_t* flags;
    Bounds bounds;
    uint16_t numOfPoints;
    Point* points;
} Glyph;


typedef struct GlyphIndexes {
    int offsetSize;
    uint32_t* offsets;
} GlyphIndexes;

typedef struct EncodingRecord {
    uint16_t platformID;
    uint16_t encodingID;
    uint32_t subtableOffset;
} EncodingRecord;

typedef struct Font {
    Hashmap tables;
    uint16_t numOfTables;

    GlyphIndexes glyphOffsets;
    uint16_t numOfGlyphs;

    Glyph* glyphs;

} Font;

void SkipBytes(int numberOfBytesToSkip, FILE* file){
    fseek(file, numberOfBytesToSkip, SEEK_CUR);
}

bool ReadUInt8(uint8_t* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    return false;
}

bool ReadInt8(int8_t* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    return false;
}

bool ReadUInt16(uint16_t* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    if(CONVERT_TO_B){
        *buffer = ((*buffer >> 8) | (*buffer << 8));
        return false;
    } else {
        return false;
    }
}

bool ReadInt16(int16_t* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    if (CONVERT_TO_B) {
        uint16_t temp = (uint16_t)*buffer;
        temp = (temp >> 8) | (temp << 8);
        *buffer = (int16_t)temp;
    }
    return false;
}

bool ReadUInt32(uint32_t* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    if(CONVERT_TO_B){
        *buffer = ((*buffer >> 24) & 0x000000FF) |
          ((*buffer >> 8)  & 0x0000FF00) |
          ((*buffer << 8)  & 0x00FF0000) |
          ((*buffer << 24) & 0xFF000000);

        return false;
    } else {
        return false;
    }
}

bool ReadInt32(int32_t* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    if (CONVERT_TO_B) {
        uint32_t temp = (uint32_t)*buffer;
        temp = ((temp >> 24) & 0x000000FF) |
               ((temp >> 8)  & 0x0000FF00) |
               ((temp << 8)  & 0x00FF0000) |
               ((temp << 24) & 0xFF000000);
        *buffer = (int32_t)temp;
    }
    return false;
}

bool ReadF2DOT14(F2DOT14* buffer, FILE* file) {
    *buffer = 0;
    size_t bytesRead = fread(buffer, sizeof(*buffer), 1, file);
    if (bytesRead != 1) {
        return true;
    }
    if (CONVERT_TO_B) {
        uint16_t temp = (uint16_t)*buffer;
        temp = (temp >> 8) | (temp << 8);
        *buffer = (int16_t)(temp * 16384.0f);
    }
    return false;
}

bool ReadTag(char* buffer, FILE* file) {
    size_t bytesRead = fread(buffer, 1, 4, file);
    if (bytesRead != 4) return true;

    return false;
}

bool JumpToTable(char* tableName, Hashmap* tables, FILE* file) {
    uint32_t tableOffset;
    GetValue(tables, tableName, &tableOffset);
    fseek(file, tableOffset, SEEK_SET);
    return false;
}

bool FlagBit(uint8_t flag, int index) {
    return ((1 << index) & flag) != 0;
}

bool FlagBit16(uint16_t flag, int index) {
    return ((1 << index) & flag) != 0;
}

bool AddImpliedPoints(Glyph* glyph, float scale) {
    int contourStart = 0;
    int newPointCount = 0;
    for (size_t i = 0; i < glyph->numOfContours; i++) {
        int endPointIndex = glyph->endPtsOfContours[i];
        int pIndex;
        for (pIndex = contourStart; pIndex <= endPointIndex; pIndex++) {
            if(glyph->points[pIndex].isOnCurve){
                break;
            }
        }
        printf("endInd %zu: %u\n", i,  glyph->endPtsOfContours[i]);
        for (size_t j = contourStart; j <= endPointIndex; j++) {
            int currIndex = j + pIndex - contourStart;
            int nextIndex = (j + pIndex + 1 - contourStart) % (endPointIndex + 1);
            if(nextIndex < currIndex) nextIndex += contourStart;
            
            Point curr = glyph->points[currIndex];
            Point next = glyph->points[nextIndex];
            printf("currIndex: %d, nextIndex: %d (%d,%d)\n", currIndex, nextIndex, curr.isOnCurve, next.isOnCurve);
            if(curr.isOnCurve == next.isOnCurve && j <= endPointIndex) {
                newPointCount++;
            }
        }
        contourStart = endPointIndex + 1;
    }
    printf("new Point Count: %d\n", newPointCount);
    size_t newPointNum = glyph->numOfPoints + newPointCount;
    Point* newPoints = realloc(glyph->points, sizeof(Point) * newPointNum);
    if (newPoints == NULL) {
        return true;
    }
    glyph->points = newPoints;

    contourStart = 0;
    int newIndex = glyph->numOfPoints;
    int insertCount = 0;

    for (size_t i = 0; i < glyph->numOfContours; i++) {
        printf("insertCount: %d\n", insertCount);
        glyph->endPtsOfContours[i] += insertCount;
        int endPointIndex = glyph->endPtsOfContours[i];
        int pIndex;
        for (pIndex = contourStart; pIndex < endPointIndex; pIndex++) {
            if(glyph->points[pIndex].isOnCurve){
                break;
            }
        }
        insertCount = 0;
        
        for (size_t j = contourStart; j <= endPointIndex; j++) {
            int currIndex = j + pIndex - contourStart;
            int nextIndex = (j + pIndex + 1 - contourStart) % (endPointIndex + 1);
            if(nextIndex < currIndex) nextIndex += contourStart;
            Point curr = glyph->points[currIndex];
            Point next = glyph->points[nextIndex];
            
            if(curr.isOnCurve == next.isOnCurve && j <= endPointIndex) {
                Point impliedPoint;
                impliedPoint.x = (curr.x + next.x) / 2.0f * scale;
                impliedPoint.y = (curr.y + next.y) / 2.0f * scale;
                impliedPoint.isOnCurve = !curr.isOnCurve;

                for (size_t k = glyph->numOfPoints - 1; k >= currIndex + 1; k--) {
                    if (k + 1 < newPointNum) {
                        glyph->points[k + 1] = glyph->points[k];
                    }
                }
                printf("newPoint: %u (%d, %d) prev: %d, next: %d\n", currIndex + 1, impliedPoint.x, impliedPoint.y, currIndex, nextIndex);
                glyph->points[currIndex + 1] = impliedPoint;
                insertCount++;
                glyph->numOfPoints++;
                
                endPointIndex++;
                
                newIndex++;
            }
            printf("endPointIndex: %d\n", endPointIndex);
        }
        glyph->endPtsOfContours[i] += insertCount;
        contourStart = glyph->endPtsOfContours[i] + 1;
    }
    for (size_t i = 0; i < glyph->numOfContours; i++) {
        printf("endInd %zu: %u\n", i,  glyph->endPtsOfContours[i]);
    }
    
    return false;
}

bool ReadCoordinates(Point* points, uint16_t numOfPoints, uint8_t* flags, FILE* file) {
    // x coordinates
    for (size_t i = 0; i < numOfPoints; i++) {
        points[i].x = (i != 0) ? points[i - 1].x : 0; // Start with the previous point's position because they are relative positions.
        if(FlagBit(flags[i], 0)) { // Is point on curve
            points[i].isOnCurve = true;
        } else {
            points[i].isOnCurve = false;
        }
        if(FlagBit(flags[i], 1)) { // See if 1 or 2 byte, if 1 byte then the sign is determined by bit 4.
            uint8_t tmpCoord;
            ReadUInt8(&tmpCoord, file); 
            int sign = FlagBit(flags[i], 4) ? 1 : -1; // Get the sign of coord.
            points[i].x += (int16_t)tmpCoord * sign;
        } else if (!FlagBit(flags[i], 4)) { // See if coord same as previous.
            int16_t tmpCoord = 0;
            ReadInt16(&tmpCoord, file);
            points[i].x += tmpCoord;
        }
    }

    // y coordinates
    // Same as x coords, just flags are on different bit.
    for (size_t i = 0; i < numOfPoints; i++) {
        points[i].y = (i != 0) ? points[i - 1].y : 0;
        if(FlagBit(flags[i], 2)){
            uint8_t tmpCoord;
            ReadUInt8(&tmpCoord, file);
            int sign = FlagBit(flags[i], 5) ? 1 : -1;
            points[i].y += (int16_t)tmpCoord * sign;
        } else if (!FlagBit(flags[i], 5)) {
            int16_t tmpCoord = 0;
            ReadInt16(&tmpCoord, file);
            points[i].y += tmpCoord;
        }
    }
    
    return false;
}

bool ReadSimpleGlyph(Glyph* glyph, FILE* file){

    ReadInt16(&glyph->bounds.xMin, file);
    ReadInt16(&glyph->bounds.yMin, file);
    ReadInt16(&glyph->bounds.xMax, file);
    ReadInt16(&glyph->bounds.yMax, file);
    if(glyph->numOfContours != 0) {

        glyph->endPtsOfContours = calloc(glyph->numOfContours, sizeof(uint16_t));
        for (size_t i = 0; i < glyph->numOfContours; i++) {
            ReadInt16(&glyph->endPtsOfContours[i], file);
        }
    } else {
        printf("number of Contours is 0\n");
        return false;
    }
    
    ReadUInt16(&glyph->instructionLength, file);
    SkipBytes(glyph->instructionLength, file);
    glyph->numOfPoints = glyph->endPtsOfContours[glyph->numOfContours - 1] + 1;
    glyph->flags = calloc(glyph->numOfPoints, sizeof(uint8_t));
    printf("number of points: %u\n", glyph->numOfPoints);
    for (size_t i = 0; i < glyph->numOfPoints; i++) {
        ReadUInt8(&glyph->flags[i], file);
        if(FlagBit(glyph->flags[i], 3)) {
            uint8_t flagToRepeat = glyph->flags[i];
            uint8_t numOfRepeat;
            ReadUInt8(&numOfRepeat, file);
            if (i + numOfRepeat < glyph->numOfPoints) {
                for (size_t j = 1; j <= numOfRepeat; j++) {
                    glyph->flags[i + j] = flagToRepeat;
                }
                i += numOfRepeat;
            }
        }
    }
    glyph->points = calloc(glyph->numOfPoints, sizeof(Point));
    ReadCoordinates(glyph->points, glyph->numOfPoints, glyph->flags, file);
    for (size_t i = 0; i < glyph->numOfContours; i++) {
        printf("endInd %zu: %u\n", i,  glyph->endPtsOfContours[i]);
    }
    AddImpliedPoints(glyph, 1);
    for (size_t i = 0; i < glyph->numOfContours; i++) {
        printf("endInd %zu: %u\n", i,  glyph->endPtsOfContours[i]);
    }
    printf("numOfPoints: %u\n", glyph->numOfPoints);
    return false;
}

//TODO: finish composite glyph flattening
bool ReadCompositeGlyph(Glyph* glyph, uint32_t currentOffset, FILE* file) {
    ReadInt16(&glyph->bounds.xMin, file);
    ReadInt16(&glyph->bounds.yMin, file);
    ReadInt16(&glyph->bounds.xMax, file);
    ReadInt16(&glyph->bounds.yMax, file);
    currentOffset += 16 * 4;
    uint16_t flags;
    do {
        ReadUInt16(&flags, file);
        uint16_t glyphIndex;
        ReadUInt16(&glyphIndex, file);
        currentOffset += 16 * 2;
        printf("glyph index: %u\n", glyph->glyphIndex);

        for (int i = 15; i >= 0; i--) {
            printf("%d", (flags >> i) & 1);
        }
        printf("\n");
        if (FlagBit16(flags, 0)) {
            FWord argument1, argument2;
            ReadInt16(&argument1, file);
            ReadInt16(&argument2, file);
        } else {
            uint8_t arg1,arg2;
            ReadInt8(&arg1, file);
            ReadInt8(&arg2, file);
        }
        if (FlagBit16(flags, 2) ) {
            F2DOT14  scale;
            ReadF2DOT14(&scale, file);
        } else if (FlagBit16(flags, 6) ) {
            F2DOT14  xscale, yscale;
            ReadF2DOT14(&xscale, file);
            ReadF2DOT14(&yscale, file);
        } else if (FlagBit16(flags, 7) ) {
            F2DOT14  xscale, scale01, scale10, yscale;
            ReadF2DOT14(&xscale, file);
            ReadF2DOT14(&scale01, file);
            ReadF2DOT14(&scale10, file);
            ReadF2DOT14(&yscale, file);
        }
    } while (FlagBit16(flags, 5));

    if (FlagBit16(flags, 8)){
        uint16_t numInstr;
        //uint8_t instr[numInstr]; optional af instructions
    }

    return false;
}

bool ReadComponentGlyph(Glyph* glyph, FILE* file) {
    return false;
};

bool ReadAnyGlyph(Font* font, uint16_t index, FILE* file) {

    printf("glyph %u at offset: %u\n",index, font->glyphOffsets.offsets[index]);
    fseek(file, font->glyphOffsets.offsets[index], SEEK_SET);
    int16_t numOfContours;
    ReadInt16(&numOfContours, file);
    font->glyphs[index].numOfContours = numOfContours;
    font->glyphs[index].glyphIndex = index;
    printf("number of contours: %d\n", font->glyphs[index].numOfContours);
    if(numOfContours < 0){
        printf("composite glyph\n");
        //ReadCompositeGlyph(&font->glyphs[index], font->glyphOffsets.offsets[index], file);
        font->glyphs[index].glyphIndex = index;
    } else {
        printf("simple glyph\n");
        ReadSimpleGlyph(&font->glyphs[index], file);
    }
    return false;
}

void DestroyGlyph(Glyph* glyph){
    if(glyph && glyph->endPtsOfContours) {
        free(glyph->endPtsOfContours);
    }
    if(glyph && glyph->flags) {
        free(glyph->flags);
    }
    if(glyph && glyph->points){
        free(glyph->points);
    }
}

void DestroyFont(Font* font){
    
    if(font->glyphs){
        for (size_t i = 0; i < font->numOfGlyphs; i++) {
            DestroyGlyph(&font->glyphs[i]);
        }
        free(font->glyphs);
    }
    DestroyHashMap(&font->tables);
    if(font->glyphOffsets.offsets) {
        free(font->glyphOffsets.offsets);
    }
}

bool GetOffsetOfAllGlyphs(Font* font, FILE* file){
    // Find the loca table glyph offset format
    JumpToTable("head", &font->tables, file);
    SkipBytes(50, file);
    int16_t indexToLocFormat;
    ReadInt16(&indexToLocFormat, file);
    uint32_t* tmpOffsets = calloc(font->numOfGlyphs, sizeof(uint32_t));
    int pos = 1;
    // Go to loca table
    JumpToTable("loca", &font->tables, file);
    // Read the first so its easier later.
    if(indexToLocFormat == 0){
        uint16_t tmp;
        ReadUInt16(&tmp, file); 
        tmpOffsets[0] = tmp;
    } else {
        ReadUInt32(&tmpOffsets[0], file);
    }
    for (size_t i = 1; i < font->numOfGlyphs + 1; i++) {
        if(indexToLocFormat == 0){
            uint16_t tmp;
            ReadUInt16(&tmp, file);
            tmpOffsets[pos] = tmp;
            if(tmpOffsets[pos] != tmpOffsets[pos - 1]){
                pos++;
            }
        } else {
            ReadUInt32(&tmpOffsets[pos], file);
            if(tmpOffsets[pos] != tmpOffsets[pos - 1]){
                pos++;
            }
        }
    }
    font->glyphOffsets.offsets = calloc(pos, sizeof(uint32_t));
    uint32_t glyftableOffset;
    
    GetValue(&font->tables, "glyf", &glyftableOffset);
    font->numOfGlyphs = pos;
    for (size_t i = 0; i < pos; i++) {
        font->glyphOffsets.offsets[i] = tmpOffsets[i] + glyftableOffset;
    }
    free(tmpOffsets);
    return false;
}

bool GetCharMapping(Hashmap* tables, FILE* file) {
    JumpToTable("cmap", tables, file);
    uint32_t tableOffset;
    GetValue(tables, "cmap", &tableOffset);
    // skip u16 version number
    SkipBytes(2, file);
    uint16_t numOfTables;
    ReadUInt16(&numOfTables, file);
    EncodingRecord* encodingRecords = calloc(numOfTables, sizeof(EncodingRecord));
    printf("number of char tables: %u\n", numOfTables);
    for (size_t i = 0; i < numOfTables; i++) {
        ReadUInt16(&encodingRecords[i].platformID, file);
        ReadUInt16(&encodingRecords[i].encodingID, file);
        ReadUInt32(&encodingRecords[i].subtableOffset, file);
        encodingRecords[i].subtableOffset += tableOffset;
        //printf("subtable index %zu: platformID: %u, encodingID: %u, offset: %u\n", i, encodingRecords[i].platformID, encodingRecords[i].encodingID, encodingRecords[i].subtableOffset);
    }

    EncodingRecord formatToUse;
    formatToUse.subtableOffset = 0;
    // encoding records are allways ordered first by platformID asc, then by encodingID asc.
    // Format 12 should be used first if not supported format 4.
    // TODO: format 14 + format 12 OR 4, if font has pID = 0 & eID = 5.
    for (size_t i = numOfTables; i >= 0; i--) {
        if(platform == 3){
            if(encodingRecords[i].platformID == 3 && encodingRecords[i].encodingID == 10){
                formatToUse = encodingRecords[i];
                break;
            }
            if(encodingRecords[i].platformID == 3 && encodingRecords[i].encodingID == 1 || encodingRecords[i].platformID == 3 && encodingRecords[i].encodingID == 0){
                formatToUse = encodingRecords[i];
                break;
            }
        } else if(platform == 0){
            if(encodingRecords[i].platformID == 0 && encodingRecords[i].encodingID == 4){
                formatToUse = encodingRecords[i];
                break;
            }
            if(encodingRecords[i].platformID == 0 && encodingRecords[i].encodingID == 3){
                formatToUse = encodingRecords[i];
                break;
            }
        } else {
            printf("Platform not supported!\n");
            return true;
        }
    }
    if(formatToUse.subtableOffset == 0){
        printf("no implemented format found.\n");
        free(encodingRecords);
        return true;
    }
    fseek(file, formatToUse.subtableOffset, SEEK_SET);
    uint16_t format;
    ReadUInt16(&format, file);
    printf("selected format: %u\n", format);
    // skip reserved 2 bytes
    SkipBytes(2, file);
    if(format == 12) {
        uint32_t length;
        ReadUInt32(&length, file);
        // skip language, some Macintosh bs.
        SkipBytes(4, file);
        uint32_t numOfGroups;
        ReadUInt32(&numOfGroups, file);
        printf("number of groups: %u\n", numOfGroups);
        uint32_t startCharCode, endCharCode, startGlyphID;
        for (size_t i = 0; i < numOfGroups; i++) {
            ReadUInt32(&startCharCode, file);
            ReadUInt32(&endCharCode, file);
            ReadUInt32(&startGlyphID, file);
            //printf("start char code: 0x%08X\n", startCharCode);
            //printf("end char code: 0x%08X\n", endCharCode);
            //printf("start glyphID: %u\n", startGlyphID);
        }
        
        

    } else if(format == 4){

    } else {
        printf("format %u not supported\n", format);
        free(encodingRecords);
        return true;
    }

    return false;
}

bool BezierInterpolation(Point p0, Point p1, Point p2, float t, Point* result){
    if(!result) return true;

    float u = 1.0f - t;
    result->x = u * u * p0.x + 2 * u * t * p1.x + t * t * p2.x;
    result->y = u * u * p0.y + 2 * u * t * p1.y + t * t * p2.y;
    return false;
}

bool CalcQuadraticRoots(float a, float b, float c, float* retRootA, float* retRootB) {
    *retRootA = NAN;
    *retRootB = NAN;
    if (a == 0) {
        if (b != 0) *retRootA = -c / b;
    } else {
        float discriminant = b * b - 4 * a * c;
        if(discriminant >= 0) {
            float s = sqrt(discriminant);
            *retRootA = (-b + s) / (2 * a);
            *retRootB = (-b - s) / (2 * a);
        }
    }
    return false;
}

bool RayCastIntersections(Point ray, Point p0, Point p1, Point p2, int* numOfIntersections) {
    float a = p0.y - 2 * p1.y + p2.y;
    float b = 2 * (p1.y - p0.y);
    float c = p0.y;
    float t0 = NAN, t1 = NAN;
    CalcQuadraticRoots(a, b, c - ray.y, &t0, &t1);
    int num = 0;

    float eps = 1e-6f;
    for (int i = 0; i < 2; i++) {
        float t = (i == 0) ? t0 : t1;
        if (isnan(t) || t < 0.0f || t > 1.0f) continue;
        Point bez;
        BezierInterpolation(p0, p1, p2, t, &bez);
        if (bez.x > ray.x + eps && bez.y != p0.y) num++;
    }
    
    *numOfIntersections += num;
    return false;
}

bool isPointInsideGlyph(Point castOrigin, Glyph glyph) {
    int numOfIntersections = 0;
    size_t lastPoint = 0;

    int contourStart = 0;
    for (size_t contourIndex = 0; contourIndex < glyph.numOfContours; contourIndex++) {
        int contourEnd = glyph.endPtsOfContours[contourIndex];
        int firstOnCurvePointIndex;
        for (firstOnCurvePointIndex = contourStart; firstOnCurvePointIndex <= contourEnd; firstOnCurvePointIndex++) {
            if(glyph.points[firstOnCurvePointIndex].isOnCurve){
                break;
            }
        }
        for (size_t i = firstOnCurvePointIndex; i <= contourEnd; i += 2) {
            int startIndex = i + firstOnCurvePointIndex - contourStart;
            int midIndex = (i + firstOnCurvePointIndex + 1 - contourStart) % (contourEnd + 1);
            int endIndex = (i + firstOnCurvePointIndex + 2 - contourStart) % (contourEnd + 1);
            if(midIndex < startIndex) midIndex += contourStart;
            if(endIndex < startIndex) endIndex += contourStart;
            Point start = glyph.points[startIndex];
            Point mid = glyph.points[midIndex];
            Point end = glyph.points[endIndex];

            RayCastIntersections(castOrigin, start, mid, end, &numOfIntersections);
        }
        

        contourStart = contourEnd + 1;
    }
    return (numOfIntersections % 2) != 0;
}

void GlyphToBitmap(uint8_t *bitmap, Glyph glyph, int scale, int16_t xMin, int16_t yMin, int16_t xMax, int16_t yMax) {
    for (size_t i = 0; i < glyph.numOfPoints; i++) {
        printf("point %zu: (%d,%d)\n", i, glyph.points[i].x, glyph.points[i].y);
    }
    
    FILE *fp = fopen("output.txt", "w");

    int width = xMax - xMin;
    int height = yMax - yMin;
    //int16_t xMax = glyph.bounds.xMax / scale;
    //int16_t yMax = glyph.bounds.yMax / scale;
    printf("bounds: (%d, %d) (%d, %d)\n", xMin, yMin, xMax, yMax);
    for (int y = yMin; y < yMax; y++) {
        for (int x = xMin; x < xMax; x++) {
            int idx = (y - yMin) * width + (x - xMin);
            Point rayOrigin;
            rayOrigin.x = x / scale;
            rayOrigin.y = y / scale;
            if(isPointInsideGlyph(rayOrigin, glyph)){
                bitmap[idx] = 255;
            }
            
            for (size_t i = 0; i < glyph.numOfPoints; i++) {
                int16_t distx = glyph.points[i].x / scale;
                int16_t disty = glyph.points[i].y / scale;
                if((distx - x)*(distx - x) + (disty - y)*(disty - y) < (50 * (i + 1))) {
                    
                    if(glyph.points[i].isOnCurve) {
                        bitmap[idx] = 100;
                    } else {
                        bitmap[idx] = 200;
                    }
                
                    break;
                }
            
            }
            
            fprintf(fp, "%u", bitmap[idx]);
            if (x < xMax - 1) fprintf(fp, ",");

        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int main(int argc, char** argv){
    FILE* file = fopen("../../Fonts/DejaVuSans.ttf", "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }
    Font font;
    
    // Offset subtable
    // 12 byte
    // skip the magic number
    SkipBytes(4, file);
    // next 2 bytes are the number of tables
    
    ReadUInt16(&font.numOfTables, file);
    printf("numOfTables: %u\n", font.numOfTables);
    // Currently skipping the rest of the table 2 bytes each
    // TODO: searchRange, entrySelector and rangeShift 
    SkipBytes(6, file);

    if(CreateHashMap(&font.tables, font.numOfTables)){
        printf("failed to create!\n");
        return 1;
    }
    // Table directory 
    // Each table has a directory
    // tag (4), checkSum (4), offset (4), length (4)
    for (size_t i = 0; i < font.numOfTables; i++) {
        char tag[5] = {0};
        uint32_t checkSum;
        uint32_t offset;
        uint32_t length;
        ReadTag(tag, file);
        ReadUInt32(&checkSum, file);
        ReadUInt32(&offset, file);
        ReadUInt32(&length, file);
        // printf("%c%c%c%c offset: %u\n", tag[0], tag[1], tag[2], tag[3], offset);
        if(InsertPair(&font.tables, tag, &offset)){
            printf("failed to insert!\n");
            return 1;
        }
    }
    // Find the number of glyphs
    JumpToTable("maxp", &font.tables, file);
    SkipBytes(4, file);
    ReadUInt16(&font.numOfGlyphs, file);
    printf("number of Glyphs: %u\n", font.numOfGlyphs);

    GetOffsetOfAllGlyphs(&font, file);
    GetCharMapping(&font.tables, file);
    
    // Go to glyph table
    JumpToTable("glyf", &font.tables, file);

    font.glyphs = calloc(font.numOfGlyphs, sizeof(Glyph));
    /*
    for (size_t i = 0; i < 40; i++) {
        ReadAnyGlyph(&font, i, file);
    }
    */
   int glyphToPrint = 62;
    ReadAnyGlyph(&font, glyphToPrint, file);
    uint8_t* bitmap;
    int scale = 1;
    int16_t xMax = -32768, yMax = -32768;
    int16_t xMin = 32767, yMin = 32767;
    for (size_t i = 0; i < font.glyphs[glyphToPrint].numOfPoints; i++){
        if(xMax < font.glyphs[glyphToPrint].points[i].x){
            xMax = font.glyphs[glyphToPrint].points[i].x;
        }
        if(yMax < font.glyphs[glyphToPrint].points[i].y){
            yMax = font.glyphs[glyphToPrint].points[i].y;
        }
        if(xMin > font.glyphs[glyphToPrint].points[i].x){
            xMin = font.glyphs[glyphToPrint].points[i].x;
        }
        if(yMin > font.glyphs[glyphToPrint].points[i].y){
            yMin = font.glyphs[glyphToPrint].points[i].y;
        }
    }
    xMax += 200;
    yMax += 200;
    xMin -= 200;
    yMin -= 200;
    bitmap = calloc((abs(yMin) + yMax)  * (abs(xMin) + xMax) / scale, sizeof(uint8_t));
    printf("glyph nOfCon: %d, nOfPoints: %u\n", font.glyphs[glyphToPrint].numOfContours, font.glyphs[glyphToPrint].numOfPoints);
    GlyphToBitmap(bitmap, font.glyphs[glyphToPrint], scale, xMin, yMin, xMax, yMax);
    free(bitmap);
    DestroyFont(&font);

    fclose(file);
    return 0;
}