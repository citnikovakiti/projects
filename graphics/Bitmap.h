#ifndef BITMAP_GUARDIAN_H
#define BITMAP_GUARDIAN_H
#include <cmath>
#include <vector>
#include "megeometry.h"
struct Pixel {
    unsigned char r, g, b;
    Pixel(): r(0), g(0), b(0){};
    Pixel(const Vec3f &pix): r(pix.z), g(pix.y), b(pix.x){}
};
void SaveBMP(const char* fname, const std::vector<Pixel> &pixels, int w, int h);
std :: vector <Vec3f> ReadBMP(const char* fname, int& width1, int& height1);
#endif
