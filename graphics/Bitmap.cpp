#include <vector>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include "megeometry.h"

using namespace std;
struct Pixel {
    unsigned char r, g, b;
    Pixel(): r(0), g(0), b(0){};
    Pixel(const Vec3f &pix): r(pix.z), g(pix.y), b(pix.x){}
};

void WriteBMP(const char* fname, Pixel* a_pixelData, int width, int height)
{
  int paddedsize = (width*height) * sizeof(Pixel);

  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};

  bmpfileheader[ 2] = (unsigned char)(paddedsize    );
  bmpfileheader[ 3] = (unsigned char)(paddedsize>> 8);
  bmpfileheader[ 4] = (unsigned char)(paddedsize>>16);
  bmpfileheader[ 5] = (unsigned char)(paddedsize>>24);

  bmpinfoheader[ 4] = (unsigned char)(width    );
  bmpinfoheader[ 5] = (unsigned char)(width>> 8);
  bmpinfoheader[ 6] = (unsigned char)(width>>16);
  bmpinfoheader[ 7] = (unsigned char)(width>>24);
  bmpinfoheader[ 8] = (unsigned char)(height    );
  bmpinfoheader[ 9] = (unsigned char)(height>> 8);
  bmpinfoheader[10] = (unsigned char)(height>>16);
  bmpinfoheader[11] = (unsigned char)(height>>24);

  std::ofstream out(fname, std::ios::out | std::ios::binary);
  out.write((const char*)bmpfileheader, 14);
  out.write((const char*)bmpinfoheader, 40);
  out.write((const char*)a_pixelData, paddedsize);
  out.flush();
  out.close();
}

void SaveBMP(const char* fname, const vector<Pixel> &pixels, int w, int h)
{
  std::vector<Pixel> pixels2(w*h);

  for (size_t i = 0; i < pixels2.size(); i++)
  {
    Pixel px;
    px.r       = pixels[i].r;
    px.g       = pixels[i].g;
    px.b       = pixels[i].b;
    pixels2[i] = px;
  }
  WriteBMP(fname, &pixels2[0], w, h);
}
//void ReadBMP(const char* fname, vector<Vec3f> &pixels, int& width1, int& height1){
 vector <Vec3f> ReadBMP(const char* fname, int& width1, int& height1){
    int size = 0, pixels_adress = 0, width = 0, height = 0;
    short int bits_per_pixel = 0;
    fstream file;
    file.open(fname, ios::in | ios::binary);
    file.seekg( 2, ios::beg);
    file.read ((char*)&size, sizeof(int));
    file.seekg( 10, ios::beg);
    file.read ((char*)&pixels_adress, sizeof(int));
    file.seekg( 18, ios::beg);
    file.read ((char*)&width, sizeof(int));
    file.read ((char*)&height, sizeof(int));
    file.seekg( 28, ios::beg);
    file.read ((char*)&bits_per_pixel, sizeof(short int));
    file.seekg( pixels_adress, ios::beg);

    unsigned int bgr;
    vector<Vec3f> pixels (width *  height);

    for (int j = 0; j<height; j++) {
        for (int i = 0; i < width; i++) {
            file.read((char *) &bgr, 4);
            pixels[i + (height - 1 - j) * width].x = ((bgr & 0x00FF0000) >> 16) / 255.;
            pixels[i + (height - 1 - j) * width].y = ((bgr & 0x0000FF00) >> 8)/ 255.;
            pixels[i + (height - 1 - j) * width].z = (bgr & 0x000000FF)/ 255.;

        }
    }
    width1 = width;
    height1 = height;
    return  pixels;
}
