
#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <utility>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <cmath>
#include <omp.h>
#include "Bitmap.h"
#include "megeometry.h"


using namespace std;


static const float kInfinity = numeric_limits<float>::max();
static const float kEpsilon = 1e-8;

vector<Vec3f> envmap;
int envmap_width;
int envmap_height;

inline
float clamp(const float &lo, const float &hi, const float &v)
{ return max(lo, min(hi, v)); }

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

struct Material {
    Material(const float r, const Vec4f &a, const Vec3f &color, const float &spec) :
            refractive_index(r), //преломление
            albedo(a),
            diffuse_color(color),
            specular_exponent(spec) {} // степень неровности поверхности

    Material() : refractive_index(1.0), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

struct Light {
    Light(const Vec3f &p, const float /*&*/i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Object {
    Object() {};
    Material material;
    virtual ~Object() {}
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;
    virtual void getSurfaceProperties(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
};


struct Options
{
    uint32_t width = 512;
    uint32_t height = 512;
    float fov = 90;
    vector <Light> lights;
    Vec3f backgroundColor;
    Matrix44f cameraToWorld;
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;
    Sphere(const Vec3f &c, const float r, const Material &m) {
        center = c;
        radius = r;
        material = m;
    }
    //Нахождение пересечения луча и сферы
    // p = origin + dir*t - уравнение прямой
    // || p - c || = r - уравнение сферы
    bool ray_intersect(const Vec3f &origin, const Vec3f &dir, float &t) const {
        const Vec3f L = center - origin;  // вектор из начала луча к центру сферы
        const float tca = L * dir; // кратчайшее растояние от луча до центра

        float d2 = L * L - tca * tca;
        float t2;
        float t2hc = radius * radius - d2;

        if (t2hc < 0) return false; // пересечения нет
        t2hc = sqrtf(t2hc);
        t = tca - t2hc;
        t2 = tca + t2hc;
        if (t < 0) t = t2;
        if (t < 0) return false;
        return true;
    }
};
struct Triangle{
    Vec3f A, B, C;
    Material material;
    Triangle(const Vec3f &a, const Vec3f &b, const Vec3f &c,const Material &m){
        A = a;
        B = b;
        C = c;
        material = m;
    }
    bool ray_intersect(const Vec3f &origin, const Vec3f &dir, float &t) const{
        const Vec3f e1 = B - A;
        const Vec3f e2 = C - A;
        const Vec3f s = origin - A;
        const Vec3f x = cross(dir, e2);  // [dir,e2]
        const Vec3f y = cross(s, e1);
        const float d = e1*x; // det(A,B,C)

        if (fabs(d) < 1e-3) return false ; // Система не имеет решения
        const float u = (s*x) * (1.0/d);
        if (u < 0 || u > 1) return false;
        const float v = (dir*y) * (1.0/d);
        if (v < 0 || v > 1 || u + v > 1) return false;
        t = (e2 * y) * (1.0/d);
        return true;
    }
};

//Нахождение пересечения луча и треугольника
// p0 + lt = (1-u-v)a + ub + vc, u,v из [0,1], u+v из [0,1]
bool rayTriangleIntersect(const Vec3f &origin, const Vec3f &dir,const Vec3f &A,const Vec3f &B,const Vec3f &C, float &t, float &u, float &v){

    const Vec3f e1 = B - A;
    const Vec3f e2 = C - A;
    const Vec3f s = origin - A;
    const Vec3f x = cross(dir, e2);  // [dir,e2]
    const Vec3f y = cross(s, e1);
    const float d = e1*x; // det(A,B,C)

    if (fabs(d) < kEpsilon) return false ; // Система не имеет решения
    u = (s*x) * (1.0/d);
    if (u < 0 || u > 1) return false;
    v = (dir*y) * (1.0/d);
    if (v < 0 || v > 1 || u + v > 1) return false;
    t = (e2 * y) * (1.0/d);
    return true;
}
struct TriangleMesh: public Object{
    TriangleMesh(
            const uint32_t nfaces, // размер массива
            const unique_ptr<uint32_t []> &faceIndex,
            const unique_ptr<uint32_t []> &vertsIndex, // массив индексов вершин
            const unique_ptr<Vec3f []> &verts, // массив вершин
            unique_ptr<Vec3f []> &normals, // координаты нормали
            unique_ptr<Vec2f []> &st): // координаты текстур
            numTris(0) // колличесвто вершин
    {
        uint32_t k= 0; // зацикливание граней
        uint32_t maxVertIndex = 0; // максимальный размер массива вершин
        //колличество вершин
        for (uint32_t i = 0; i < nfaces; ++i) {
            numTris += faceIndex[i] - 2; // число треугольников всегда равно числу вершин в многоугольнике минус 2
            for (uint32_t j = 0; j < faceIndex[i]; ++j)
                if (vertsIndex[k + j] > maxVertIndex)
                    maxVertIndex = vertsIndex[k + j];
            k += faceIndex[i];
        }
        maxVertIndex += 1;
        //Положение вершин сетки
        P = unique_ptr<Vec3f []>(new Vec3f[maxVertIndex]);
        for (uint32_t i = 0; i < maxVertIndex; ++i) {
            P[i] = verts[i];
        }

        // хранение индексов треугольников
        trisIndex = unique_ptr<uint32_t []>(new uint32_t [numTris * 3]); // для каждого треугольника 3 индекса
        uint32_t l = 0;
        N = unique_ptr<Vec3f []>(new Vec3f[numTris * 3]);
        texCoordinates = unique_ptr<Vec2f []>(new Vec2f[numTris * 3]);
        for (uint32_t i = 0, k = 0; i < nfaces; ++i) {
            for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) { // для каждого треугольника на грани
                trisIndex[l] = vertsIndex[k];
                trisIndex[l + 1] = vertsIndex[k + j + 1];
                trisIndex[l + 2] = vertsIndex[k + j + 2];
                N[l] = normals[k];
                N[l + 1] = normals[k + j + 1];
                N[l + 2] = normals[k + j + 2];
                texCoordinates[l] = st[k];
                texCoordinates[l + 1] = st[k + j + 1];
                texCoordinates[l + 2] = st[k + j + 2];
                l += 3;
            }
            k += faceIndex[i];
        }
    }
    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &triIndex, Vec2f &uv) const
    {
        uint32_t j = 0;
        bool isect = false;
        for (uint32_t i = 0; i < numTris; ++i) {
            const Vec3f &v0 = P[trisIndex[j]];
            const Vec3f &v1 = P[trisIndex[j + 1]];
            const Vec3f &v2 = P[trisIndex[j + 2]];
            float t = kInfinity, u, v;
            if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v) && t < tNear) {
                tNear = t;
                uv.x = u;
                uv.y = v;
                triIndex = i;
                isect = true;
            }
            j += 3;
        }

        return isect;
    }
    void getSurfaceProperties(
            const Vec3f &hitPoint,
            const Vec3f &viewDirection,
            const uint32_t &triIndex,
            const Vec2f &uv,
            Vec3f &hitNormal,
            Vec2f &hitTextureCoordinates) const
    {
        //Нормаль
        const Vec3f &v0 = P[trisIndex[triIndex * 3]];
        const Vec3f &v1 = P[trisIndex[triIndex * 3 + 1]];
        const Vec3f &v2 = P[trisIndex[triIndex * 3 + 2]];
        hitNormal = cross(v1 - v0,v2 - v0);
        hitNormal.normalize();

        // координаты текстур
        const Vec2f &st0 = texCoordinates[triIndex * 3];
        const Vec2f &st1 = texCoordinates[triIndex * 3 + 1];
        const Vec2f &st2 = texCoordinates[triIndex * 3 + 2];
        hitTextureCoordinates = st0 * (1 - uv.x - uv.y)  + st1 * uv.x   +  st2 * uv.y;

    }

    uint32_t numTris;
    unique_ptr<Vec3f []> P;
    unique_ptr<uint32_t []> trisIndex;
    unique_ptr<Vec3f []> N;
    unique_ptr<Vec2f []> texCoordinates;
};
// Отражение вектора относительно нормали N
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f){
    float cosi = - max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t);
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k));
}

//Считывание данных из файла
TriangleMesh* loadPolyMeshFromFile(const char *file)
{
    ifstream ifs;
    try {
        ifs.open(file);
        if (ifs.fail()) throw;
        stringstream ss;
        ss << ifs.rdbuf();
        uint32_t numFaces;
        ss >> numFaces; //1 число в файле - колличество граней в сетке
        unique_ptr<uint32_t []> faceIndex(new uint32_t[numFaces]);
        uint32_t vertsIndexArraySize = 0;

        for (uint32_t i = 0; i < numFaces; ++i) {
            ss >> faceIndex[i];
            vertsIndexArraySize += faceIndex[i];
        }


        unique_ptr<uint32_t []> vertsIndex(new uint32_t[vertsIndexArraySize]);
        uint32_t vertsArraySize = 0;

        // Массив индекса вершин
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> vertsIndex[i];
            if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
        }
        vertsArraySize += 1;

        // чтение вершин
        std::unique_ptr<Vec3f []> verts(new Vec3f[vertsArraySize]);
        for (uint32_t i = 0; i < vertsArraySize; ++i) {
            ss >> verts[i].x >> verts[i].y >> verts[i].z;
        }

        // чтение нормалей
        std::unique_ptr<Vec3f []> normals(new Vec3f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> normals[i].x >> normals[i].y >> normals[i].z;
        }

        // чтение координат текстур
        std::unique_ptr<Vec2f []> st(new Vec2f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> st[i].x >> st[i].y;
        }

        return new TriangleMesh(numFaces, faceIndex, vertsIndex, verts, normals, st);
    }
    catch (...) {
        ifs.close();
    }
    ifs.close();

    return nullptr;
}

bool trace(
        const Vec3f &orig, const Vec3f &dir,
        const vector<unique_ptr<Object>> &objects,
        float &tNear, uint32_t &index, Vec2f &uv, Object **hitObject){
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearTriangle = kInfinity;
        uint32_t indexTriangle;
        Vec2f uvTriangle;
        if (objects[k]->intersect(orig, dir, tNearTriangle, indexTriangle, uvTriangle) && tNearTriangle < tNear) {
            *hitObject = objects[k].get();
            tNear = tNearTriangle;
            index = indexTriangle;
            uv = uvTriangle;
        }
    }

    return (*hitObject != nullptr);
}

Vec3f castRay(
        const Vec3f &orig, const Vec3f &dir,
        const vector<std::unique_ptr<Object>> &objects,
         Vec3f hitColor, const vector<Light> lights)
{
    float tnear = kInfinity;
    Vec2f uv;
    uint32_t index = 0;
    Object *hitObject = nullptr;

    Material material(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    if (trace(orig, dir, objects, tnear, index, uv, &hitObject)){
        Vec3f hitPoint = orig + dir * tnear;
        Vec3f hitNormal;
        Vec2f hitTexCoordinates;
        hitObject->getSurfaceProperties(hitPoint, dir, index, uv, hitNormal, hitTexCoordinates);
        float NdotView = max(0.f, hitNormal*(-dir));
        const int M = 5;
        float checker = (fmod(hitTexCoordinates.x * M, 1.0) > 0.5) ^ (fmod(hitTexCoordinates.y * M, 1.0) < 0.5);
        float c = 0.2 * (1 - checker) + 0.8 * checker;

        hitColor = c * NdotView;
        material.diffuse_color = hitColor;
        float diffuse_light_intensity = 0, specular_light_intensity = 0;
        for (size_t i=0; i<lights.size(); i++) {
           Vec3f light_dir = (lights[i].position - hitPoint).normalize();
           diffuse_light_intensity += lights[i].intensity * max(0.f, light_dir * hitNormal);
           specular_light_intensity += powf(max(0.f, -reflect(-light_dir, hitNormal) * dir), material.specular_exponent) * lights[i].intensity;
       }
        Vec3f p = material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1];
        float max = std::max(p[0], std::max(p[1], p[2]));
        if (max > 1) p = p * (1. / max);
        return p;
    }
    return hitColor;
}

vector<Pixel> render(
        const Options &options,
        const std::vector<std::unique_ptr<Object>> &objects)
{

    float scale = tan(deg2rad(options.fov * 0.5));
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);

    //auto timeStart = std::chrono::high_resolution_clock::now();

    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    vector<Pixel> image(options.width * options.height);
#pragma omp parallel for
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();
            Vec3f px = castRay(orig, dir, objects, envmap[i+j*options.width], lights);
            image[i + (options.height - 1 - j) * options.width].b = (char) (255 * clamp(0, 1,px.x));
            image[i + (options.height - 1 - j) * options.width].g = (char) (255 * clamp(0, 1, px.y));
            image[i + (options.height - 1 - j) * options.width].r = (char) (255 * clamp(0, 1, px.z));
        }
        //fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
    }
   // auto timeEnd = chrono::high_resolution_clock::now();
   // auto passedTime = chrono::duration<double, std::milli>(timeEnd - timeStart).count();
    //fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

    return image;
}


bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material, const vector<Triangle> &triangles) {
    float spheres_dist = numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir * dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    float triangles_dist = numeric_limits<float>::max();
    for (size_t i=0; i < triangles.size(); i++) {
        float dist_i;
        float u,v;
        if(rayTriangleIntersect(orig, dir,triangles[i].A, triangles[i].B, triangles[i].C, dist_i,u,v) && dist_i < triangles_dist && dist_i <spheres_dist){
            material = triangles[i].material;
            material.diffuse_color.x = u* triangles[i].material.diffuse_color.x;
            material.diffuse_color.y = v* triangles[i].material.diffuse_color.y;
            material.diffuse_color.z = (1-u-v)* triangles[i].material.diffuse_color.z;
            N = cross(triangles[i].B - triangles[i].A, triangles[i].C - triangles[i].A).normalize();
            float d = triangles[i].A * N;
            //float d = triangles[i].C * N;
            float t = (orig * N + d) / (dir * N);
            Vec3f hit = orig + dir * t;
            triangles_dist = dist_i;
        }

    }
    float checkerboard_dist = numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
         float d = -(orig.y+4)/dir.y;
         Vec3f pt = orig + dir*d;
         if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheres_dist && d <triangles_dist) {
             checkerboard_dist = d;
             hit = pt;
             N = Vec3f(0,1,0);
             material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(221/255., 174/255., 206/255.) : Vec3f(209/255., 174/255., 206/255.);
         }
     }
    // Бесконечная плоскость
    /*float plane_dist = numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
        float dist_i = - (orig.y + 4) / dir.y;
        if (dist_i>0 && dist_i<spheres_dist) {
            plane_dist = dist_i;
            hit = orig + dir * dist_i;
            N = Vec3f(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.0, .0, .0) : Vec3f(.1, .1, .1);

            // material.diffuse_color = Vec3f(0,0,0.4);
        }
    }*/
  //  return min(triangles_dist, min(spheres_dist, plane_dist))<1000;
    return min(triangles_dist, min(spheres_dist, checkerboard_dist)) <1000;

}
Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const vector<Sphere> &spheres, const vector<Light> &lights, const vector<Triangle> &triangles, size_t depth=0) {
    Vec3f point, N;
    Material material;

    if (depth>4 || !scene_intersect(orig, dir, spheres, point, N, material, triangles)) {
        Sphere env(Vec3f(0,0,0), 100, Material());
        float dist = 0;
        env.ray_intersect(orig, dir, dist);
        Vec3f p = orig+dir*dist;
        int a = (atan2(p.z, p.x)/(2*M_PI) + .5)*envmap_width;
        int b = acos(p.y/100)/M_PI*envmap_height;
        return envmap[a+b*envmap_width];
    }
    // Отражение.
    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights,triangles,depth + 1);

    //Преломление
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights,triangles, depth + 1);
    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();

        //Тени
        float light_distance = (lights[i].position - point).norm();
        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // сдвигаем lights[i] точку в направлении нормали
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        // Есои источник пересекает объект, пропускаем его.
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial,triangles) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        //Модель Ламберта - диффузная модель
        diffuse_light_intensity  += lights[i].intensity * max(0.f, light_dir*N);
        // Модель Фонга
        specular_light_intensity += powf(max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;

    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

int main(int argc, char **argv)
{
    unordered_map<string, string> cmdLineParams;

    for(int i=0; i<argc; i++)
    {
        string key(argv[i]);

        if(key.size() > 0 && key[0]=='-')
        {
            if(i != argc-1) // not last argument
            {
                cmdLineParams[key] = argv[i+1];
                i++;
            }
            else
                cmdLineParams[key] = "";
        }
    }

    string outFilePath = "zout.bmp";
    if(cmdLineParams.find("-out") != cmdLineParams.end())
        outFilePath = cmdLineParams["-out"];

    int sceneId = 0;
    if(cmdLineParams.find("-scene") != cmdLineParams.end())
        sceneId = atoi(cmdLineParams["-scene"].c_str());

    if(sceneId == 2) {
        Options options;
        Matrix44f tmp = Matrix44f(0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295,
                                  0.624695, 0, -1.63871, -5.747777, -40.400412, 1);
        options.cameraToWorld = tmp.inverse();
        options.fov = 50.0393;

        vector <std::unique_ptr<Object>> objects;
        TriangleMesh *mesh = loadPolyMeshFromFile("../cow.geo");
        if (mesh != nullptr)
            objects.push_back(unique_ptr<Object>(mesh));
        envmap = ReadBMP("../background.bmp", envmap_width, envmap_height);
        SaveBMP(outFilePath.c_str(), render(options, objects), options.width, options.height);
        cout << "end.";
    }
    else if(sceneId == 1){
        Material      ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
        Material      glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
        Material red_rubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);

        Material sienna(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(139/255., 69/255., 55/255.),   50.);

        Material     mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);
        Material      grey(1.0, Vec4f(0.9,  0.3, 0.1, 0.0), Vec3f(0.2, 0.2, 0.2),  10.);
        Material      rubber(1.5, Vec4f(0.9,  0.5, 0.0, 0.9), Vec3f(0.9, 0.1, 0.1),   100.);
        Material      orange(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1, 0.4, 0.2),   50.);



        vector <Triangle> triangles;
        triangles.push_back(Triangle(Vec3f(0.4, -0.2, -1),Vec3f(0.5,  -0.2, -0.7),Vec3f(0.5, 0, -1), orange));


        std::vector<Sphere> spheres;
        spheres.push_back(Sphere(Vec3f(-3,    -1,   -16), 2.6,      sienna));
        spheres.push_back(Sphere(Vec3f(-4.3,    -3,   -14.5), 1,      sienna));
        spheres.push_back(Sphere(Vec3f(-1.7,    -3,   -14.5), 1,      sienna));
        spheres.push_back(Sphere(Vec3f(-0.7,    -0.5,   -14.5), 0.9,      sienna));
        spheres.push_back(Sphere(Vec3f(-5.2,    -0.5,   -14.5), 0.9,      sienna));
        spheres.push_back(Sphere(Vec3f(-3,    2.2,   -15), 1.7,      sienna));
        //Уши
        spheres.push_back(Sphere(Vec3f(-5,    3.8,   -17), 0.9,      sienna));
        spheres.push_back(Sphere(Vec3f(-1.7,    3.8,   -16.8), 0.9,      sienna));

        //Глаза
        spheres.push_back(Sphere(Vec3f(-3.4,    2.7,   -13.8), 0.5,      grey));
        spheres.push_back(Sphere(Vec3f(-2.6,    2.7,   -13.8), 0.5,      grey));

        spheres.push_back(Sphere(Vec3f( 0,    5,   -25), 7,     mirror));
        spheres.push_back(Sphere(Vec3f(3, -2.3, -12), 1,      glass));
        spheres.push_back(Sphere(Vec3f( 3.6, -1, -16), 2,      red_rubber));
        spheres.push_back(Sphere(Vec3f(1,    -2.7,   -15), 1.4,      orange));

        std::vector<Light>  lights;
        lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
        lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
        lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

        const int width    = 1024;
        const int height   = 724;
        const int fov      = M_PI/3.;

        vector<Pixel> image(width*height);
        envmap = ReadBMP("../environment.bmp", envmap_width, envmap_height);

#pragma omp parallel for
        for (uint32_t j = 0; j < height; ++j)
            for (uint32_t i = 0; i < width; ++i){
                float a = (2 * (i + 0.5) / (float) width - 1) * tan(fov / 2.) * width / (float) height;
                float b = -(2 * (j + 0.5) / (float) height - 1) * tan(fov / 2.);
                Vec3f dir = Vec3f(a, b, -1).normalize();
                Vec3f temp = cast_ray(Vec3f(0, 0, 0), dir, spheres, lights, triangles);

                Vec3f &c = temp;
                float max = std::max(c[0], std::max(c[1], c[2]));
                if (max > 1) c = c * (1. / max);

                for (size_t j = 0; j < 3; j++) {
                    temp[j] = 255 * std::max(0.f, min(1.f, temp[j]));
                }

                Pixel px(temp);
                image[i + (height - 1 - j) * width] = px;
        }
        SaveBMP(outFilePath.c_str(), image, width, height);
        cout << "end.";
    }
    else
        cout << "\nSorry, there is no such scene! \n";
    return 0;
}
