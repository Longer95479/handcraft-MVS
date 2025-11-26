#ifndef FRAG_SHADER_H
#define FRAG_SHADER_H

#include <iostream>
#include <cmath>

void fragShader(double T, int frame_size, int w, int h, const auto& lambda);

struct vec4;
struct vec3;
struct vec2;


//-----------------------------------------------------------------------------//

double abs(double x) { return fabs(x); }

struct vec4 {
  double x, y, z, w;

  vec4(): x(.0), y(.0), z(.0), w(.0) {}
  vec4(double _x, double _y, double _z, double _w): x(_x), y(_y), z(_z), w(_w) {}

};

vec4 operator + (const vec4& a, const vec4& b) { return vec4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w); }
vec4 operator + (const vec4& a, double s) { return vec4(a.x + s, a.y + s, a.z + s, a.w + s); }
vec4 operator - (const vec4& a, const vec4& b) { return vec4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w); }
vec4 operator - (double s, const vec4& a) { return vec4(s - a.x, s - a.y, s - a.z, s - a.w); }
vec4 operator * (const vec4& a, double s) { return vec4(a.x * s, a.y * s, a.z * s, a.w * s); }
vec4 operator * (double s, const vec4& a) { return vec4(a.x * s, a.y * s, a.z * s, a.w * s); }
vec4 operator / (const vec4& a, const vec4& b) { return vec4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w); }
vec4 operator / (const vec4& a, double s) { return vec4(a.x / s, a.y / s, a.z / s, a.w / s); }
vec4& operator += (vec4& a, const vec4& b) { a = a + b; return a; }
vec4 sin(const vec4& a) { return vec4(sin(a.x), sin(a.y), sin(a.z), sin(a.w)); }
vec4 cos(const vec4& a) { return vec4(cos(a.x), cos(a.y), cos(a.z), cos(a.w)); }
vec4 exp(const vec4& a) { return vec4(exp(a.x), exp(a.y), exp(a.z), exp(a.w)); }
vec4 tanh(const vec4& a) { return vec4(tanh(a.x), tanh(a.y), tanh(a.z), tanh(a.w)); }


struct vec3 {
  double x, y, z;

  vec3(): x(.0), y(.0), z(.0) {}
  vec3(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}

  vec2 xz() const;
  vec3 zxy() const;
  vec3 xyy() const;
  vec3 yxy() const;
};

vec3 operator + (const vec3& a, const vec3& b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
vec3 operator + (const vec3& a, double s) { return vec3(a.x + s, a.y + s, a.z + s); }
vec3 operator - (const vec3& a, const vec3& b) { return vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
vec3 operator * (const vec3& a, const vec3& b) { return vec3(a.x * b.x, a.y * b.y, a.z * b.z); }
vec3 operator * (const vec3& a, double s) { return vec3(a.x * s, a.y * s, a.z * s); }
vec3 operator * (double s, const vec3& a) { return vec3(a.x * s, a.y * s, a.z * s); }
vec3 operator / (const vec3& a, const vec3& b) { return vec3(a.x / b.x, a.y / b.y, a.z / b.z); }
vec3 operator / (const vec3& a, double s) { return vec3(a.x / s, a.y / s, a.z / s); }
vec3& operator += (vec3& a, double s) { a = a + s; return a; }
vec3& operator += (vec3& a, const vec3& b) { a = a + b; return a; }
vec3& operator -= (vec3& a, const vec3& b) { a = a - b; return a; }
double dot(const vec3& a, const vec3& b) { return a.x * b.x +  a.y * b.y + a.z * b.z; }
vec3 cross(const vec3& a, const vec3& b) { return vec3(-a.z * b.y + a.y * b.z, a.z * b.x - a.x * b.z, -a.y * b.x + a.x * b.y); }
vec3 sin(const vec3& a) { return vec3(sin(a.x), sin(a.y), sin(a.z)); }
vec3 cos(const vec3& a) { return vec3(cos(a.x), cos(a.y), cos(a.z)); }
vec3 exp(const vec3& a) { return vec3(exp(a.x), exp(a.y), exp(a.z)); }
vec3 tanh(const vec3& a) { return vec3(tanh(a.x), tanh(a.y), tanh(a.z)); }
vec3 mix(const vec3& a, const vec3& b, double t) { return a * (1 - t) + b * t; }
vec3 normalize(const vec3& a) { vec3 tmp = a * a; return a / sqrt(tmp.x + tmp.y + tmp.z); }
vec3 round(const vec3& a) { return vec3(round(a.x), round(a.y), round(a.z)); }
vec3 abs(const vec3& a) { return vec3(abs(a.x), abs(a.y), abs(a.z)); }
double length(const vec3& a) { return sqrt(dot(a, a)); }


struct vec2 {
  double x, y;

  vec2(): x(.0), y(.0) {}
  vec2(double _x, double _y): x(_x), y(_y) {}

  vec2 yx() const;
  vec3 xyy() const;
  vec4 xyyx() const;

};

vec2 operator + (const vec2& a, const vec2& b) { return vec2(a.x + b.x, a.y + b.y); }
vec2 operator + (const vec2& a, double s) { return vec2(a.x + s, a.y + s); }
vec2 operator - (const vec2& a, const vec2& b) { return vec2(a.x - b.x, a.y - b.y); }
vec2 operator - (double s, const vec2& a) { return vec2(s - a.x, s - a.y); }
vec2 operator * (const vec2& a, const vec2& b) { return vec2(a.x * b.x, a.y * b.y); }
vec2 operator * (const vec2& a, double s) { return vec2(a.x * s, a.y * s); }
vec2 operator * (const double s, const vec2& a) { return vec2(a.x * s, a.y * s); }
vec2 operator / (const vec2& a, double s) { return vec2(a.x / s, a.y / s); }
vec2& operator += (vec2& a, const vec2& b) { a = a + b; return a; }
vec2& operator += (vec2& a, double s) { a = a + s; return a; }
double dot(const vec2& a, const vec2& b) { return a.x * b.x +  a.y * b.y; }
vec2 sin(const vec2& a) { return vec2(sin(a.x), sin(a.y)); }
vec2 cos(const vec2& a) { return vec2(cos(a.x), cos(a.y)); }
vec2 exp(const vec2& a) { return vec2(exp(a.x), exp(a.y)); }
vec2 tanh(const vec2& a) { return vec2(tanh(a.x), tanh(a.y)); }
double length(const vec2& a) { return sqrt(dot(a, a)); }


vec2 vec2::yx() const { return vec2(y, x); }
vec3 vec2::xyy() const { return vec3(x, y, y); }
vec4 vec2::xyyx() const { return vec4(x, y, y, x); }

vec2 vec3::xz() const { return vec2(x, z); }
vec3 vec3::zxy() const { return vec3(z, x, y); }
vec3 vec3::xyy() const { return vec3(x, y, y); }
vec3 vec3::yxy() const { return vec3(y, x, y); }

void fragShader(double T, int frame_size, int w, int h, const auto& lambda)
{
  char buf[256];
  for (int j = 0; j < frame_size; j++) {
    snprintf(buf, sizeof(buf), "output%02d.ppm", j);
    const char* output_path = buf;
    FILE* f = fopen(output_path, "wb");
  
    fprintf(f, "P6\n");
    fprintf(f, "%d %d\n", w, h);
    fprintf(f, "255\n");

    double t = (double) j * T / frame_size;
    vec2 r((double)w, (double)h);

    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        vec2 FC((double)x, (double)y);
        vec3 rgb((double)x, (double)y, j);
        vec4 o;

        lambda(t, j, r, FC, rgb, o);
        
        fputc(o.x * 255, f);
        fputc(o.y * 255, f);
        fputc(o.z * 255, f);
      }
    }
  
    fclose(f);
    printf("generated ppm file: %s\n", output_path);
  }
}

#endif

