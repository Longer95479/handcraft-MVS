// mkdir output
// mkdir build
// cd output
// g++ -O3 -o ../build/shader ../src/shader2.cpp
// ../build/shader
// ffmpeg -i output%02d.ppm -r 60 output.mp4

#include <iostream>
#include <cmath>


double abs(double x) { return fabs(x); }


struct vec4;
struct vec3;
struct vec2;


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

  vec2 xz();
  vec3 zxy();
  vec3 xyy();
  vec3 yxy();
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

  vec2 yx();
  vec3 xyy();
  vec4 xyyx();

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


vec2 vec2::yx() { return vec2(y, x); }
vec3 vec2::xyy() { return vec3(x, y, y); }
vec4 vec2::xyyx() { return vec4(x, y, y, x); }

vec2 vec3::xz() { return vec2(x, z); }
vec3 vec3::zxy() { return vec3(z, x, y); }
vec3 vec3::xyy() { return vec3(x, y, y); }
vec3 vec3::yxy() { return vec3(y, x, y); }


int main(int argc, char**argv)
{
  char buf[256];
  for (int j = 0; j < 60; j++) {
    snprintf(buf, sizeof(buf), "output%02d.ppm", j);
    const char* output_path = buf;
    FILE* f = fopen(output_path, "wb");
  
    int w = 16 * 60;
    int h = 9 * 60;
  
    fprintf(f, "P6\n");
    fprintf(f, "%d %d\n", w, h);
    fprintf(f, "255\n");

    double t = (double) j / 60;
    vec2 r((double)w, (double)h);

    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        vec2 FC((double)x, (double)y);
        vec3 rgb((double)x, (double)y, 122.0);
        vec4 o;

        // plasma
        // vec2 p=(FC*2.-r)/r.y,l,i,v=p*(l+=4.-4.*abs(.7-dot(p,p)));for(;i.y++<8.;o+=(sin(v.xyyx())+1.)*abs(v.x-v.y))v+=cos(v.yx()*i.y+i+t)/i.y+.7;o=tanh(5.*exp(l.x-4.-p.y*vec4(-1,1,2,0))/o);
        
        // cyberspace 1
        // for(double i=0,z=0,d=0;z+i++<7e1;o+=vec4(z,1,9,1)/d)
        // {
        //   vec3 p=abs(z*normalize(rgb*2.-r.xyy()));
        //   p.z+=t*5.;p+=sin(p+p);
        //   for(d=0.;d++<9.;p+=.4*cos(round(.2*d*p)+.2*t).zxy());
        //   z+=d=.1*sqrt(length(p.xyy()*p.yxy()));
        // }
        // o=tanh(o/7e3);

        // cyberspace 2
        for(float i=0,z=0,d=0;z+i++<8e1;o+=vec4(z,4,1,1)/d)
        {vec3 p=z*normalize(rgb*2.-r.xyy());p.z+=t/.1;p-=sin(p+p);for(d=0.;d++<9.;p+=.7*cos(round(.2*d*p)+t*.5).zxy());z+=d=.1*sqrt(length(p.xyy()*p.yxy()));}
        o=tanh(o/5e3);

        fputc(o.x * 255, f);
        fputc(o.y * 255, f);
        fputc(o.z * 255, f);
      }
    }
  
    fclose(f);
    printf("generated ppm file: %s\n", output_path);

  }
  return 0;
}

