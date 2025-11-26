// mkdir output
// mkdir build
// cd output
// g++ -O3 -fconcepts -o ../build/shader ../src/shaders/shader_plasma.cpp && ../build/shader
// ffmpeg -i output%02d.ppm -r 60 output.mp4
// mpv output.mp4

#include "frag_shader.h"

int main(int argc, char**argv)
{
    auto sd = [&](double t, int j, 
                  const vec2& r,
                  const vec2& FC,
                  const vec3& rgb,
                  vec4& o) {

        // plasma
        vec2 p=(FC*2.-r)/r.y,l,i,v=p*(l+=4.-4.*abs(.7-dot(p,p)));for(;i.y++<8.;o+=(sin(v.xyyx())+1.)*abs(v.x-v.y))v+=cos(v.yx()*i.y+i+t)/i.y+.7;o=tanh(5.*exp(l.x-4.-p.y*vec4(-1,1,2,0))/o);
        
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
        // for(float i=0,z=0,d=0;z+i++<8e1;o+=vec4(z,4,1,1)/d)
        // {vec3 p=z*normalize(rgb*2.-r.xyy());p.z+=t/.1;p-=sin(p+p);for(d=0.;d++<9.;p+=.7*cos(round(.2*d*p)+t*.5).zxy());z+=d=.1*sqrt(length(p.xyy()*p.yxy()));}
        // o=tanh(o/5e3);
    };

    fragShader(1.0, 60, 16*20, 9*20, sd);
    
    return 0;
}

