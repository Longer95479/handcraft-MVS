// g++ -O3 -fconcepts -o main shader_geo_frag.cpp && ./main

#include <vector>
#include <memory>

#include "frag_shader.h"
#include "../matrixlib/my_matrix_lib.h"

using namespace Loong;

#define INF 1e8

class Object {
  public:
    Object() {}
    Object(double transparency, double reflection,
           const Matrix<double, 3, 1>& surface_color,
           const Matrix<double, 3, 1>& emission_color,
           const Matrix<double, 3, 1>& center)
        : transparency_(transparency), reflection_(reflection),
          surface_color_(surface_color), emission_color_(emission_color),
          center_(center) {}

    virtual bool intersect(const Matrix<double, 3, 1>& ray_ori, 
                           const Matrix<double, 3, 1>& ray_dir,
                           double* t0) = 0;

    double transparency_, reflection_;
    Matrix<double, 3, 1> surface_color_, emission_color_;
    Matrix<double, 3, 1> center_;
};


typedef std::unique_ptr<Object> ObjectUniquePtr;
typedef std::shared_ptr<Object> ObjectSharedPtr;


class Triangle: public Object {
  public:
    Triangle() {}
    Triangle(const Matrix<double, 3, 1>& pt0,
             const Matrix<double, 3, 1>& pt1,
             const Matrix<double, 3, 1>& pt2,
             const Matrix<double, 3, 1>& norm_vec)
        : pt0_(pt0), pt1_(pt1), pt2_(pt2),
          Object( 0, 0, Matrix<double, 3, 1>({{0},{0},{0}}), 
                        Matrix<double, 3, 1>({{0},{0},{0}}),
                        (pt0 + pt1 + pt2)/3)
    {
        norm_vec_ = norm_vec / norm_vec.norm();
    }


    friend std::ostream& operator << (std::ostream& os, const Triangle& tg)
    {
        os << "tg.pt0" << std::endl << tg.pt0_ << std::endl;
        os << "tg.pt1" << std::endl << tg.pt1_ << std::endl;
        os << "tg.pt2" << std::endl << tg.pt2_ << std::endl;
        os << "tg.norm_vec" << std::endl << tg.norm_vec_ << std::endl;
        return os;
    }

    virtual bool intersect(const Matrix<double, 3, 1>& ray_ori, 
                           const Matrix<double, 3, 1>& ray_dir,
                           double* t0) override 
    {
        //TODO
        return true;
    }


    Matrix<double, 3, 1> pt0_;
    Matrix<double, 3, 1> pt1_;
    Matrix<double, 3, 1> pt2_;
    Matrix<double, 3, 1> norm_vec_;
};


class Sphere: public Object {
  public:
      Sphere() {}
      Sphere(const Matrix<double, 3, 1>& center, double radius) 
          : radius_(radius), radius2_(radius * radius),
            Object( 0, 0, Matrix<double, 3, 1>({{0},{0},{0}}), 
                          Matrix<double, 3, 1>({{0},{0},{0}}),
                          center) {}

      Sphere(const Matrix<double, 3, 1>& center, double radius,
             double transparency, double reflection,
             const Matrix<double, 3, 1>& surface_color) 
          : radius_(radius), radius2_(radius * radius),
            Object( transparency, reflection, 
                    surface_color, Matrix<double, 3, 1>({{0},{0},{0}}),
                    center) {}

      Sphere(const Matrix<double, 3, 1>& center, double radius,
             double transparency, double reflection,
             const Matrix<double, 3, 1>& surface_color, 
             const Matrix<double, 3, 1>& emission_color) 
          : radius_(radius), radius2_(radius * radius),
            Object( transparency, reflection, 
                    surface_color, emission_color,
                    center) {}

      virtual bool intersect(const Matrix<double, 3, 1>& ray_ori, 
                             const Matrix<double, 3, 1>& ray_dir,
                             double* t0) override 
      {
          Matrix<double, 3, 1> l = center_ - ray_ori;
          double t_mid = ( l.transpose() * ray_dir )(0, 0);
          if (t_mid < 0) return false;
          double d2 = (l.transpose() * l)(0, 0) - t_mid * t_mid;
          if (d2 > radius2_) return false;
          *t0 = t_mid - sqrt(radius2_ - d2);
          return true;
      }


      double radius_, radius2_;

};

void generate_ray(const vec2& uv, const vec2& wh, 
                  Matrix<double, 3, 1>* ray_ori,
                  Matrix<double, 3, 1>* ray_dir)
{
    double h_meter = 0.09, w_meter = 0.16;
    double t = h_meter/2, d = -h_meter/2;
    double l = -w_meter/2, r = w_meter/2;
    double f = 100, n = 0.1;

    Matrix<double, 2, 2> S {
        {{2.0/wh.x, 0.0},
         {0.0, -2.0/wh.y}}
    };

    Matrix<double, 4, 1> pt_n_ndc {
        {{uv.x}, {uv.y}, {-1.0}, {1.0}}
    };

    pt_n_ndc.block<2, 1>(0, 0) = S * pt_n_ndc.block<2, 1>(0, 0) + Matrix<double, 2, 1>({{-(1-1/wh.x)},{1-1/wh.y}});

    Matrix<double, 4, 4> P_ndc2cam {
        {{r/n, 0.0, 0.0, 0.0},
         {0.0, t/n, 0.0, 0.0},
         {0.0, 0.0, 0.0, -1.0},
         {0.0, 0.0, (n-f)/(2*f*n), (f+n)/(2*f*n)}}
    };

    Matrix<double, 4, 1> pt_n_cam = P_ndc2cam * pt_n_ndc;
    pt_n_cam = pt_n_cam / pt_n_cam(3, 0); 
    // std::cout << "(" << uv.x << "," << uv.y << "): "<< pt_n_cam.transpose();

    (*ray_ori) = Matrix<double, 3, 1>({{0},{0},{0}});
    (*ray_dir) = pt_n_cam.block<3, 1>(0, 0);
    (*ray_dir) = (*ray_dir) / (*ray_dir).norm();
}

vec3 trace(const Matrix<double, 3, 1>& ray_ori, 
           const Matrix<double, 3, 1>& ray_dir,
           const std::vector<ObjectSharedPtr>& objs)
{
    Matrix<double, 3, 1> surface_color;
    double bias = 1e-4;

    // get hit point
    ObjectSharedPtr obj_ptr = nullptr;
    double t_min = INF;
    double t0;
    for (int i = 0; i < objs.size(); i++) {
        if ( objs[i]->intersect(ray_ori, ray_dir, &t0) ) {
            if (t0 < t_min) {
                t_min = t0;
                obj_ptr = objs[i];
            }
        }
    }

    if ( obj_ptr == nullptr ) return vec3(0.2, 0.2, 0.2);

    // light or shadow
    // std::cout << "pos: " << obj_ptr->center_.transpose() << std::endl;

    Matrix<double, 3, 1> pt_hit = ray_ori + ray_dir * t_min; 
    Matrix<double, 3, 1> norm_hit = pt_hit - obj_ptr->center_;
    norm_hit = norm_hit / norm_hit.norm();

    for (auto light: objs) {
        if ( light->emission_color_(0, 0) > 0 ) {
            // it's light

            Matrix<double, 3, 1> shadow_ray_dir = light->center_ - pt_hit;
            shadow_ray_dir = shadow_ray_dir / shadow_ray_dir.norm();

            // is light or shadow
            bool is_in_shadow = false;
            for (int i = 0; i < objs.size(); i++) {
                if (objs[i] != light) {
                    double t0_tmp;
                    if ( objs[i]->intersect(pt_hit + norm_hit*bias, shadow_ray_dir, &t0_tmp) ) {
                        is_in_shadow = true;
                        break;
                    }
                }
            }
            if (is_in_shadow) continue;

            double ratio = std::max(0.0, (norm_hit.transpose()*shadow_ray_dir)(0,0) );
            // std::cout << "ratio: " << ratio << std::endl;
            Matrix<double, 3, 1> light_strength = light->emission_color_ * ratio;
            surface_color = surface_color + 
                            Matrix<double, 3, 1>({{light_strength(0,0)*obj_ptr->surface_color_(0,0)},
                                                  {light_strength(1,0)*obj_ptr->surface_color_(1,0)},
                                                  {light_strength(2,0)*obj_ptr->surface_color_(2,0)}});
        }
    }
    surface_color = surface_color + obj_ptr->emission_color_;

    return vec3(surface_color(0, 0), surface_color(1, 0), surface_color(2, 0));
}


int main(int argc, char**argv)
{
    Triangle tg{ 
        Matrix<double, 3, 1>({{1.0}, {0.0}, {0.0}}),
        Matrix<double, 3, 1>({{-1.0}, {0.0}, {0.0}}),
        Matrix<double, 3, 1>({{0.0}, {1.0}, {0.0}}),
        Matrix<double, 3, 1>({{0.0}, {0.0}, {1.0}}),
    };

    // std::cout << tg;

    std::vector<ObjectSharedPtr> objs;
    // center  radius  transparency  reflection  surface_color  emission_color
    objs.push_back( std::make_shared<Sphere>(
                        Matrix<double, 3, 1>({{0.0},{-10004},{-20}}), 
                        10000, 0.0, 0.0,
                        Matrix<double, 3, 1>({{0.2},{0.4},{0.2}})
                    ) 
                  );

    objs.push_back( std::make_shared<Sphere>(
                        Matrix<double, 3, 1>({{0.0},{0.0},{-20}}), 
                        4.0, 0.0, 0.0,
                        Matrix<double, 3, 1>({{1.0},{0.32},{0.36}})
                    ) 
                  );

    objs.push_back( std::make_shared<Sphere>(
                        Matrix<double, 3, 1>({{5.0},{-1.0},{-15}}), 
                        2.0, 0.0, 0.0,
                        Matrix<double, 3, 1>({{0.9},{0.76},{0.46}})
                    ) 
                  );

    objs.push_back( std::make_shared<Sphere>(
                        Matrix<double, 3, 1>({{5.0},{0.0},{-25}}), 
                        3.0, 0.0, 0.0,
                        Matrix<double, 3, 1>({{0.65},{0.77},{0.97}})
                    ) 
                  );

    objs.push_back( std::make_shared<Sphere>(
                        Matrix<double, 3, 1>({{-5.5},{0.0},{-15}}), 
                        3.0, 0.0, 0.0,
                        Matrix<double, 3, 1>({{0.9},{0.9},{0.9}})
                    ) 
                  );
    // light
    objs.push_back( std::make_shared<Sphere>(
                        Matrix<double, 3, 1>({{0.0},{20.0},{-5}}), 
                        1.0, 0.0, 0.0,
                        Matrix<double, 3, 1>({{0.0},{0.0},{0.0}}), 
                        Matrix<double, 3, 1>({{1.0},{1.0},{1.0}})
                    ) 
                  );


    
    auto sd = [&](double t, int j, 
                  const vec2& r,
                  const vec2& FC,
                  const vec3& rgb,
                  vec4& o)
    {

        try {
            Matrix<double, 3, 1> ray_ori;
            Matrix<double, 3, 1> ray_dir;
            vec3 sc;
            generate_ray(vec2(FC.x, FC.y), vec2(r.x, r.y), &ray_ori, &ray_dir); // fk, vec2(r.x, r.y)
            sc = trace(ray_ori, ray_dir, objs);
            o = vec4(sc.x, sc.y, sc.z, 0.0);
            // std::cout << "(" << FC.x << "," << FC.y << "): "<< ray_dir.transpose();
        }
        catch (const char* e) {
            std::cout << e << std::endl;
        }
    };

    // TimeRange  FrameSize  Width  Height
    // fragShader(1.0, 1, 16, 9, sd);
    fragShader(1.0, 1, 16*60, 9*60, sd);
    // fragShader(1.0, 60, 16*20, 9*20, sd);

    return 0;
}

