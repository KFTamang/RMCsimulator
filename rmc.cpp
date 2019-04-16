#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>

#include "visualizer.hh"

#include <TH1.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

const bool debug = true;
// const bool debug = false;

class source{
private:
  double x,y,z;
public:
  source(double _x, double _y,double _z):
    x(_x),y(_y),z(_z)
    {}
  double get_sourse_pos_x(){return x;};
  double get_sourse_pos_y(){return y;};
  double get_sourse_pos_z(){return z;};

};



class ray{
private:
  source sor;
  double pos_x = 0.03;
  double pos_y = 0.0;
  double pos_z = 1000.0;
  double dir_x = 0.0015;
  double dir_y = 0.0;
  double dir_z = -1.0;
  //  TRandom3 rnd;
  std::mt19937 mt;
  std::uniform_real_distribution<double> rnd;
public:
  ray(source _sor):
    sor(_sor),
    mt(100),
    rnd(0.0,1.0)
    {}

  void gen_new_ray();
  void gen_new_ray(double);
  void gen_new_ray_parallel(double);
  void gen_new_ray_parallel_withangle(double);
  void get_cross_point_zplane(double*, double);
  void set_pos(double* _pos){
    pos_x = _pos[0];
    pos_y = _pos[1];
    pos_z = _pos[2];
  };
  void set_pos_with_random_width(double* _pos,double width){
    pos_x = _pos[0] + width*(2*rnd(mt)-1.0);
    pos_y = _pos[1] + width*(2*rnd(mt)-1.0);
  };
  void set_dir(double* dir){
    dir_x = dir[0];
    dir_y = dir[1];
    dir_z = dir[2];
  };
  void calc_pos(){};
  void print_dir(){
    std::cout << "dir x:" << dir_x << ", y:" << dir_y << ", z:" << dir_z << std::endl; 
  }
  void print_pos(){
    std::cout << "pos x:" << pos_x << ", y:" << pos_y << ", z:" << pos_z << std::endl; 
  }
};
void ray::gen_new_ray(){
  pos_x =  sor.get_sourse_pos_x();
  pos_y =  sor.get_sourse_pos_y();
  pos_z =  sor.get_sourse_pos_z();
  // dir_x =  rnd.Uniform(0,1.0);
  // dir_y =  rnd.Uniform(0,1.0);
  // dir_z =  rnd.Uniform(0,1.0);
  dir_x =  rnd(mt);
  dir_y =  rnd(mt);
  dir_z =  rnd(mt);
}
void ray::gen_new_ray(double cos_max){
  pos_x =  sor.get_sourse_pos_x();
  pos_y =  sor.get_sourse_pos_y();
  pos_z =  sor.get_sourse_pos_z();
  if( cos_max < -1.0 ){
    std::cout << "ERROR:cos_max is too small" << std::endl;
    return;
  }else if( cos_max>1.0 ){
    std::cout << "ERROR:cos_max is too large" << std::endl;
    return;
  }

  dir_z =  (cos_max+1)*rnd(mt)-1; //rnd.Uniform(-1.0,cos_max);
  double r = pow(1.0 - dir_z*dir_z,0.5);
  double theta = 2.*M_PI*rnd(mt); //rnd.Uniform(0,2*M_PI);
  dir_x =  r*cos(theta);
  dir_y =  r*sin(theta);

  // for debug
//  std::cout << dir_x << ":" << dir_y << ":" << dir_z << " " << r << " " << theta << std::endl;
}

void ray::gen_new_ray_parallel(double radius){
  pos_x = 4.0*radius*(rnd(mt)-0.5); //rnd.Uniform(-2*radius,2*radius);
  pos_y = 4.0*radius*(rnd(mt)-0.5); //rnd.Uniform(-2*radius,2*radius);

  // for debug
//  std::cout << dir_x << ":" << dir_y << ":" << dir_z << " " << r << " " << theta << std::endl;
}

void ray::gen_new_ray_parallel_withangle(double angle ){
  // double radius = 1.0;
  // pos_x = 4.0*radius*(rnd(mt)-0.5); //rnd.Uniform(-2*radius,2*radius);
  // pos_y = 4.0*radius*(rnd(mt)-0.5); //rnd.Uniform(-2*radius,2*radius);
  dir_x = angle;
  dir_y = 0.0;
  dir_z = -1.0;
  // for debug
//  std::cout << dir_x << ":" << dir_y << ":" << dir_z << " " << r << " " << theta << std::endl;
}



// calculating the positin of the cross point of the z=z_0 plane and the strainght line along the x-ray
void ray::get_cross_point_zplane(double* pos, double z){
  double a = (z-pos_z)/dir_z;
  pos[0] = a*dir_x + pos_x;  
  pos[1] = a*dir_y + pos_y;  
  pos[2] = z;
  if(debug){
    std::cout << "a:" <<  a  << ", pos[0]:" << pos[0];
    std::cout << ", pos[1]:" << pos[1] << ", pos[2]:" << pos[2] << std::endl;
  } 
}




class collimator{
  private:

  protected:
    double pos;

  public:
    collimator(double _pos):
      pos(_pos)
      {};
    virtual bool blocks(ray ray0){ // returns true if the collimator blocks the line of sight of the x-ray
      return false;
    }
};

class rotationModulationCollimator : public collimator{
  private:
    static constexpr int num_strips = 5; //32;
    static constexpr double width_collimator = 4.5;
    static constexpr double width_strip = width_collimator/num_strips/2 ; // 0.15;
    static constexpr double width_gap = (2.0*width_collimator - num_strips * width_strip)/(num_strips-1);
//    static constexpr double width_gap = 0.02;
//    static constexpr int num_strips = ;

    double angle = 0.0; // in degree
    void rotate_2d(double* pos, double angle){
      double rotated_pos0 =  pos[0]*cos(angle)-pos[1]*sin(angle);
      double rotated_pos1 =  pos[0]*sin(angle)+pos[1]*cos(angle);
      pos[0] = rotated_pos0;
      pos[1] = rotated_pos1;
    }

  public:
    double colli_image_x[100];
    double colli_image_y[100];
    double colli_image_z[100*100];
    TH2D* image = new TH2D("colli_image","colli_image",400,-width_collimator,width_collimator,400,-width_collimator,width_collimator); // debug

    rotationModulationCollimator(double _pos): collimator(_pos){
      //dubug
      std::cout << "width_strip:" << width_strip << " width_collimator:" << width_collimator << " width_gap:" << width_gap << std::endl;
    }

    bool clock_collimator(double clock_angle){
      if(clock_angle>=0.0 && clock_angle<360){
        angle = clock_angle;
        return true;
      }else{
        return false;
      }
    }

    virtual bool blocks(ray ray0){ // returns true if the RMC blocks the line of sight of the x-ray
      double cp_pos[3];
      ray0.get_cross_point_zplane(cp_pos,pos);
      rotate_2d(cp_pos,-angle);
//      std::cout << cp_pos[0] << " " << cp_pos[1] << " " << cp_pos[2] << std::endl; //debug
      for(int i=0;i<num_strips;++i){
        double strip_left =      i*width_strip + i*width_gap - width_collimator;
        double strip_right = (i+1)*width_strip + i*width_gap - width_collimator;
//        std::cout << i << " " << strip_left << " " << strip_right << std::endl; //debug
        if(cp_pos[0]>strip_left && cp_pos[0]<strip_right && cp_pos[1]<width_collimator && cp_pos[1]>-width_collimator){
          image->Fill(cp_pos[0],cp_pos[1]);
          rotate_2d(cp_pos,angle);
          return true;
        }
      }
      return false;
    }


};




class detector{
  private:
    double pos;
    double radius = 9.0/2;
  public:
    detector(double _pos):
      pos(_pos)
      {}
    bool detects(ray);
    double getRadius(){return radius;}

    TH2D* image = new TH2D("det_image","det_image",100,-radius,radius,100,-radius,radius);

};
bool detector::detects(ray ray0){
  double point[3];
  ray0.get_cross_point_zplane(point, pos);
  double dist = pow(point[0]*point[0]+point[1]*point[1],0.5);
//  std::cout << point[0] << ":" << point[1] << ":" << point[2] << " " << dist << std::endl;
  if(dist<radius){
    image->Fill(point[0],point[1]);
    return true; 
  }else{
    return false;
  }
}


class rmc_system {
  public:
    rmc_system ();
    ~rmc_system ();
    void generate_response();
    static constexpr int N  = 1000;
    double response_x[N];
    double response_y[N];
    bool run_single_step_parallel(){
      ray0.gen_new_ray_parallel(det.getRadius());
      if( col0.blocks(ray0) || col1.blocks(ray0) || !det.detects(ray0) ){
        return false;
      }else{
        return true;
      }
    };
    bool run_single_step_parallel_with_angle(double angle){
      ray0.gen_new_ray_parallel_withangle(angle);
      if( col0.blocks(ray0) || col1.blocks(ray0) || !det.detects(ray0) ){
        return false;
      }else{
        return true;
      }
    };
    void set_source_pos(){
      double pos[3] = {};
      ray0.get_cross_point_zplane(pos,1000.0);
      ray0.set_pos(pos);
    };
    void set_source_pos_with_width(double width){
      double pos[3] = {};
      ray0.get_cross_point_zplane(pos,1000.0);
      ray0.set_pos_with_random_width(pos,width);
    };
    void set_source_pos_with_angle(double angle){
      double dir[3] = {angle,0.0,-1.0};
      ray0.set_dir(dir);
      set_source_pos();
    };
    void set_source_pos_randomly_with_angle(double angle,double source_width){
      double dir[3] = {angle,0.0,-1.0};
      ray0.set_dir(dir);
      set_source_pos_with_width(source_width);
    };
    void set_source_pos_randomly(double source_width){
      set_source_pos_with_width(source_width);
    };
    void clock_collimator(){
      ++collimator_angle;
      col0.clock_collimator(collimator_angle);
      col1.clock_collimator(collimator_angle);

      };

    rotationModulationCollimator get_col0(){
      return col0;
    };
    rotationModulationCollimator get_col1(){
      return col1;
    };

  private:
    rotationModulationCollimator col0;
    rotationModulationCollimator col1;
    detector det;
    source sor;
    ray ray0;
    int collimator_angle;



};

rmc_system ::rmc_system () : 
  col0(0.1),
  col1(155.0),
  det(0.0),
  sor(0.03,0.0,1000.0),
  ray0(sor),
  collimator_angle(0) {
}

rmc_system ::~rmc_system () {
}

void rmc_system::generate_response(){
  int iter = 1;
  for(int i=0;i<N;++i){
    int count = 0;
    double angle  = 0.2*i/N;
    set_source_pos_with_angle(angle);
    if(debug)    ray0.print_dir();
    if(debug)    ray0.print_pos();
    for(int m=0;m<iter;++m){
      if(run_single_step_parallel_with_angle(angle)){
        ++count;
      }
      set_source_pos_randomly(det.getRadius());
      if(debug) ray0.print_pos();
      if(debug) ray0.print_dir();
      if(debug) std::cout << std::endl;
      response_x[i] = angle;
      response_y[i] = 1.0*count/iter;

    }
  }
}

// class image_reconstructor {
//   public:
//     image_reconstructor(collimator,collimator,detector);
//     ~image_reconstructor();

//   private:
//     double response[1000][2];  
// };

// image_reconstructor::image_reconstructor(collimator col0,collimator col1, detector det) {
//   source sor(0.003,0.0,1000.0);
//   ray ray0(sor);
  
// }

// image_reconstructor::~image_reconstructor() {
// }


int main(){

  visualizer vis(400,800);

  rmc_system rmc_system0;

//  TH1D* h = new TH1D("h","h",360,-0.5,360-0.5);
  double count_rate_data_x[360] = {};
  double count_rate_data_y[360] = {};
  const  int iteration = 10;
//   for(int m=0;m<360;++m){
//     int cnt = 0;
//     for(int i=0;i<iteration;++i){
//       if(rmc_system0.run_single_step_parallel()){
//         ++cnt;
//       }
//     }
//     // rmc_system0.clock_collimator();
//     std::cout << m << " " << cnt << std::endl;
// //    h->Fill(m,1.0*cnt/iteration);
//     count_rate_data_x[m] = m;
//     count_rate_data_y[m] = 1.0*cnt/iteration;

//   }
  
  // vis.add1d(count_rate_data_x,count_rate_data_y,360);

  vis.add2dhist(rmc_system0.get_col0().image);

  rmc_system0.generate_response();
  vis.add1d(rmc_system0.response_x,rmc_system0.response_y,rmc_system0.N);



  // debug
//   TCanvas* c1 = new TCanvas();
// //  h->Draw("hist");
//   c1->SaveAs("colli_image.pdf(");
//   col0.image->Draw("colz");
//   c1->SaveAs("colli_image.pdf");
//   col1.image->Draw("colz");
//   c1->SaveAs("colli_image.pdf");
//   det.image->Draw("colz");
//   c1->SaveAs("colli_image.pdf)");
//   c1->Close();

  vis.showQue();
  vis.saveCanvas("rmc_image.pdf");
  return 0;
}
