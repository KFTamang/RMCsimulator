#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <TH1.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

class source{
private:
  double x,y,z;
public:
  source(double _x, double _y,double _z):
    x(_x),y(_y),z(_z)
    {}
};



class ray{
private:
  source sor;
  double pos_x = 0.013;
  double pos_y = 0.0;
  double pos_z = 1000.0;
  double dir_x = 0.0015;
  double dir_y = 0.0;
  double dir_z = -1.0;
  TRandom3 rnd;
public:
  ray(source _sor):
    sor(_sor)
    {}

  void gen_new_ray();
  void gen_new_ray(double);
  void gen_new_ray_parallel(double);
  void get_cross_point_zplane(double*, double);
};
void ray::gen_new_ray(){
  // pos_x =  sor.get_sourse_pos_x();
  // pos_y =  sor.get_sourse_pos_y();
  // pos_z =  sor.get_sourse_pos_z();
  dir_x =  rnd.Uniform(0,1.0);
  dir_y =  rnd.Uniform(0,1.0);
  dir_z =  rnd.Uniform(0,1.0);
}
void ray::gen_new_ray(double cos_max){
  // pos_x =  sor.get_sourse_pos_x();
  // pos_y =  sor.get_sourse_pos_y();
  // pos_z =  sor.get_sourse_pos_z();
  if( cos_max < -1.0 ){
    std::cout << "ERROR:cos_max is too small" << std::endl;
    return;
  }else if( cos_max>1.0 ){
    std::cout << "ERROR:cos_max is too large" << std::endl;
    return;
  }

  dir_z =  rnd.Uniform(-1.0,cos_max);
  double r = pow(1.0 - dir_z*dir_z,0.5);
  double theta = rnd.Uniform(0,2*M_PI);
  dir_x =  r*cos(theta);
  dir_y =  r*sin(theta);

  // for debug
//  std::cout << dir_x << ":" << dir_y << ":" << dir_z << " " << r << " " << theta << std::endl;
}

void ray::gen_new_ray_parallel(double radius){
  pos_x =  rnd.Uniform(-2*radius,2*radius);
  pos_y =  rnd.Uniform(-2*radius,2*radius);

  // for debug
//  std::cout << dir_x << ":" << dir_y << ":" << dir_z << " " << r << " " << theta << std::endl;
}

// calculating the positin of the cross point of the z=z_0 plane and the strainght line along the x-ray
void ray::get_cross_point_zplane(double* pos, double z){
  double a = (z-pos_z)/dir_z;
  pos[0] = a*dir_x + pos_x;  
  pos[1] = a*dir_y + pos_y;  
  pos[2] = z;
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
    static constexpr int num_strips = 32;
    static constexpr double width_strip = 0.15;
    static constexpr double width_collimator = 4.5;
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
          rotate_2d(cp_pos,angle);
          image->Fill(cp_pos[0],cp_pos[1]);
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



int main(){

  
 rotationModulationCollimator col0(0.1);
 rotationModulationCollimator col1(155.0);
  // collimator col0(1);
  // collimator col1(2);

  source sor(0,0,1000.0);
  ray ray0(sor);

  detector det(0.0);

    col0.clock_collimator(0);
    col1.clock_collimator(0);

  TH1D* h = new TH1D("h","h",360,-0.5,360-0.5);

  const  int iteration = 100000;
  for(int m=0;m<360;++m){
    int cnt = 0;
    col0.clock_collimator(m);
    col1.clock_collimator(m);
    for(int i=0;i<iteration;++i){
//      ray0.gen_new_ray(-pow(100.0/101.0,0.5));
      ray0.gen_new_ray_parallel(det.getRadius());
//      ray0.gen_new_ray(-1.0+0.045*0.045/2);
      if( col0.blocks(ray0) || col1.blocks(ray0) || !det.detects(ray0) ){
        continue;
      }else{
//        std::cout << i << " hit" << std::endl;
        ++cnt;
      }
    }
    std::cout << m << " " << cnt << std::endl;
    h->Fill(m,1.0*cnt/iteration);
  }
  

  // debug
  TCanvas* c1 = new TCanvas();
  h->Draw("hist");
  c1->SaveAs("colli_image.pdf(");
  col0.image->Draw("colz");
  c1->SaveAs("colli_image.pdf");
  col1.image->Draw("colz");
  c1->SaveAs("colli_image.pdf");
  det.image->Draw("colz");
  c1->SaveAs("colli_image.pdf)");
  c1->Close();

  return 0;
}
