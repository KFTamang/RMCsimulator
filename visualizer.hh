#include <iostream>
#include <TRandom3.h>
//#include <TH.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>


class visualizer {
    private:
        static constexpr int N = 10;
        TH1D* h1[N];
        TH2D* h2[N];
        int current;
        TCanvas* c1;
        TH1* hque_1d[N];
        TH2* hque_2d[N];
        int que_current_1d;
        int que_current_2d;
        std::string opt[2*N];
    public:
        visualizer():
            current(0),
            que_current_1d(0),
            que_current_2d(0)
            {
                c1 = new TCanvas("c1","c1",400,600);
                std::cout << "Initialize visualizer instance" << std::endl;
                std::cout << "current : " << current << std::endl;
            }
        visualizer(int width, int height):
            current(0),
            que_current_1d(0),
            que_current_2d(0)
            {
                c1 = new TCanvas("c1","c1",width,height);
                std::cout << "Initialize visualizer instance" << std::endl;
                std::cout << "current : " << current << std::endl;
            }
        ~visualizer(){}
        bool plot1d(double* x,double* y,int N){
            std::cout << "plot1d" << std::endl;
            double min = x[0];
            double max = x[N-1];
            double deltaX = (x[N-1]-x[0])/(N-1)/2;
            h1[current] = new TH1D(Form("h%d",current),Form("h%d",current),N,min-deltaX,max+deltaX);
            for(int i=0;i<N;++i){
                std::cout << i << " " <<  x[i] << ":" << y[i] << std::endl;
                h1[current]->Fill(x[i],y[i]);
            }
            h1[current]->Draw();
            ++current;
            return true;
        }

        bool plot2d(double* x,double* y,double* z,int Nx,int Ny){
            std::cout << "plot2d" << std::endl;
            double minX = x[0];
            double maxX = x[Nx-1];
            double deltaX = (x[Nx-1]-x[0])/(Nx-1)/2;
            double minY = y[0];
            double maxY = y[Ny-1];
            double deltaY = (y[Ny-1]-y[0])/(Ny-1)/2;
            h2[current] = new TH2D(Form("h%d",current),Form("h%d",current),
                                    Nx,minX-deltaX,maxX+deltaX,Ny,minY-deltaY,maxY+deltaY);
            for(int i=0;i<Nx;++i){
                for(int m=0;m<Ny;++m){
                    std::cout << i << ":" << m << " " << x[i] << ":" << y[i] << ":" << z[m] << std::endl;
                    h2[current]->Fill(x[i],y[m],z[i+m*Nx]);
                }
            }
            h2[current]->Draw("colz");
            ++current;
            return true;
        }
        int add1d(double* x,double* y,int N){
            std::cout << "add1d" << std::endl;
            double min = x[0];
            double max = x[N-1];
            double deltaX = (x[N-1]-x[0])/(N-1)/2;
            hque_1d[que_current_1d] = new TH1D(Form("h%d_1d",que_current_1d),Form("h%d_1d",que_current_1d),N,min-deltaX,max+deltaX);
            for(int i=0;i<N;++i){
                std::cout << i << " " <<  x[i] << ":" << y[i] << std::endl;
                hque_1d[que_current_1d]->Fill(x[i],y[i]);
            }
            opt[que_current_1d] = "";
//            hque_1d[que_current_1d]->Draw();
            ++que_current_1d;
            std::cout << "que_current_1d:" << que_current_1d << std::endl;
            return que_current_1d-1;
        }

        int add2d(double* x,double* y,double* z,int Nx,int Ny){
            std::cout << "add2d" << std::endl;
            double minX = x[0];
            double maxX = x[Nx-1];
            double deltaX = (x[Nx-1]-x[0])/(Nx-1)/2;
            double minY = y[0];
            double maxY = y[Ny-1];
            double deltaY = (y[Ny-1]-y[0])/(Ny-1)/2;
            hque_2d[que_current_2d] = new TH2D(Form("h%d_2d",que_current_2d),Form("h%d_2d",que_current_2d),
                                    Nx,minX-deltaX,maxX+deltaX,Ny,minY-deltaY,maxY+deltaY);
            for(int i=0;i<Nx;++i){
                for(int m=0;m<Ny;++m){
                    std::cout << i << ":" << m << " " << x[i] << ":" << y[i] << ":" << z[m] << std::endl;
                    hque_2d[que_current_2d]->Fill(x[i],y[m],z[i+m*Nx]);
                }
            }
            opt[que_current_2d+N] = "colz";
//            hque_2d[que_current_2d]->Draw("colz");
            ++que_current_2d;
            return que_current_2d-1;
        }
        int add1dhist(TH1D* h){
            std::cout << "add1dhist" << std::endl;
            hque_1d[que_current_1d] = h;
            opt[que_current_1d] = "";
//            hque_1d[que_current_1d]->Draw();
            ++que_current_1d;
            std::cout << "que_current_1d:" << que_current_1d << std::endl;
            return que_current_1d-1;
        }

        int add2dhist(TH2D* h){
            std::cout << "add2dhist" << std::endl;
            hque_2d[que_current_2d] = h;
            opt[que_current_2d+N] = "colz";
//            hque_2d[que_current_2d]->Draw("colz");
            ++que_current_2d;
            return que_current_2d-1;
        }
        int showQue(){
            c1->Clear();
            c1->Divide(1,que_current_1d+que_current_2d);
            for(int i=0;i<que_current_1d;++i){
                c1->cd(i+1);
                hque_1d[i]->Draw(opt[i].data());
            }
            for(int i=0;i<que_current_2d;++i){
                c1->cd(i+1+que_current_1d);
                hque_2d[i]->Draw(opt[i+N].data());
            }
            return que_current_1d+que_current_2d;
        }
        bool saveCanvas(std::string name){
            c1->SaveAs(name.data());
            return true;
        }
};

void testVisualizer(){
    visualizer vis(400,1200);
    double x[] = {0,1,2,3,4,5,6,7,8,9};
    double y[] = {0,1,2,3,4,5,6,7,8,9};
    double z[] =   {3,4,5,6,5,4,3,2,1,0,
                    4,5,6,5,4,3,2,1,0,3,
                    5,6,5,4,3,2,1,0,3,4,
                    6,5,4,3,2,1,0,3,4,5,
                    5,4,3,2,1,0,3,4,5,6,
                    4,3,2,1,0,3,4,5,6,5,
                    3,2,1,0,3,4,5,6,5,4,
                    2,1,0,3,4,5,6,5,4,3,
                    1,0,3,4,5,6,5,4,3,2,
                    0,3,4,5,6,5,4,3,2,1};
//    vis.plot1d(x,y,10);
//    vis.plot2d(x,y,z,10,10);
    vis.add1d(x,y,10);
    vis.add1d(x,z+10,10);
    vis.add2d(x,y,z,10,10);
    TH2D* h = new TH2D("h","h",100,0,100,100,0,100);
    h->Fill(50,50,20);
    vis.add2dhist(h);
    vis.showQue();
    vis.saveCanvas("vis_test.pdf");
}



