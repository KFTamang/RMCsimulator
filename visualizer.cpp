#include <iostream>
#include <TH1.h>
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
        TH1* hque[N];
        int que_current;
        std::string opt[N];
    public:
        visualizer():
            current(0),
            que_current(0)
            {
                c1 = new TCanvas("c1","c1",600,500);
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
            hque[que_current] = new TH1D(Form("h%d",que_current),Form("h%d",que_current),N,min-deltaX,max+deltaX);
            for(int i=0;i<N;++i){
                std::cout << i << " " <<  x[i] << ":" << y[i] << std::endl;
                hque[que_current]->Fill(x[i],y[i]);
            }
            opt[que_current] = "";
//            hque[que_current]->Draw();
            ++que_current;
            std::cout << "que_current:" << que_current << std::endl;
            return que_current-1;
        }

//         int add2d(double* x,double* y,double* z,int Nx,int Ny){
//             std::cout << "add2d" << std::endl;
//             double minX = x[0];
//             double maxX = x[Nx-1];
//             double deltaX = (x[Nx-1]-x[0])/(Nx-1)/2;
//             double minY = y[0];
//             double maxY = y[Ny-1];
//             double deltaY = (y[Ny-1]-y[0])/(Ny-1)/2;
//             hque[que_current] = new TH2D(Form("h%d",que_current),Form("h%d",que_current),
//                                     Nx,minX-deltaX,maxX+deltaX,Ny,minY-deltaY,maxY+deltaY);
//             for(int i=0;i<Nx;++i){
//                 for(int m=0;m<Ny;++m){
//                     std::cout << i << ":" << m << " " << x[i] << ":" << y[i] << ":" << z[m] << std::endl;
//                     hque[que_current]->Fill(x[i],y[m],z[i+m*Nx]);
//                 }
//             }
//             opt[que_current] = "colz";
// //            hque[que_current]->Draw("colz");
//             ++que_current;
//             return que_current-1;
//         }
        int showQue(){
            c1->Clear();
            c1->Divide(1,que_current);
            for(int i=0;i<que_current;++i){
                c1->cd(i+1);
                hque[i]->Draw(opt[i].data());
            }
            return que_current;
        }
        bool saveCanvas(std::string name){
            c1->SaveAs(name.data());
            return true;
        }
};

void testVisualizer(){
    visualizer vis;
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
    vis.add1d(x,z+20,10);
    vis.showQue();
    vis.saveCanvas("vis_test.pdf");
}

int main(){
    testVisualizer();
    return 0;
}

