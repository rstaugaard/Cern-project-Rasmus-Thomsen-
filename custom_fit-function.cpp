#include "TMinuit.h"

auto rng = new TRandom();
auto h1 = new TH1F("h1","Histogram",40,10,30);
   
// Fitting function   
Double_t func(float x,Double_t *par)
 {
  Double_t value = par[0]*TMath::Gaus(x,par[1],par[2]);
  return value;
 }

Int_t numbins = 0;

void fcn(Int_t &npar, double *gin, double &f, double *par, int iflag)
{
   double chi2 = 0;

   Double_t xval;
   Double_t yval;
   Double_t err;
   Double_t dx;
   numbins = 0;
   for(unsigned int i = 0 ; i < h1->GetNbinsX(); ++i)
  {
       xval = h1->GetBinCenter(i);
       yval = h1->GetBinContent(i);
       err = h1->GetBinError(i);
       if (yval > 0)
       {
          dx = (func(xval,par) - yval)/err;
          chi2 += (dx*dx);
	  numbins++;
       }
  }
       std::cout << chi2 << std::endl;     
       f = chi2 ;
}

// Starting parameters of fit 
int numpar = 3;
double arglist[3] = {10000,20,2};


void myminimizer()
{
    TMinuit *gMinuit2 = new TMinuit(10);
    gMinuit2->SetFCN(fcn);
    Int_t ierflg = 0 ;
    gMinuit2->mnexcm("SET ERR", arglist ,1,ierflg);

    gMinuit2->mnparm(0, "N", 10000, 0.1, 0, 0, ierflg);
    gMinuit2->mnparm(1, "mu", 21, 0.01,0, 0, ierflg);
    gMinuit2->mnparm(2, "sigma",    2, 0.001, 0, 0, ierflg);

    gMinuit2->mnexcm("MINOS", arglist , 2,ierflg);
}



void new_fit()
{     
   for(int i = 0; i < 10000; i++)
   {
       h1->Fill(rng->Gaus(20,2)+rng->Gaus(0,0.5));
   }

   myminimizer();
   h1->Draw();
   std::cout << numbins << std::endl;

// Draw fit on top of histogram 
   double x[400], y[400];
   double pars[3]  = {9.62532e+02,1.99917e+01,2.06660e+00};
   for (int i=12;i< 213;i++) 
   {
      x[i] = 12+i*0.1;
      y[i] = func(x[i],pars);
   }

    auto g = new TGraph(200,x,y);
    g->Draw("SAME");
}
