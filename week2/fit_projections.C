Int_t numbins = 0;

Double_t gauss(Double_t x,Double_t *par)
{
  Double_t value = par[0]*TMath::Gaus(x,par[1],par[2])+par[3];
  return value;  
}

//  Opening file

TFile *nfile = new TFile("Projections.root","read");
TNtuple* my_tuple;
    
Float_t minimize_func(Double_t *par)
{   
    Float_t BinY[400], NY[400], errY[400], Ne[400] ;
    
    Float_t chi2 = 0;
    
    Float_t *row_content;
    
    Int_t numbins = 0;
  

    for (int i = 0; i<my_tuple->GetEntries(); i++)
    {
       my_tuple->GetEntry(i);
       row_content = my_tuple->GetArgs();
       
       BinY[i] = row_content[3];
       NY[i] = row_content[4];
       errY[i] = row_content[5];
       Ne[i] = 0;
       
       
       if (NY[i] > 0 && errY[i] > 0 && BinY[i] > 0.68 && BinY[i] < 0.95)
       {
          Float_t dx = (gauss(BinY[i],par) - NY[i])/errY[i];
	  std::cout << par[0] << std::endl;
          chi2 += (dx*dx);
          numbins++;
       }
       
    }
    
    std::cout << chi2 << std::endl;
    return chi2 ;
   
}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   f = minimize_func(par);
}

int numpar = 3;
double arglist[4] = {8.45083e+02,7.36461e-01,7.43977e-02,3.25635e+02};


void myminimizer()
{

    TMinuit *gMinuit2 = new TMinuit(10);
    gMinuit2->SetFCN(fcn);
    Int_t ierflg = 0 ;
    gMinuit2->mnexcm("SET ERR", arglist ,1,ierflg);

    gMinuit2->mnparm(0, "N", arglist[0], 0.0001, 800, 1500, ierflg);
    gMinuit2->mnparm(1, "mu", arglist[1], 0.0001,0.6, 1, ierflg);
    gMinuit2->mnparm(2, "sigma",arglist[2], 0.0001, 0.02, 0.3, ierflg);
    gMinuit2->mnparm(3, "c",arglist[3], 0.0001, 50, 400, ierflg);

    gMinuit2->mnexcm("MINOS", arglist , 2,ierflg);
}


void fit_projections()
{

    nfile->GetObject("tuple1", my_tuple);

    Float_t BinY[400], NY[400], errY[400], Ne[400] ;
    
    Float_t chi2 = 0;
    
    Float_t *row_content;


    for (int i = 0; i<my_tuple->GetEntries(); i++)
    {
         my_tuple->GetEntry(i);
         row_content = my_tuple->GetArgs();
       
         BinY[i] = row_content[3];
         NY[i] = row_content[4];
         errY[i] = row_content[5];
         Ne[i] = 0;
     }
   
   auto gr = new TGraphErrors(200,BinY,NY,Ne,errY);
   
   gr->Draw();

   myminimizer();

   double x[800], y[800];
   double pars[4]  = {8.45083e+02,7.36461e-01,7.43977e-02,3.25635e+02};
   for (int i=0;i< 60 ;i++)
   {
       x[i] = 0.5 + i*0.01;
       y[i] = gauss(x[i],pars);
   }

    auto g = new TGraph(60,x,y);
    g->SetLineColor(2);
    g->Draw("SAME");

}
