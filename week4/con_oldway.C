#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/old/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

TH1F *h4 = new TH1F("h4","h4", 200,0.2,1.4);
TH1F *h5 = new TH1F("h5","h5", 100,0.4,3);


     //  Handy functions     ----------------------------------------------------
Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

Float_t calc_E(Float_t px,Float_t py,Float_t pz)
{
    return TMath::Sqrt(mpi*mpi+px*px+py*py+pz*pz);
}

Double_t calc_InvM(Double_t *p1, Double_t *p2)
{
    Double_t E1 = calc_E(p1[0],p1[1],p1[2]);
    Double_t E2 = calc_E(p2[0],p2[1],p2[2]);
    
    return TMath::Sqrt(TMath::Power(E1+E2,2) - TMath::Power(p1[0]+p2[0],2)- TMath::Power(p1[1]+p2[1],2)- TMath::Power(p1[2]+p2[2],2));
}

Double_t calc_fourmass(Double_t *p1, Double_t *p2,Double_t *p3,Double_t *p4)
{
    Double_t E1 = calc_E(p1[0],p1[1],p1[2]);
    Double_t E2 = calc_E(p2[0],p2[1],p2[2]);
    Double_t E3 = calc_E(p3[0],p3[1],p3[2]);
    Double_t E4 = calc_E(p4[0],p4[1],p4[2]);
    
    return TMath::Sqrt(TMath::Power(E1+E2+E3+E4,2) - TMath::Power(p1[0]+p2[0]+p3[0]+p4[0],2)- TMath::Power(p1[1]+p2[1]+p3[1]+p4[1],2)- TMath::Power(p1[2]+p2[2]+p3[2]+p4[2],2));
}

Int_t count = 0;

     //  fit functions         ----------------------------------------------------

Double_t gauss(float x,Double_t *par)
{
    Double_t value = par[0]*TMath::Gaus(x,par[1],par[2]);
    return value;
}
 
 
Double_t background(Float_t x,Double_t *par)
{
    Double_t value = par[0]*TMath::Power(x-par[1],par[2])*TMath::Exp(-x*par[3]);
    return value;
}

Double_t totalfit(Float_t x,Double_t *par)
{
   Double_t bp[4] = { 6.06557e+04,2.55558e-01,9.25427e-01,3.92736e+00};
   //Double_t pars1[4] = {par[0],par[1],par[2],par[3]};
   Double_t pars2[3] = {par[0],par[1],par[2]};
   Double_t pars3[3] = {par[3],par[4],par[5]};
   
   Double_t value = background(x,bp)+gauss(x,pars2)+gauss(x,pars3);
   return value;
}



Int_t numbins;
 

void fcn(Int_t &npar, double *gin, double &f, double *par, int iflag)
{
   double chi2 = 0;

   Double_t xval;
   Double_t yval;
   Double_t err;
   Double_t dx;
   numbins = 0;
   for(unsigned int i = 0 ; i < h4->GetNbinsX(); ++i)
  {    
       xval = h4->GetBinCenter(i);
       yval = h4->GetBinContent(i);
       err = h4->GetBinError(i);
       if (yval > 0 && xval > 0.3)
       {
          dx = (totalfit(xval, par) - yval)/err;
          chi2 += (dx*dx);
	  numbins++;
       }
   }
       //std::cout << chi2 << std::endl;     
       f = chi2 ;
}


const Int_t numpar = 6;
Int_t num_entries = 0;

void myminimizer(Double_t *par, Double_t *err)
{
    Double_t arglist[6] = {1.34362e+03,7.49886e-01,6.04308e-02,1.25633e+03,4.97597e-01,1.80171e-02};
    
    TMinuit *gMinuit2 = new TMinuit(10);
    gMinuit2->SetFCN(fcn);
    
    Int_t ierflg = 0 ;
    gMinuit2->mnexcm("SET ERR", arglist ,1,ierflg);
    //gMinuit2->mnparm(0, "A", arglist[0], 0.0001, 0, 0, ierflg);
    //gMinuit2->mnparm(1, "B", arglist[1], 0.0001,0, 0, ierflg);
    //gMinuit2->mnparm(2, "C",arglist[2], 0.0001, 0, 0, ierflg);
    //gMinuit2->mnparm(3, "D",arglist[3], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(0, "N", arglist[0], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(1, "mu", arglist[1], 0.0001,0, 0, ierflg);
    gMinuit2->mnparm(2, "sigma",arglist[2], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(3, "N2", arglist[3], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(4, "mu2", arglist[4], 0.0001,0, 0, ierflg);
    gMinuit2->mnparm(5, "sigma2",arglist[5], 0.0001, 0, 0, ierflg);
   

    gMinuit2->mnexcm("MIGRAD", arglist , 2,ierflg);
    
    for (int j = 0; j < 6; j++)
    {
       gMinuit2->GetParameter(j,par[j],err[j]);
    }
}


//                     ----------------------------------------------------

void con_oldway()
{
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,0.2,2.5,200,0.2,2.5);
    
    TH1F *h1 = new TH1F("h1","px",200,-3,3);   
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-4,4);
    
    TH2F *d5 = new TH2F("d5","dxy / dxy" , 100,-2,2,100,-2,2);
    TH2F *d6 = new TH2F("d6","tdxy / dxy" , 100,-2,2,100,-2,2);
    
    
    static Float_t trk_pt[220], trk_eta[220], trk_phi[220], trk_dedx[220], trk_p[220];
    static Int_t trk_isK[220], trk_isPi[220], trk_isP[220];
    static Int_t trk_q[220];
    static Int_t ntrk;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q); 
    my_tree->SetBranchAddress("trk_isPi",&trk_isPi);
    my_tree->SetBranchAddress("trk_isP",&trk_isP);
    my_tree->SetBranchAddress("trk_isK",&trk_isK);
    my_tree->SetBranchAddress("trk_dedx",&trk_dedx);
    my_tree->SetBranchAddress("trk_p",&trk_p);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
 
    int numEnt = my_tree->GetEntries();

    for (int irow = 0; irow< numEnt; irow++)
    {
    
    if (mpi < 0.1)
    {
       std::cout << "Entry:" << irow << std::endl;
    }
    
    
        my_tree->GetEntry(irow);
        if (ntrk == 4)
	{
	    Float_t px[4];
            Float_t py[4];
	    Float_t pz[4];
	    Float_t pt[4];
	    Float_t q[4];
	    Float_t dedx[4];
	    Float_t eta[4];
	    Float_t phi[4];
	    Float_t p[4];
	    
	    Int_t isPi[4];
	    Int_t isP[4];
	    Int_t isK[4];
	
	    for (int i = 0; i < 4; i++)
	    {   
	        eta[i] = trk_eta[i];
		phi[i] = trk_phi[i];
		pt[i] = trk_pt[i];
	        px[i] = pt[i]*TMath::Cos(phi[i]);
	        py[i] = pt[i]*TMath::Sin(phi[i]);
	        pz[i] = pt[i]*TMath::SinH(eta[i]);
		q[i] = trk_q[i];
		dedx[i] = trk_dedx[i];
		
		isPi[i] = trk_isPi[i];
	        isP[i] = trk_isP[i];
		isK[i] = trk_isK[i];
	    }
	    
	    h1->Fill(px[0]+px[1]+px[2]+px[3]);
	    h2->Fill(py[0]+py[1]+py[2]+py[3]);
	    h3->Fill(pz[0]+pz[1]+pz[2]+pz[3]);
		
	    
	    if (q[0]+q[1]+q[2]+q[3] != 0) 
	    {
	        continue;
	    }
	    
	    Double_t p1[3] = {px[0],py[0],pz[0]};
            Double_t p2[3] = {px[1],py[1],pz[1]};
	    Double_t p3[3] = {px[2],py[2],pz[2]};
            Double_t p4[3] = {px[3],py[3],pz[3]};
	    
	    Float_t invM1, invM2,invM3, invM4;
	    
	    if (isP[0] != 2 && isP[1] != 2 && isP[2] != 2 && isP[3] != 2)
	    {
	      
	       if ((isK[0] != 2 || dedx[0] < 4) && (isK[1] != 2 || dedx[1] < 4)&& (isK[2] != 2 || dedx[2] < 4) && (isK[3] != 2 || dedx[3] < 4))
	       {
	        
	           if (dedx[0]> 0.5 && dedx[1] > 0.5 && dedx[2] > 0.5 && dedx[3] > 0.5)
	           {
		       if (q[0]+q[1] == 0)
	               {   
	                   if (q[0] + q[2] == 0)
		           {
		               invM1 = calc_InvM(p1,p2);
                               invM2 = calc_InvM(p3,p4);
		
	                       invM3 = calc_InvM(p1,p3);
                               invM4 = calc_InvM(p2,p4);

		           }
		
		           else
		           {
	                       invM1 = calc_InvM(p1,p2);
                               invM2 = calc_InvM(p3,p4);
		
		               invM3 = calc_InvM(p1,p4);
	                       invM4 = calc_InvM(p2,p3);
		            }
		
	                }
	     
	                else
	                {
		               invM1 = calc_InvM(p1,p3);
	                       invM2 = calc_InvM(p2,p4);
		
                               invM3 = calc_InvM(p1,p4);
                               invM4 = calc_InvM(p2,p3);
	                 }
		   
		   
		   
		   
	                d1->Fill(invM1,invM2);
		        if (TMath::Abs(invM1-mrho) < 1.25 && pt[0] < 1.2 && pt[1] < 1.2 && pt[2] < 1.2 && pt[3] < 1.2)
		        {
		             h4->Fill(invM2);
			     num_entries += 1;
			     if (TMath::Abs(invM2-0.745) < 3*0.07)
			     {
			         h5->Fill(calc_fourmass(p1,p2,p3,p4));
			     }
			    
		         }
	         
	                d1->Fill(invM3,invM4);
		        if (TMath::Abs(invM3-mrho) < 1.25 && pt[0] < 1.2 && pt[1] < 1.2 && pt[2] < 1.2 && pt[3] < 1.2)
		        {
		             h4->Fill(invM4);
			     num_entries += 1;
			     if (TMath::Abs(invM4-0.745) < 3*0.07)
			     {
			         Double_t p1[3] = {px[0],py[0],pz[0]};
				 Double_t p2[3] = {px[1],py[1],pz[1]};
				 Double_t p3[3] = {px[2],py[2],pz[2]};
				 Double_t p4[3] = {px[3],py[3],pz[3]};
				   
			         h5->Fill(calc_fourmass(p1,p2,p3,p4));
			     }
		        }
		      
	           }
		   
	       
	       }
	       
	   }
	  
		
	 }    


     }
      
   
      
    TCanvas *c1 = new TCanvas("c1","c1",1800,800);
    c1->Divide(3,1);
    d1->SetTitle("Inv. Mass / Inv. Mass ; Inv. Mass [GeV] ; Inv. Mass [GeV]");
    
    gStyle->SetPalette(kCividis);
    gStyle->SetOptStat(false);
    c1->cd(1); h1->Draw();
    c1->cd(2); h2->Draw();
    c1->cd(3); h3->Draw();
    
    const char* name1 = "X-projection";
    const char* name2 = "Y-projection";
    
    auto Xhist = d1->ProjectionX(name1);
    auto Yhist = d1->ProjectionY(name2);
    Yhist->Sumw2();
    Xhist->Sumw2();
     
    Xhist->SetTitle("X-projection; Inv. Mass [GeV] ; ");
    Yhist->SetTitle("Y-projection; Inv. Mass [GeV] ; ");
    
    TCanvas *hist2d = new TCanvas("hist2d","hist2d",800,600);
    hist2d->DrawFrame(0,0,3,3);
    d1->Draw("Colz");
    TLatex Tl;
    Tl.SetTextAlign(12);
    Tl.SetTextSize(0.075);
    Tl.DrawLatex(1.95,2.15,"#font[22]{CMS}");
    TGaxis *A1 = new TGaxis(0,3,0,3,0.5,10,510);
    A1->SetTitle("#sqrt(s) = 13 TeV");
    A1->SetTitleSize(0.03);
    A1->Draw("SAME");
    
    TCanvas *projs = new TCanvas("projs","projs",1200,600);
    projs->Divide(2,1);
    projs->cd(1); Xhist->Draw();
    projs->cd(2); Yhist->Draw();
    
    
    
    TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
    c2->Divide(2,2);
    c2->cd(1); d1->Draw("Colz");
    c2->cd(2); Xhist->Draw("E1");
    c2->cd(3); Yhist->Draw("E1");
    c2->cd(4); h4->Draw("E1");
    
    Xhist->SetLineColor(kRed);
    h1->SetLineColor(kBlue);
    Yhist->SetLineColor(kRed);
    h4->SetLineColor(kBlack);
   
   
    Double_t params[6];
    Double_t errs[6];
    myminimizer(params,errs);
    Double_t x[1400], y[1400];
    
    //background fit
    //params[0] =   6.06643e+04;
    //params[1] = 2.55552e-01  ;
    //params[2] = 9.25497e-01;
    //params[3] =  3.92749e+00 ;
   
    for (Int_t i=0;i< 1400;i++) 
    {
        x[i] = 0.3+i*0.01;
        y[i] = totalfit(x[i],params);
     }

    auto g2 = new TGraph(1400,x,y);
    g2->SetLineWidth(2);
    g2->SetLineColor(kRed);
    g2->Draw("SAME");
    
    Double_t x1[300], y1[300];
    for (Int_t i=0;i< 300;i++) 
    {
        Double_t params1[3] = {params[3],params[4],params[5]};
        x1[i] = 0.35+i*0.001;
        y1[i] = gauss(x1[i],params1);
     }

    auto g3 = new TGraph(300,x1,y1);
    g3->SetLineWidth(2);
    g3->SetLineColor(60);
    g3->Draw("SAME");
    
    
    Double_t x2[500], y2[500];
    for (Int_t i=0;i< 500;i++) 
    {
        x2[i] = 0.5+i*0.001;
	Double_t params2[3] = {params[0],params[1],params[2]};
        y2[i] = gauss(x2[i],params2);
     }

    auto g4 = new TGraph(500,x2,y2);
    g4->SetLineWidth(2);
    g4->SetLineColor(60);
    g4->Draw("SAME");
    
    for (Int_t i=0;i< 300;i++) 
    {
        Double_t back_params[4] = {6.06557e+04,2.55558e-01,9.25427e-01,3.92736e+00};
        y[i] = background(x[i],back_params);
     }

    auto g5 = new TGraph(1400,x,y);
    g5->SetLineWidth(2);
    g5->SetLineColor(kOrange);
    g5->Draw("SAME");
    
    std::cout << "Number of bins: " << "  "  << numbins << std::endl;
    
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    
    h5->Draw();
    

    
    
     

}



