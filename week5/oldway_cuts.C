#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
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

Double_t breitwigner(float x,Double_t *par)
{
    Double_t value = par[0]*TMath::BreitWigner(x,par[1],par[2]);
    return value;
}
 
 
Double_t background(Float_t x,Double_t *par)
{
    Double_t value = par[0]*TMath::Power(x-par[1],par[2])*TMath::Exp(-x*par[3]);
    return value;
}

Double_t totalfit(Float_t x,Double_t *par)
{
   Double_t bp[4] = {7.37157e+04 , 2.47286e-01,1.13211e+00, 4.88666e+00 };
   Double_t pars1[4] = {par[0],par[1],par[2],par[3]};
   Double_t pars2[3] = {par[4],par[5],par[6]};
   Double_t pars3[3] = {par[7],par[8],par[9]};
   
   Double_t value = background(x,pars1)+gauss(x,pars2)+breitwigner(x,pars3);
   return value;
}


Double_t dist(Double_t dxy1, Double_t dz1, Double_t dxy2, Double_t dz2)
{
    return TMath::Sqrt(TMath::Power(dxy1-dxy2,2)+TMath::Power(dz1-dz2,2));
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
       if (yval > 0 && xval > 0.25)
       {
          dx = (totalfit(xval, par) - yval)/err;
          chi2 += (dx*dx);
	  numbins++;
       }
   }
       //std::cout << chi2 << std::endl;     
       f = chi2 ;
}


const Int_t numpar = 10;
Int_t num_entries = 0;

void myminimizer(Double_t *par, Double_t *err)
{
    Double_t arglist[10] = {7.37157e+04 , 2.47286e-01,1.13211e+00, 4.88666e+00,5.70120e+02, 7.40403e-01, 6.65430e-02,2.68407e+01,4.97216e-01, 2.18391e-02};
    
    TMinuit *gMinuit2 = new TMinuit(10);
    gMinuit2->SetFCN(fcn);
    
    Int_t ierflg = 0 ;
    gMinuit2->mnexcm("SET ERR", arglist ,1,ierflg);
    gMinuit2->mnparm(0, "A", arglist[0], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(1, "B", arglist[1], 0.0001,0, 0, ierflg);
    gMinuit2->mnparm(2, "C",arglist[2], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(3, "D",arglist[3], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(4, "N", arglist[4], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(5, "mu", arglist[5], 0.0001,0, 0, ierflg);
    gMinuit2->mnparm(6, "sigma",arglist[6], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(7, "k", arglist[7], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(8, "mu2", arglist[8], 0.0001,0, 0, ierflg);   
    gMinuit2->mnparm(9, "gamma", arglist[9], 0.0001,0, 0, ierflg); 

    gMinuit2->mnexcm("MIGRAD", arglist , 0,ierflg);
    
    for (int j = 0; j < numpar; j++)
    {
       gMinuit2->GetParameter(j,par[j],err[j]);
    }
}


//                     ----------------------------------------------------

void oldway_cuts()
{
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,0.2,2.5,200,0.2,2.5);
    
    TH1F *h1 = new TH1F("h1","px",200,-3,3);   
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-4,4);
    
    TH2F *d3 = new TH2F("d3","InvM1 / InvM2",200,0.2,2.5,200,0.2,2.5);
    TH2F *d4 = new TH2F("d4","invM / pT" , 100,0,3,100,0,3);
    TH2F *d5 = new TH2F("d5","dxy / dxy" , 100,-2,2,100,-2,2);
    TH2F *d6 = new TH2F("d6","tdxy / dxy" , 100,-2,2,100,-2,2);
    
    
    static Float_t trk_pt[220], trk_eta[220], trk_phi[220], trk_dedx[220], trk_p[220], trk_dxy[220], trk_dz[220];
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
    my_tree->SetBranchAddress("trk_dxy",&trk_dxy);
    my_tree->SetBranchAddress("trk_dz",&trk_dz);
    
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
	    Float_t dxy[4];
	    Float_t dz[4];
	    
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
		dxy[i] = trk_dxy[i];
		dz[i] = trk_dz[i];
		
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
	    
	    Float_t invM1, invM2,invM3, invM4, dist1, dist2,dist3,dist4;
	    
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
			       
			       dist1 = TMath::Abs(dxy[0]-dxy[1]);
			       dist2 = TMath::Abs(dxy[2]-dxy[3]);
			       dist3 = TMath::Abs(dxy[0]-dxy[2]);
			       dist4 = TMath::Abs(dxy[1]-dxy[3]);
       

		           }
		
		           else
		           {
	                       invM1 = calc_InvM(p1,p2);
                               invM2 = calc_InvM(p3,p4);
		
		               invM3 = calc_InvM(p1,p4);
	                       invM4 = calc_InvM(p2,p3);
			       
			       dist1 = TMath::Abs(dxy[0]-dxy[1]);
			       dist2 = TMath::Abs(dxy[2]-dxy[3]);
			       dist3 = TMath::Abs(dxy[0]-dxy[3]);
			       dist4 = TMath::Abs(dxy[1]-dxy[2]);
			      
		            }
		
	                }
	     
	                else
	                {
		               invM1 = calc_InvM(p1,p3);
	                       invM2 = calc_InvM(p2,p4);
		
                               invM3 = calc_InvM(p1,p4);
                               invM4 = calc_InvM(p2,p3);
			       
			       dist1 = TMath::Abs(dxy[0]-dxy[2]);
			       dist2 = TMath::Abs(dxy[1]-dxy[3]);
			       dist3 = TMath::Abs(dxy[0]-dxy[3]);
			       dist4 = TMath::Abs(dxy[1]-dxy[2]);

	                 }
		   
		        d4->Fill(pt[0]+pt[1]+pt[2]+pt[3],dist1);
			d4->Fill(pt[0]+pt[1]+pt[2]+pt[3],dist2);
			d4->Fill(pt[0]+pt[1]+pt[2]+pt[3],dist3);
			d4->Fill(pt[0]+pt[1]+pt[2]+pt[3],dist4);
			
		   
		       
			
	                d1->Fill(invM1,invM2);
		        if (TMath::Abs(invM1-0.74) < 0.5 && dist1 < 0.3 & dist2 < 0.3 && pt[0]+pt[1]+pt[2]+pt[3] < 2.2 && pt[0]+pt[1]+pt[2]+pt[3] > 0.8)
		        {
		             h4->Fill(invM2);
			     d3->Fill(invM1,invM2);
			     
			     
			     num_entries += 1;
			     if (TMath::Abs(invM2-0.74) < 0.068)
			     {
			         h5->Fill(calc_fourmass(p1,p2,p3,p4));
			     }
			    
		         }
	         
	                d1->Fill(invM3,invM4);
		        if (TMath::Abs(invM3-0.74) < 0.5 && dist3 < 0.3 & dist4 < 0.3 && pt[0]+pt[1]+pt[2]+pt[3] < 2.2 && pt[0]+pt[1]+pt[2]+pt[3] > 0.8)
		        {
		             h4->Fill(invM4);
			     d3->Fill(invM3,invM4);
			     
			     num_entries += 1;
			     if (TMath::Abs(invM4-0.74) < 0.068)
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
    d1->SetTitle(" ; Inv. Mass [GeV] ; Inv. Mass [GeV]");
    
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
    d1->Draw("Colz");
    TLatex Tl;
    Tl.SetTextSize(0.075);
    Tl.DrawLatex(1.95,2.15,"#font[22]{CMS}");
    
    TLatex s1;
    
    s1.SetTextSize(0.034);
    s1.DrawLatex(2.18,2.51, "#sqrt{s} = 13 #font[22]{TeV}");
    //TGaxis *A1 = new TGaxis(0.2,3,2.5,3,0.5,1,23);
    //A1->SetTitle("#sqrt(s) = 13 TeV");
    //A1->SetTitleSize(0.03);
    //A1->Draw();
    
    hist2d->Update();
    hist2d->Draw();
    
    TCanvas *projs = new TCanvas("projs","projs",1200,600);
    projs->Divide(2,1);
    projs->cd(1); Xhist->Draw();
    projs->cd(2); Yhist->Draw();
    
    
    auto Xhist2 = d3->ProjectionX("X");
    auto Yhist2 = d3->ProjectionY("Y");
    Yhist2->Sumw2();
    Xhist2->Sumw2();
     
    Xhist2->SetTitle("X-projection; Inv. Mass [GeV] ; ");
    Yhist2->SetTitle("Y-projection; Inv. Mass [GeV] ; ");
    
    TCanvas *hist2d2 = new TCanvas("hist2d2","hist2d",800,600);
    d3->Draw("Colz");
    TLatex Tl3;
    Tl3.SetTextAlign(12);
    Tl3.SetTextSize(0.075);
    Tl3.DrawLatex(1.95,2.15,"#font[22]{CMS}");
    //TGaxis *A2 = new TGaxis(0,3,0,3,0.5,10,510);
    //A2->SetTitle("#sqrt(s) = 13 TeV");
    //A2->SetTitleSize(0.03);
    //A2->Draw("SAME");
    
    TCanvas *projs2 = new TCanvas("projs2","projs2",1200,600);
    projs2->Divide(2,1);
    projs2->cd(1); Xhist2->Draw();
    projs2->cd(2); Yhist2->Draw();
    
    
    
    TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
    c2->Divide(2,2);
    c2->cd(1); d1->Draw("Colz");
    c2->cd(2); Xhist->Draw("E1");
    c2->cd(3); Yhist->Draw("E1");
    
    Xhist->SetLineColor(kRed);
    Yhist->SetLineColor(kRed);
     
    TCanvas *fitgraph = new TCanvas("fitgraph","fitgraph",800,600); 
    h4->Draw("E1");
    h4->SetTitle(" ; Inv. Mass [GeV] ; ");
    TLatex *Tl4 = new TLatex();
    Tl4->SetTextSize(0.075);
    Tl4->DrawLatex(1.4*0.78,2300*0.86,"#font[22]{CMS}");
    
     
    TLatex s2;
    s2.SetTextSize(0.034);
    s2.DrawLatex(0.87*1.4,2301, "#sqrt{s} = 13 #font[22]{TeV}");

    h1->SetLineColor(kBlue);
    h4->SetLineColor(39);
   
   
    Double_t params[10];
    Double_t errs[10];
    myminimizer(params,errs);
    Double_t x[1400], y[1400];
    
    //background fit
    //params[0] =   6.06643e+04;
    //params[1] = 2.55552e-01  ;
    //params[2] = 9.25497e-01;
    //params[3] =  3.92749e+00 ;
   
    for (Int_t i=0;i< 1400;i++) 
    {
        x[i] = 0.26+i*0.005;
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
        x1[i] = 0.40+i*0.00065;
        y1[i] = breitwigner(x1[i],params1);
     }

    auto g3 = new TGraph(300,x1,y1);
    g3->SetLineWidth(2);
    g3->SetLineColor(kGreen);
    g3->Draw("SAME");
    
    
    Double_t x2[500], y2[500];
    for (Int_t i=0;i< 500;i++) 
    {
        x2[i] = 0.5+i*0.00091;
	Double_t params2[3] = {params[0],params[1],params[2]};
        y2[i] = gauss(x2[i],params2);
     }

    auto g4 = new TGraph(500,x2,y2);
    g4->SetLineWidth(2);
    g4->SetLineColor(kGreen);
    g4->Draw("SAME");
    
    for (Int_t i=0;i< 300;i++) 
    {
        Double_t back_params[4] = {7.37157e+04 , 2.47286e-01,1.13211e+00, 4.88666e+00 };
        y[i] = background(x[i],back_params);
     }

    auto g5 = new TGraph(1400,x,y);
    g5->SetLineWidth(2);
    g5->SetLineColor(kBlue);
    g5->Draw("SAME");
    
    auto legend = new TLegend(0.64,0.26,0.99,0.78);
    //legend->SetHeader("Fit","C");
    legend->AddEntry(h4,"Data","lep");
    legend->AddEntry(g5,"Background: A(x-B)^{C}exp(-Dx)","l");
    legend->AddEntry(g3,"#splitline{Kaon peak:}{ #splitline{#mu = 0.497 #pm 0.001 GeV}{#gamma = 0.021 #pm 0.001 GeV}}","l");
    legend->AddEntry(g4,"#splitline{Rho peak:}{ #splitline{#mu = 0.741 #pm 0.001 GeV}{#sigma = 0.064 #pm 0.001 GeV}} ","l");
    legend->AddEntry(g2,"Total fit: Chi2 / NDof : 475 / 175","l");
    legend->SetTextSize(0.028);
    legend->Draw();
    
    fitgraph->Update();
    fitgraph->Draw();
    
    std::cout << "Number of bins: " << "  "  << numbins << std::endl;
    
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    
    h5->Draw();
    
    TCanvas *c4 = new TCanvas("c4","c4",800,600);
    TLatex Tl2;
    Tl2.SetTextAlign(12);
    Tl2.SetTextSize(0.075);
    Tl2.DrawLatex(1,1,"#font[22]{CMS}");
    
    d4->Draw("Colz");
    d4->SetTitle("Inv. Mass / Total pT; Inv. Mass [GeV] ; Total pT [GeV]");
    
   
    

    
    
     

}



