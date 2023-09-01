#include <tuple>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <vector>
#include <utility>

struct VarData {
    int i;
    int j;
    Double_t var;
    Double_t distance;
};

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

TH1F *h4 = new TH1F("h4","h4", 200,0.2,1.4);
TH1F *h5 = new TH1F("h5","h5", 100,0.4,3);


     //  Handy functions     ----------------------------------------------------
Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

Float_t calc_E(Double_t *p1)
{
    return TMath::Sqrt(mpi*mpi+p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
}

Double_t calc_InvM(Double_t *p1, Double_t *p2)
{
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    
    return TMath::Sqrt(TMath::Power(E1+E2,2) - TMath::Power(p1[0]+p2[0],2)- TMath::Power(p1[1]+p2[1],2)- TMath::Power(p1[2]+p2[2],2));
}

Double_t calc_threemass(Double_t *p1, Double_t *p2,Double_t *p3)
{
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    Double_t E3 = calc_E(p3);
    
    return TMath::Sqrt(TMath::Power(E1+E2+E3,2) - TMath::Power(p1[0]+p2[0]+p3[0],2)- TMath::Power(p1[1]+p2[1]+p3[1],2)- TMath::Power(p1[2]+p2[2]+p3[2],2));
}



Double_t calc_fourmass(Double_t *p1, Double_t *p2,Double_t *p3,Double_t *p4)
{
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    Double_t E3 = calc_E(p3);
    Double_t E4 = calc_E(p4);
    
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

Double_t dist(Double_t dxy1, Double_t dz1, Double_t dxy2, Double_t dz2)
{
    return TMath::Sqrt(TMath::Power(dxy1-dxy2,2)+TMath::Power(dz1-dz2,2));
}

Double_t angle(Double_t *p1, Double_t *p2)
{
    return (((p1[0]*p2[0]+p1[1]*p2[1])/(TMath::Sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])*TMath::Sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]))));
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

void five_track()
{
    TH2F *d1 = new TH2F("d1","two-mass / two-mass",100,0.2,2,100,0.2,2);
    
    TH1F *h1 = new TH1F("h1","px",200,-4,4);   
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-5,5);
    
    TH2F *his4 = new TH2F("his4","Delta dxy",100,0,4,100,0.2,1.6);   
    TH2F *his5 = new TH2F("his5","Delta dz",100,0,4,100,0.2,1.6);
    TH2F *his6 = new TH2F("his6","Dist dxyz",100,0,4,100,0.2,1.6);
    
    TH1F *h6 = new TH1F("h6","two-mass", 100,0.2,1.4);
    TH1F *h6ex = new TH1F("h6ex","two-mass", 100,0.2,1.4);
    
    TH1F *h7 = new TH1F("h7","Angle", 100,-1,1);
    TH2F *d2 = new TH2F("d2","two-mass / three-mass",100,0,3,100,0,4);
    TH2F *d3 = new TH2F("d3","three-mass / pt" , 100,0,3,100,0,4);
    TH2F *d4 = new TH2F("d4","three-mass / pt" , 100,-2,2,100,-2,2);
    TH2F *d5 = new TH2F("d5","dxy / dxy" , 100,0,4,100,0,4);
    TH2F *d6 = new TH2F("d6","dxy / dxy" , 100,0,4,100,0,4);
    TH2F *d7 = new TH2F("d7","Angle / two-mass" , 100,0,3,100,-1,1);
    TH2F *d8 = new TH2F("d8","Dist (dxy+dz) / two-mass" , 100,0,4,100,0,3);
    TH2F *d9 = new TH2F("d9","E1 / E2" , 100,-1,1,100,0,3);
    
    
    static Float_t trk_pt[420], trk_eta[420], trk_phi[420], trk_dedx[420], trk_p[420], trk_dxy[420], trk_dz[420];
    static Int_t trk_isK[420], trk_isPi[420], trk_isP[420];
    static Float_t ThxR, ThxL, ThyR, ThyL, zPV;
    static Int_t trk_q[420];
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
    my_tree->SetBranchAddress("zPV",&zPV);
    
    my_tree->SetBranchAddress("ThxR",&ThxR);
    my_tree->SetBranchAddress("ThxL",&ThxL);
    my_tree->SetBranchAddress("ThyR",&ThyR);
    my_tree->SetBranchAddress("ThyL",&ThyL);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
 
    int numEnt = my_tree->GetEntries();
    Int_t track = 5;

    for (int irow = 0; irow< numEnt; irow++)
    {    
        my_tree->GetEntry(irow);
        if (ntrk == track)
	{
	    Float_t px[track];
            Float_t py[track];
	    Float_t pz[track];
	    Float_t pt[track];
	    Float_t q[track];
	    Float_t dedx[track];
	    Float_t eta[track];
	    Float_t phi[track];
	    Float_t p[track];
	    Float_t dxy[track];
	    Float_t dz[track];
	    
	    Int_t isPi[track];
	    Int_t isP[track];
	    Int_t isK[track];
	
	    for (int i = 0; i < track; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		pt[i] = trk_pt[i];
		p[i] = trk_p[i];
		q[i] = trk_q[i];
		dedx[i] = trk_dedx[i];
		dxy[i] = trk_dxy[i];
		dz[i] = trk_dz[i];
		
		isPi[i] = trk_isPi[i];
	        isP[i] = trk_isP[i];
		isK[i] = trk_isK[i];
	    }
	    
	    h1->Fill(px[0]+px[1]+px[2]+px[3]+px[4]-6500*(ThxL+ThxR));
	    h2->Fill(py[0]+py[1]+py[2]+py[3]+py[4]+6500*(ThyL+ThyR));
	    h3->Fill(pz[0]+pz[1]+pz[2]+pz[3]+pz[4]);
		
	    
	    Double_t invM1, invM2,invM3, invM4;
	    Double_t tmas1, tmas2,tmas3,tmas4;
	    Double_t p1[3] = {px[0],py[0],pz[0]};
	    Double_t p2[3] = {px[1],py[1],pz[1]};
            Double_t p3[3] = {px[2],py[2],pz[2]};
	    Double_t p4[3] = {px[3],py[3],pz[3]};
	    Double_t p5[3] = {px[4],py[4],pz[4]};
	    
	    if (TMath::Abs(q[0]+q[1]+q[2]+q[3]+q[4]) == 3)
            {
	    	 continue;
            }
		    
	    
	    std::vector<VarData> massList;
	    Int_t count = 0; 
	    for (Int_t i = 0; i < track; i++)
	    {    
	        
	        for (Int_t j = i+1; j < track;j++)
		{
		    Double_t ps[3] = {px[i],py[i],pz[i]};
	            Double_t ps2[3] = {px[j],py[j],pz[j]};
		    
		    if (q[i]+q[j] == 0)
		    {
		         Double_t var = calc_InvM(ps,ps2);
			 Double_t distance = dist(dxy[i],dz[i],dxy[j],dz[j]);
			 
		         his4->Fill(TMath::Abs(dxy[i]-dxy[j]),calc_InvM(ps,ps2));
			 his5->Fill(TMath::Abs(dz[i]-dz[j]),calc_InvM(ps,ps2));
			 his6->Fill(dist(dxy[i],dz[i],dxy[j],dz[j]),calc_InvM(ps,ps2));
			 
			 d9->Fill(angle(ps,ps2),dist(dxy[i],dz[i],dxy[j],dz[j]));
			 
			 d8->Fill(var,dist(dxy[i],dz[i],dxy[j],dz[j]));
			 h6->Fill(var);
			 d7->Fill(var,angle(ps,ps2));
			 h7->Fill(angle(ps,ps2));
			 
			 if (dist(dxy[i],dz[i],dxy[j],dz[j]) > 0.4)
			 {
			     continue;
			 }
			 
			 Double_t total_pt = 0;
			 for (int l = 0; l < track; l++)
			 {
			    total_pt += pt[l];
			 }
			
			 
			 if (angle(ps,ps2) < -0.1 && total_pt < 2)
			 {   
			     massList.push_back({i, j, var, distance});
			     h6ex->Fill(var);  
			     count += 1;
			     
			 }
			 
			 else
			 {
			     massList.push_back({i, j, 0, distance});
			 }
			 
			 
		    }
		 }
	    
	    }
	    
	    
	    for (size_t idx1 = 0; idx1 < massList.size(); idx1++)
            {
                for (size_t idx2 = idx1 + 1; idx2 < massList.size(); idx2++)
                {
                    Double_t var1 = massList[idx1].var;
                    Double_t var2 = massList[idx2].var;
                    int i1 = massList[idx1].i;
                    int j1 = massList[idx1].j;
	            int i2 = massList[idx2].i;
                    int j2 = massList[idx2].j;
		    Double_t dist1 = massList[idx1].distance;
		    Double_t dist2 = massList[idx1].distance;
		   
		   
		//   const Double_t threshold = 0.3; 

                   if (i1 != i2 && j1 != j2 && var1 != 0 && var2 != 0)
                   {

                      if (dist1 < 0.4 && dist2 < 0.4)
                      {
		      
		      
		      
	              }     
                //        for (Int_t t = 0 ; t < track; t++)
                //         {
		//	      std::cout << dist1 << " " << dist2 << std::endl;
                //              if (t != i1 && t!= j1 && t != i2 && t != j2)
                //              {
                //                  Double_t oDist1 = dist(dxy[t], dz[t], dxy[i1], dz[i1]);
		//		  Double_t oDist2 = dist(dxy[t], dz[t], dxy[j1], dz[j1]);
				 
		//		  Double_t oDist3 = dist(dxy[t], dz[t], dxy[i2], dz[i2]);
		//		  Double_t oDist4 = dist(dxy[t], dz[t], dxy[j2], dz[j2]);
		//		  
                //                  if ((oDist1+oDist2)/2 > threshold  && (oDist3+oDist4)/2 > threshold)
                //                  {
		//		      d1->Fill(var1, var2);
                //                  }
                //               }
                //          }


		 
                //      }	    
                //    }
                // }
	     //}
	       
		 
		
	 }    


     }
      
   
      
    TCanvas *c1 = new TCanvas("c1","c1",1800,800);
    c1->Divide(3,1);
    d1->SetTitle("Two-mass / Two-mass ; Invariant Mass1 [GeV] ; Invariant Mass2 [GeV]");
    
    gStyle->SetPalette(kCividis);
    gStyle->SetOptStat(false);
    c1->cd(1); h1->Draw();
    c1->cd(2); h2->Draw();
    c1->cd(3); h3->Draw();
    
    TCanvas *c10 = new TCanvas("c10","c10",1800,800);
    c10->Divide(3,1);

    c10->cd(1); his4->Draw("Colz");
    c10->cd(2); his5->Draw("Colz");
    c10->cd(3); his6->Draw("Colz");
    
    const char* name1 = "X-projection";
    const char* name2 = "Y-projection";
    
    auto Xhist = d1->ProjectionX(name1);
    auto Yhist = d1->ProjectionY(name2);
    Yhist->Sumw2();
    Xhist->Sumw2();
     
    Xhist->SetTitle("X-projection; Energy [GeV] ; #Counts");
    Yhist->SetTitle("Y-projection; Energy [GeV] ; #Counts");
    
    TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
    c2->Divide(2,2);
    c2->cd(1); d1->Draw("Colz");
    c2->cd(2); Xhist->Draw("E1");
    c2->cd(3); Yhist->Draw("E1");
    
    Xhist->SetLineColor(kRed);
    Yhist->SetLineColor(kRed);
   
    std::cout << "Number of bins: " << "  "  << numbins << std::endl;    
    
    TCanvas *c5 = new TCanvas("c5","c5",800,600);
    h6->Draw("E1");
    h6ex->Draw("SAMEE1");
    h6ex->SetLineColor(kRed);
    
    
    TCanvas *c6 = new TCanvas("c6","c6",800,600);
    d7->Draw("Colz");
    
    TCanvas *c7 = new TCanvas("c7","c7",800,600);
    d8->Draw("Colz");
    
    TCanvas *c8 = new TCanvas("c8","c8",800,600);
    d9->Draw("Colz");
    
    TCanvas *c9 = new TCanvas("c9","c9",800,600);
    h7->Draw("E1");

}



