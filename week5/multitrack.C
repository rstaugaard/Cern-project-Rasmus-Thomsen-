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

double calculateMean(const std::vector<int>& values) {
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / values.size();
}

double calc_correlation(const std::vector<int>& x, const std::vector<int>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::invalid_argument("Both vectors must have the same non-zero size.");
    }

    double meanX = calculateMean(x);
    double meanY = calculateMean(y);

    double numerator = 0.0;
    double denomX = 0.0;
    double denomY = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        numerator += (x[i] - meanX) * (y[i] - meanY);
        denomX += std::pow(x[i] - meanX, 2);
        denomY += std::pow(y[i] - meanY, 2);
    }

    double denominator = std::sqrt(denomX) * std::sqrt(denomY);

    if (denominator == 0.0) {
        return 0.0; // If denominator is 0, correlation is not defined, return 0.
    }

    return numerator / denominator;
}



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


void rank_arr(Int_t *index, Double_t *arr, Int_t len) 
{
    for (int i = 0; i < len; i++) 
    {
        int indx = -1;
        int rank = 0;
        for (int j = 0; j < len; j++) 
	{
            if (arr[j] < arr[i]) 
	    {
                rank++;
            }
        }

        index[i] = rank;
    }
}



//                     ----------------------------------------------------

void multitrack()
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
    TH2F *d2 = new TH2F("d2","dz used to pair",100,0.2,2,100,0.2,2);
    TH2F *d3 = new TH2F("d3","three-mass / pt" , 5,0,4,5,0,4);
    TH2F *d4 = new TH2F("d4","dxy used to pair " , 100,0.2,2,100,0.2,2);
    TH2F *d5 = new TH2F("d5","dxy / dxy" , 100,0,4,100,0,4);
    TH2F *d6 = new TH2F("d6","dxy / dxy" , 100,0,4,100,0,4);
    TH2F *d7 = new TH2F("d7","Angle / two-mass" , 100,0,3,100,-1,1);
    TH2F *d8 = new TH2F("d8","Dist (dxy+dz) / two-mass" , 100,0,4,100,0,3);
    TH2F *d9 = new TH2F("d9","E1 / E2" , 100,-1,1,100,0,3);
    
    
    static Float_t trk_pt[420], trk_eta[420], trk_phi[420], trk_dedx[420], trk_p[420], trk_dxy[420], trk_dz[420];
    static Int_t trk_isK[420], trk_isPi[420], trk_isP[420];
    
    std::vector<int> dxyR; std::vector<int> dzR;
    
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
    Int_t track = 6;

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
	    
	    Int_t rank_dxy[track];
	    Int_t rank_dz[track];
	    Double_t dxys[track];
	    Double_t dzs[track];
	    
	    for (int i = 0; i < track; i++)
	    {
	        dxys[i] =  TMath::Abs(dxy[i]);
		dzs[i] =  TMath::Abs(dz[i]);
	    }
	    
	    rank_arr(rank_dxy, dxys, track);
	    rank_arr(rank_dz, dzs, track);
	    
	    for (int k = 0; k < track; k++)
	    {
	        dxyR.push_back(rank_dxy[k]);
		dzR.push_back(rank_dz[k]);
		
		//std::cout << "---" << rank_dxy[0]<< " , " << rank_dxy[1] << " , "<< rank_dxy[2] << " , " << rank_dxy[3] << std::endl;
		//std::cout << dxy[0]<< " , " << dxy[1] << " , "<< dxy[2] << " , " << dxy[3] << std::endl;
		d3->Fill(rank_dxy[k],rank_dz[k]);
	    }
	    
	    std::vector<VarData> modifmassList;
	    std::vector<VarData> modifmassList2;
	    std::vector<VarData> massList;
	    
	    for (Int_t i = 0; i < track; i++)
	    {    
	        
	        for (Int_t j = i+1; j < track;j++)
		{
		    Double_t ps[3] = {px[i],py[i],pz[i]};
	            Double_t ps2[3] = {px[j],py[j],pz[j]};
	
		    if (q[i]+q[j] == 0)
		    {
		         Double_t var = calc_InvM(ps,ps2);

			 Double_t total_pt = 0;
			 for (int l = 0; l < track; l++)
			 {
			    total_pt += pt[l];
			 }
			
			 
			 if (zPV < -1 || zPV > 1)
			 {   
			     massList.push_back({i, j, var});
			     
			 }
			 
			 else
			 {
			     massList.push_back({i, j, 0});
			 }
			 
			 
		    }
		    
		    else
	            {
		        massList.push_back({i, j, 0});
	            }
		    
		 }
	    
	    }
	    
	 
	    std::set<int> usedIndices;
	    for (int i1 = 0; i1 < track; i1++)
	    {
	        bool is_fill = false;
	        for (int i2 = i1+1; i2 < track; i2++)
		{
		    if (q[rank_dz[i1]] + q[rank_dz[i2]] == 0 &&  !is_fill &&
                usedIndices.find(i1) == usedIndices.end() && usedIndices.find(i2) == usedIndices.end())
		    {   
		        
		        Double_t pm[3] = {px[rank_dz[i1]],py[rank_dz[i1]],pz[rank_dz[i1]]};
	                Double_t pm2[3] = {px[rank_dz[i2]],py[rank_dz[i2]],pz[rank_dz[i2]]};
		        Double_t mass = calc_InvM(pm,pm2);
			
			if ((zPV < -1 || zPV > 1))
			{
			    modifmassList.push_back({i1, i2, mass});
			}
			
			else
			{
			    modifmassList.push_back({i1, i2, 0});
			}
			usedIndices.insert(i1);
                        usedIndices.insert(i2);
			is_fill = true;
		    
		    }
		    
		    else
		    {
		        modifmassList.push_back({i1, i2, 0});
		    }
		    
               } 
	    } 
		    
		    
		    
            std::set<int> usedIndices2;
	    for (int i1 = 0; i1 < track; i1++)
	    {
	        bool is_fill2 = false;
	        for (int i2 = i1+1; i2 < track; i2++)
		{
		    if (q[rank_dxy[i1]] + q[rank_dxy[i2]] == 0 &&  !is_fill2 &&
                usedIndices2.find(i1) == usedIndices2.end() && usedIndices2.find(i2) == usedIndices2.end())
		    {   
		        
		        Double_t pm[3] = {px[rank_dxy[i1]],py[rank_dxy[i1]],pz[rank_dxy[i1]]};
	                Double_t pm2[3] = {px[rank_dxy[i2]],py[rank_dxy[i2]],pz[rank_dxy[i2]]};
		        Double_t mass = calc_InvM(pm,pm2);
			
			if ((zPV < -1 || zPV > 1))
			{
			    modifmassList2.push_back({i1, i2, mass});
			}
			
			else
			{
			    modifmassList2.push_back({i1, i2, 0});
			}
			usedIndices2.insert(i1);
                        usedIndices2.insert(i2);
			is_fill2 = true;
		    
		    }
		    
		    else
		    {
		        modifmassList2.push_back({i1, i2, 0});
		    }
		   
		    
		    
		    
		
		}
	    
	    } 
	   
 
	    for (size_t idx1 = 0; idx1 < massList.size(); idx1++)
            {
                for (size_t idx2 = idx1 + 1; idx2 <  massList.size(); idx2++)
                {
                   Double_t ovar1 = massList[idx1].var;
                   Double_t ovar2 = massList[idx2].var;
		   int oi1 = massList[idx1].i;
                   int oj1 = massList[idx1].j;
	           int oi2 = massList[idx2].i;
                   int oj2 = massList[idx2].j;
		   
		   Double_t var1 = modifmassList[idx1].var;
                   Double_t var2 = modifmassList[idx2].var;
                   int i1 = modifmassList[idx1].i;
                   int j1 = modifmassList[idx1].j;
	           int i2 = modifmassList[idx2].i;
                   int j2 = modifmassList[idx2].j;
		   
		   Double_t nvar1 = modifmassList2[idx1].var;
                   Double_t nvar2 = modifmassList2[idx2].var;
                   int ni1 = modifmassList2[idx1].i;
                   int nj1 = modifmassList2[idx1].j;
	           int ni2 = modifmassList2[idx2].i;
                   int nj2 = modifmassList2[idx2].j;
		   
                   if (i1 != i2 && i1 != j2 && j1 != j2 && j1 != i2 && var1 != 0 && var2 != 0)
                   {
		      d2->Fill(var1, var2);
                   }
		   
		   if (ni1 != ni2 && ni1 != nj2 && nj1 != nj2 && nj1 != ni2 && nvar1 != 0 && nvar2 != 0)
                   {
		      d4->Fill(nvar1, nvar2);
                   }
		   
		   
		   if (oi1 != oi2 && oj1 != oj2 && oi1 != oj2 && oi2 != oj1 && ovar1 != 0 && ovar2 != 0)
                   {
		      d1->Fill(ovar1, ovar2);
                   }
		   
		   
                 }
		 
		 
	     }
	     
	     //std::cout << "Next" << std::endl;
 
		
	 }    


     }
      
   
    double correlation = calc_correlation(dxyR, dzR);
    
    std::cout << "Spearman correlation coefficient: " << correlation << std::endl;

      
    TCanvas *c1 = new TCanvas("c1","c1",1800,800);
    c1->Divide(3,1);
    d1->SetTitle("Two-mass / Two-mass ; Invariant Mass1 [GeV] ; Invariant Mass2 [GeV]");
    
    gStyle->SetPalette(kCividis);
    gStyle->SetOptStat(false);
    c1->cd(1); h1->Draw();
    c1->cd(2); h2->Draw();
    c1->cd(3); h3->Draw();
    
    TCanvas *c4 = new TCanvas("c4","c4",1600,1200);
    c4->Divide(2,2);
    c4->cd(1); d1->Draw("Colz");
    c4->cd(2); d2->Draw("Colz");
    c4->cd(3); d4->Draw("Colz");
    
    const char* name1 = "X-projection";
    const char* name2 = "Y-projection";
    
    auto Xhist = d1->ProjectionX(name1);
    auto Yhist = d1->ProjectionY(name2);
    Yhist->Sumw2();
    Xhist->Sumw2();
     
    Xhist->SetTitle("X-projection; Energy [GeV] ; ");
    Yhist->SetTitle("Y-projection; Energy [GeV] ; ");
    
    auto Xhist2 = d2->ProjectionX("X");
    auto Yhist2 = d2->ProjectionY("Y");
    Yhist2->Sumw2();
    Xhist2->Sumw2();
     
    Xhist2->SetTitle("X-projection; Energy [GeV] ; ");
    Yhist2->SetTitle("Y-projection; Energy [GeV] ; ");
    
    auto Xhist3 = d4->ProjectionX("X2");
    auto Yhist3 = d4->ProjectionY("Y2");
    Yhist3->Sumw2();
    Xhist3->Sumw2();
     
    Xhist3->SetTitle("X-projection; Energy [GeV] ; ");
    Yhist3->SetTitle("Y-projection; Energy [GeV] ; ");
    
    TCanvas *c3 = new TCanvas("c3","c3",1600,600);
    c3->Divide(2,1);
    c3->cd(1); Xhist->Draw("E1");
    Xhist2->Draw("SAMEE1");
    Xhist3->Draw("SAMEE1");
    
    
    auto legend = new TLegend(0.64,0.4,0.99,0.78);
    legend->AddEntry(Xhist,"All combinations included","lep");
    legend->AddEntry(Xhist2,"Pairs combined using dz","lep");
    legend->AddEntry(Xhist3,"Pairs combined using dxy","lep");
    legend->SetTextSize(0.028);
    legend->Draw();
    
    
    c3->cd(2); Yhist->Draw("E1");
    Yhist2->Draw("SAMEE1");
    Yhist3->Draw("SAMEE1");
    
    auto legend2 = new TLegend(0.64,0.4,0.99,0.78);
    legend2->AddEntry(Yhist,"All combinations included","lep");
    legend2->AddEntry(Yhist2,"Pairs combined using dz","lep");
    legend2->AddEntry(Yhist3,"Pairs combined using dxy","lep");
    legend2->SetTextSize(0.028);
    legend2->Draw();
    
    Xhist2->SetLineColor(kRed);
    Yhist2->SetLineColor(kRed);
    
    Xhist3->SetLineColor(kGreen);
    Yhist3->SetLineColor(kGreen);
    
   
    TCanvas *c5 = new TCanvas("c5","c5",800,600);
    d3->Draw("Colz");
}



