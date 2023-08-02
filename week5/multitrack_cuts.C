#include <tuple>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <vector>
#include <utility>
#include <algorithm>

struct VarData {
    int i;
    int j;
    Double_t var;
    bool isFour;
};

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

TH1F *h6 = new TH1F("h6","h6", 200,0.2,1.4);
TH1F *h7 = new TH1F("h7","two-mass", 200,0.24,1.6);


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
   Double_t bp[4] = { 1.23189e+04,2.47533e-01,8.05796e-01,3.77877e+00};
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



bool isElementFarFromOthers(const std::vector<double>& data,const std::vector<double>& comp, double threshold) 
{
    bool isFarFromOthers = false;
    for (size_t i = 0; i < data.size(); ++i) 
    {
        double currentElement = data[i];
        for (size_t j = 0; j < comp.size(); ++j) 
	{
            double difference = TMath::Abs(currentElement - comp[j]);
            if (difference <= threshold) 
            {
                isFarFromOthers = true;
                break;
            }
            
        }

        
    }
    
    if  (isFarFromOthers == true)
    { 
         return true;
    }
  
    else
    {
         return false;
    }
  
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
   for(unsigned int i = 0 ; i < h6->GetNbinsX(); ++i)
  {    
       xval = h7->GetBinCenter(i);
       yval = h7->GetBinContent(i);
       err = h7->GetBinError(i);
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
    Double_t arglist[numpar] = {1.97250e+02,6.97620e-01,1.25601e-01,100,0.5,0.01};
    
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
    
    for (int j = 0; j < numpar; j++)
    {
       gMinuit2->GetParameter(j,par[j],err[j]);
    }
}


void rank_arr(std::vector<int>& index, const std::vector<double>& arr) 
{
    int len = arr.size();
    index.resize(len);
    for (int i = 0; i < len; i++) 
    {
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

void invertArray(std::vector<int>& arr) 
{
    std::vector<int> narr(arr);
    for (int i = 0; i < arr.size(); i++)
    {
        int count = 0;
	for (int j = 0; j < arr.size(); j++)
        {
            if (narr[j] == i)
	    {
	        arr[i] = count;
		break;
	    }
	    
	    else
	    {
	        count += 1;
	    }
        }    
   
    }
}



//                     ----------------------------------------------------

void multitrack_cuts()
{
    TH2F *d1 = new TH2F("d1","two-mass / two-mass",100,0.2,2,100,0.2,2);
    
    TH1F *h1 = new TH1F("h1","px",200,-4,4);   
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-5,5);
    
    TH2F *his4 = new TH2F("his4","Delta dxy",100,0,4,100,0.2,1.6);   
    TH2F *his5 = new TH2F("his5","Delta dz",100,0,4,100,0.2,1.6);
    TH2F *his6 = new TH2F("his6","Dist dxyz",100,0,4,100,0.2,1.6);
    
    
    TH1F *h8 = new TH1F("h8","all two-mass", 200,0.24,1.6);
    
    
    TH2F *d2 = new TH2F("d2","two-mass / two-mass",100,0.2,2,100,0.2,2);
    TH2F *d3 = new TH2F("d3","two-mass / pT" , 80,0.2,1.4,80,-15,15);
    TH2F *d4 = new TH2F("d4","all pairs " , 100,0.2,2,100,0.2,2);
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
	    
	    Double_t netcharge = 0;
	    for(int c : q)
	    {
	       netcharge += c;
	    }
	    
	    if (netcharge != 0)
	    {
		continue;
	    }
	    
	    h1->Fill(px[0]+px[1]+px[2]+px[3]+px[4]+px[5]-6500*(ThxL+ThxR));
	    h2->Fill(py[0]+py[1]+py[2]+py[3]+py[4]+py[5]+6500*(ThyL+ThyR));
	    h3->Fill(pz[0]+pz[1]+pz[2]+pz[3]+pz[4]+pz[5]);
		
	    
	    Double_t invM1, invM2,invM3, invM4;
	    Double_t tmas1, tmas2,tmas3,tmas4;
	    Double_t p1[3] = {px[0],py[0],pz[0]};
	    Double_t p2[3] = {px[1],py[1],pz[1]};
            Double_t p3[3] = {px[2],py[2],pz[2]};
	    Double_t p4[3] = {px[3],py[3],pz[3]};
	    Double_t p5[3] = {px[4],py[4],pz[4]};
	    
	    std::vector <int> rank_dxy;
	    std::vector <int> rank_dz;
	    std::vector <int> rank_dist;
	    
	    std::vector <double> dxys;
	    std::vector <double> dzs;
	    std::vector <double> dds;
	    

	   
	    
	    for (int i = 0; i < track; i++)
	    {
	        dxys.push_back(TMath::Abs(dxy[i]));
		dzs.push_back(TMath::Abs(dz[i]));
		dds.push_back(TMath::Sqrt(dz[i]*dz[i]+dxy[i]*dxy[i]));
	    }
	    
	    rank_arr(rank_dxy, dxys);
	    rank_arr(rank_dz, dzs);
	    rank_arr(rank_dist, dds);
	    
	    
	    
	    std::vector <int> Rmap_dxy = rank_dxy;
	    std::vector <int> Rmap_dz=rank_dz;
	    std::vector <int> Rmap_dist=rank_dist;
	    
	    invertArray(Rmap_dxy);
	    invertArray(Rmap_dz);
	    invertArray(Rmap_dist);
	    
	    std::vector<VarData> modifmassList;
	    
	    //std::cout << Rmap_dxy[0] << " , " << Rmap_dxy[1] << " , " << Rmap_dxy[2] << " , " << Rmap_dxy[3] << " , " << Rmap_dxy[4] << " , " << Rmap_dxy[5] << std::endl;
	    
	    //std::cout << dxy[0] << " , " << dxy[1] << " , " << dxy[2] << " , " << dxy[3] << " , " << dxy[4] << " , " << dxy[5] << std::endl;
	    

	    std::set<int> usedIndices;
	    bool isFour = false;
	    for (int i1 = 0; i1 < track; i1++)
	    {
	        bool is_fill = false;
	        for (int i2 = i1+1; i2 < track; i2++)
		{
		    std::cout << i1 << " , " << i2 << std::endl;
		    std::cout << dz[Rmap_dxy[i1]] << " , " << dz[Rmap_dxy[i2]] << std::endl;
		    if (q[i1] + q[i2] == 0)
		    {
		        Double_t ps[3] = {px[i1],py[i1],pz[i1]};
	                Double_t ps2[3] = {px[i2],py[i2],pz[i2]};
			Double_t all_mass = calc_InvM(ps,ps2);
			h8->Fill(all_mass);
		    }
		    
		    if (q[Rmap_dxy[i1]] + q[Rmap_dxy[i2]] == 0 &&  !is_fill &&
                    usedIndices.find(i1) == usedIndices.end() && usedIndices.find(i2) == usedIndices.end())
		    {   
		        
		        Double_t pm[3] = {px[Rmap_dxy[i1]],py[Rmap_dxy[i1]],pz[Rmap_dxy[i1]]};
	                Double_t pm2[3] = {px[Rmap_dxy[i2]],py[Rmap_dxy[i2]],pz[Rmap_dxy[i2]]};
		        Double_t mass = calc_InvM(pm,pm2);
			
			Double_t total_pt = 0;
		        for (int l = 0; l < track; l++)
		        {
		            total_pt += pt[l];
		        }  
			
			d3->Fill(mass,zPV);

			
			std::vector<double> cur_dxy = {TMath::Abs(dxy[Rmap_dxy[i1]]),TMath::Abs(dxy[Rmap_dxy[i2]])};
			std::vector<double> cur_dz = {TMath::Abs(dz[Rmap_dxy[i1]]),TMath::Abs(dz[Rmap_dxy[i2]])};
			
			//std::cout << "Dxy:  " << cur_dxy[0] << " , " << cur_dxy[1] << std::endl;
			//std::cout << "Dz:  " << cur_dz[0] << " , " << cur_dz[1] << std::endl;
			 
			
			std::vector<double> dxyv_cut;
	                std::vector<double> dzv_cut;
	    
	                for (int i = 0; i < track; i++)
	                {
	                    if (i != rank_dz[i1] && i != rank_dz[i2])
			    {
		                dxyv_cut.push_back(TMath::Abs(dxy[i]));
		                dzv_cut.push_back(TMath::Abs(dz[i]));
			    }
	                }
			
			Double_t threshold = 0.6;
			
			bool isFour1;
	                bool isFour2;
			
			if (TMath::Abs(cur_dxy[0]-cur_dxy[1]) < threshold)
			{
			    isFour1 = false;
			}
			
			else
			{
			    isFour1 = isElementFarFromOthers(cur_dxy,dxys,threshold);
			}
			
			
			if (TMath::Abs(cur_dz[0]-cur_dz[1]) < threshold)
			{
			    isFour2 = false;
			}
			
			else
			{
			    isFour2 = isElementFarFromOthers(cur_dz,dzs,threshold*2);
			}
			
	    
	                if (isFour1 == true || isFour2 == true)
	                {
	                    isFour = true;
	                }
			
			//if ((zPV < -1 || zPV > 1))
			{
			    modifmassList.push_back({i1, i2, mass,isFour});
			}
			
			//else
			{
			//    modifmassList.push_back({i1, i2, 0,0});
			}
			usedIndices.insert(i1);
                        usedIndices.insert(i2);
			is_fill = true;
		    
		    }
		    
		    else
		    {
			modifmassList.push_back({i1, i2, 0,isFour});
		    }
		    
               } 
	    }
	    
	    std::cout << "Next" << std::endl; 
		    
		   
	   
            std::vector<double> o_vars;
	    std::vector<double> m_vars;
	    std::set<int> Indices;
	    for (size_t idx1 = 0; idx1 < modifmassList.size(); idx1++)
            {
		 Double_t var1 = modifmassList[idx1].var;
                 int i1 = modifmassList[idx1].i;
                 int j1 = modifmassList[idx1].j;
		 bool isFour = modifmassList[idx1].isFour;
		 
		 if (var1 != 0 && Indices.find(i1) == Indices.end() && Indices.find(j1) == Indices.end())
		 {  
		     Indices.insert(i1);
                     Indices.insert(j1);
		     o_vars.push_back(var1); 
		     
		     if (isFour == false)
		     {
		        m_vars.push_back(var1);      
	             }


	         }
	     }
	     
	     Double_t total_pt = 0;
             for (int l = 0; l < track; l++)
	     {
	         total_pt += pt[l];
	     } 
	     
	     //if (total_pt < 1.6 || total_pt > 2.6)
	     //{
	     //    continue;
	     //}
	     
	     
	     if (o_vars.size() == 2)
	     {
	         d1->Fill(o_vars[0],o_vars[1]);
	     
	         if (TMath::Abs(o_vars[0]-0.74) < 1)
		 {
	             h6->Fill(o_vars[1]);
		 }
		 
		 if (TMath::Abs(o_vars[1]-0.74) < 1)
		 {
	             h6->Fill(o_vars[0]);
		 }
	     
	     }
	     

	     if (o_vars.size() == 3)
	     {
		 

		 //d3->Fill(o_vars[0],TMath::MaxElement(6,pts));
		 //d3->Fill(o_vars[1],TMath::MaxElement(6,pts));
		 //d3->Fill(o_vars[2],TMath::MaxElement(6,pts));
		 
		 
		 
		 
		 if (TMath::Abs(o_vars[0]-0.74) < 1)
		 {
		     d1->Fill(o_vars[1],o_vars[2]);
		     h6->Fill(o_vars[1]);
		     h6->Fill(o_vars[2]);

		 }
		 
		 if (TMath::Abs(o_vars[1]-0.74) < 1)
		 {
		     d1->Fill(o_vars[0],o_vars[2]);
		     h6->Fill(o_vars[0]);
		     h6->Fill(o_vars[2]);

		 }
		 
		 if (TMath::Abs(o_vars[2]-0.74) < 1)
		 {
		     d1->Fill(o_vars[0],o_vars[1]);
		     
		     h6->Fill(o_vars[0]);
		     h6->Fill(o_vars[1]);

		 }
		 
		 
		 
	     }
	     
	     
	     
	     if (m_vars.size() == 2)
	     {

		 d2->Fill(m_vars[0],m_vars[1]);
	         h7->Fill(m_vars[1]);
	         h7->Fill(m_vars[0]);
		 
	     
	     }
	     

	     if (m_vars.size() == 3)
	     {
		 
		 
		 //h7->Fill(m_vars[0]);
		 //h7->Fill(m_vars[1]);
		 //h7->Fill(m_vars[2]);
		 
		 
		 
		 if (TMath::Abs(m_vars[0]-0.740) < 1)
		 {
		     d2->Fill(m_vars[1],m_vars[2]);
		     h7->Fill(m_vars[1]);
		     h7->Fill(m_vars[2]);

		 }
		 
		 if (TMath::Abs(m_vars[1]-0.740) < 1)
		 {
		     d2->Fill(m_vars[0],m_vars[2]);
		     
		     h7->Fill(m_vars[0]);
		     h7->Fill(m_vars[2]);
		 }
		 
		 if (TMath::Abs(m_vars[2]-0.740) < 1)
		 {
		     d2->Fill(m_vars[0],m_vars[1]);
		     
		     h7->Fill(m_vars[0]);
		     h7->Fill(m_vars[1]);

		 }
		 
		 
		 
	     }
	     
	     
	    
	     
	     
	     
 
		
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
    
    TCanvas *c4 = new TCanvas("c4","c4",1600,600);
    c4->Divide(2,1);
    c4->cd(1); d1->Draw("Colz");
    c4->cd(2); d2->Draw("Colz");
    
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
    
    TCanvas *c3 = new TCanvas("c3","c3",1600,600);
    c3->Divide(2,1);
    c3->cd(1); Xhist->Draw("E1");
    Xhist2->Draw("SAMEE1");
    
    auto legend = new TLegend(0.64,0.4,0.99,0.78);
    legend->AddEntry(Xhist,"With outliers","lep");
    legend->AddEntry(Xhist2,"Without outliers","lep");
    legend->SetTextSize(0.028);
    legend->Draw();
    
    
    c3->cd(2); Yhist->Draw("E1");
    Yhist2->Draw("SAMEE1");
    
    
    Xhist2->SetLineColor(kRed);
    Yhist2->SetLineColor(kRed);
   
    TCanvas *c5 = new TCanvas("c5","c5",800,600);
    h6->SetTitle("Outliers; Inv. Mass [GeV] ; ");
    h6->Draw("E1");
    h7->Draw("SAMEE1");
    
    h6->SetLineColor(kBlue);
    h7->SetLineColor(kRed);
    
    auto legend2 = new TLegend(0.64,0.4,0.99,0.78);
    legend2->AddEntry(Yhist,"With outliers","lep");
    legend2->AddEntry(Yhist2,"Without outliers","lep");
    legend2->SetTextSize(0.028);
    legend2->Draw();
    
    
    
    
    TCanvas *c6 = new TCanvas("c6","c6",800,600);
    h7->Draw("SAMEE1");
    
    
    //h6->SetLineColor(kRed);
   
    
    
    h8->SetTitle("Six-track: Pairings; Invariant Mass1 [GeV] ; ");
    
    Double_t params[10];
    Double_t errs[10];
    myminimizer(params,errs);
    Double_t x[1400], y[1400];
    
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
        y1[i] = gauss(x1[i],params1);
     }

    auto g3 = new TGraph(300,x1,y1);
    g3->SetLineWidth(2);
    g3->SetLineColor(kGreen);
    g3->Draw("SAME");
    
    
    Double_t x2[500], y2[500];
    for (Int_t i=0;i< 500;i++) 
    {
        x2[i] = 0.4+i*0.0015;
	Double_t params2[3] = {params[0],params[1],params[2]};
        y2[i] = gauss(x2[i],params2);
     }

    auto g4 = new TGraph(500,x2,y2);
    g4->SetLineWidth(2);
    g4->SetLineColor(kGreen);
    g4->Draw("SAME");
    
    for (Int_t i=0;i< 300;i++) 
    {
        Double_t back_params[4] = {1.23189e+04,2.47533e-01,8.05796e-01,3.77877e+00};
        y[i] = background(x[i],back_params);
     }

    auto g5 = new TGraph(1400,x,y);
    g5->SetLineWidth(2);
    g5->SetLineColor(kBlue);
    g5->Draw("SAME");
    
     auto legend6 = new TLegend(0.64,0.26,0.99,0.78);
    //legend->SetHeader("Fit","C");
    legend6->AddEntry(h7,"Data","lep");
    legend6->AddEntry(g5,"Background: A(x-B)^{C}exp(-Dx)","l");
    legend6->AddEntry(g3,"#splitline{Kaon peak:}{#mu = 0.50 #pm 0.02 GeV}","l");
    legend6->AddEntry(g4,"#splitline{Rho peak:}{ #mu = 0.70 #pm 0.05 GeV} ","l");
    legend6->AddEntry(g2,"Total fit: Chi2 / NDof : 377 / 179","l");
    legend6->SetTextSize(0.028);
    legend6->Draw();
   
    
    
    
     h7->SetTitle("; Inv Mass [GeV] ; ");
    
     std::cout << "Number of bins: " << "  "  << numbins << std::endl;
    
    
}



