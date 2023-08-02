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

double calculateMean(const std::vector<double>& values) 
{
    double sum = 0.0;
    for (double value : values) 
    {
        sum += value;
    }
    return sum / values.size();
}

double calculateStandardDeviation(const std::vector<double>& values, double mean) 
{    
    double sumSquaredDifferences = 0.0;
    for (double value : values) 
    {
        double difference = value - mean;
        sumSquaredDifferences += difference * difference;
    }
    double variance = sumSquaredDifferences / values.size();
    return std::sqrt(variance);
}

double calculateMedian(std::vector<double>& values) 
{
    size_t n = values.size();
    std::sort(values.begin(), values.end());
    if (n % 2 == 0) 
    {
        return (values[n / 2 - 1] + values[n / 2]) / 2.0;
    } else {
        return values[n / 2];
    }
}

double calculateMAD(std::vector<double>& values, double median)
 {
    std::vector<double> absoluteDeviations;
    for (double value : values) 
    {
        absoluteDeviations.push_back(std::abs(value - median));
    }
    return calculateMedian(absoluteDeviations);
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


const Int_t numpar = 4;
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

void multitrack_cuts()
{
    TH2F *d1 = new TH2F("d1","two-mass / two-mass",100,0.2,2,100,0.2,2);
    
    TH1F *h1 = new TH1F("h1","px",200,-4,4);   
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-5,5);
    
    TH2F *his4 = new TH2F("his4","Delta dxy",100,0,4,100,0.2,1.6);   
    TH2F *his5 = new TH2F("his5","Delta dz",100,0,4,100,0.2,1.6);
    TH2F *his6 = new TH2F("his6","Dist dxyz",100,0,4,100,0.2,1.6);
    
    TH1F *h6 = new TH1F("h6","two-mass", 200,0.24,1.4);
    TH1F *h7 = new TH1F("h7","two-mass", 200,0.24,1.4);
    
    TH1F *h8 = new TH1F("h8","all two-mass", 200,0.24,1.4);
    
    
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
	    
	    Int_t rank_dxy[track];
	    Int_t rank_dz[track];
	    Int_t rank_dist[track];
	    Double_t dxys[track];
	    Double_t dzs[track];
	    Double_t dds[track];
	    
	    
	    std::vector<double> dxyv;
	    std::vector<double> dzv;
	    std::vector<double> ddv;
	    
	    for (int i = 0; i < track; i++)
	    {
	        dxys[i] =  TMath::Abs(dxy[i]);
		dzs[i] =  TMath::Abs(dz[i]);
		dds[i] = TMath::Sqrt(dz[i]*dz[i]+dxy[i]*dxy[i]);
		dxyv.push_back(TMath::Abs(dxy[i]));
		dzv.push_back(TMath::Abs(dz[i]));
		ddv.push_back(TMath::Sqrt(dz[i]*dz[i]+dxy[i]*dxy[i]));
			    
	    }

	    rank_arr(rank_dxy, dxys, track);
	    rank_arr(rank_dz, dzs, track);
	    rank_arr(rank_dist, dds, track);

	    
	    std::vector<VarData> modifmassList;

	    std::set<int> usedIndices;
	    bool isFour = false;
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
			
			Double_t total_pt = 0;
		        for (int l = 0; l < track; l++)
		        {
		            total_pt += pt[l];
		        }  
			
			d3->Fill(mass,zPV);
			
		 
		        //if (total_pt > 1.6 && total_pt < 2.7)
		        {
			   h8->Fill(mass);
			}
			
			std::vector<double> cur_dxy = {TMath::Abs(dxy[rank_dz[i1]]),TMath::Abs(dxy[rank_dz[i2]])};
			std::vector<double> cur_dz = {TMath::Abs(dz[rank_dz[i1]]),TMath::Abs(dz[rank_dz[i2]])};
			
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
			    isFour1 = isElementFarFromOthers(cur_dxy,dxyv,threshold);
			}
			
			
			if (TMath::Abs(cur_dz[0]-cur_dz[1]) < threshold*2)
			{
			    isFour2 = false;
			}
			
			else
			{
			    isFour2 = isElementFarFromOthers(cur_dz,dzv,threshold*2);
			}
			
			
			
			std::cout << isFour2 << "| " << cur_dz[0] << " , " << cur_dz[1] << "| " << dzv_cut[0] <<" , " << dzv_cut[1] << " , " <<dzv_cut[2] <<" , " << dzv_cut[3] << std::endl; 
	    
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
	     
	     if (total_pt < 1.6 || total_pt > 2.6)
	     {
	         continue;
	     }
	     
	     
	     if (o_vars.size() == 2)
	     {
	         d1->Fill(o_vars[0],o_vars[1]);
	     
	         if (TMath::Abs(o_vars[0]-0.74) < 0.4)
		 {
	             h6->Fill(o_vars[1]);
		 }
		 
		 if (TMath::Abs(o_vars[1]-0.74) < 0.4)
		 {
	             h6->Fill(o_vars[0]);
		 }
	     
	     }
	     

	     if (o_vars.size() == 3)
	     {
		 
		 
		 Double_t total_pt = 0;
		 for (int l = 0; l < track; l++)
		 {
		     total_pt += pt[l];
		 } 
		 
		 //h6->Fill(o_vars[0]);
		 //h6->Fill(o_vars[1]);
		 //h6->Fill(o_vars[2]);
		 
		 
		 
		 
		 
		 //d3->Fill(o_vars[0],TMath::MaxElement(6,pts));
		 //d3->Fill(o_vars[1],TMath::MaxElement(6,pts));
		 //d3->Fill(o_vars[2],TMath::MaxElement(6,pts));
		 
		 
		 
		 
		 if (TMath::Abs(o_vars[0]-0.74) < 0.4)
		 {
		     d1->Fill(o_vars[1],o_vars[2]);
		     h6->Fill(o_vars[1]);
		     h6->Fill(o_vars[2]);

		 }
		 
		 if (TMath::Abs(o_vars[1]-0.74) < 0.4)
		 {
		     d1->Fill(o_vars[0],o_vars[2]);
		     h6->Fill(o_vars[0]);
		     h6->Fill(o_vars[2]);

		 }
		 
		 if (TMath::Abs(o_vars[2]-0.74) < 0.4)
		 {
		     d1->Fill(o_vars[0],o_vars[1]);
		     
		     h6->Fill(o_vars[0]);
		     h6->Fill(o_vars[1]);

		 }
		 
		 
		 
	     }
	     
	     
	     
	     if (m_vars.size() == 2)
	     {

		 d2->Fill(m_vars[0],m_vars[1]);
		 
		 if (TMath::Abs(m_vars[0]-0.74) < 0.4)
		 {
	             h7->Fill(m_vars[1]);
		 }
		 
		 if (TMath::Abs(m_vars[1]-0.74) < 0.4)
		 {
	             h7->Fill(m_vars[0]);
		 }
	     
	     }
	     

	     if (m_vars.size() == 3)
	     {
		 
		 
		 //h7->Fill(m_vars[0]);
		 //h7->Fill(m_vars[1]);
		 //h7->Fill(m_vars[2]);
		 
		 
		 
		 if (TMath::Abs(m_vars[0]-0.740) < 0.4)
		 {
		     d2->Fill(m_vars[1],m_vars[2]);
		     h7->Fill(m_vars[1]);
		     h7->Fill(m_vars[2]);

		 }
		 
		 if (TMath::Abs(m_vars[1]-0.740) < 0.4)
		 {
		     d2->Fill(m_vars[0],m_vars[2]);
		     
		     h7->Fill(m_vars[0]);
		     h7->Fill(m_vars[2]);
		 }
		 
		 if (TMath::Abs(m_vars[2]-0.740) < 0.4)
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
    
    auto legend2 = new TLegend(0.64,0.4,0.99,0.78);
    legend2->AddEntry(Yhist,"With outliers","lep");
    legend2->AddEntry(Yhist2,"Without outliers","lep");
    legend2->SetTextSize(0.028);
    legend2->Draw();
    
    Xhist2->SetLineColor(kRed);
    Yhist2->SetLineColor(kRed);
   
    TCanvas *c5 = new TCanvas("c5","c5",800,600);
    
    d3->Draw("Colz");
    
    
    
    
    TCanvas *c6 = new TCanvas("c6","c6",800,600);
    h6->Draw("E1");
    h7->Draw("SAMEE1");
    
    auto legend6 = new TLegend(0.64,0.4,0.99,0.78);
    legend6->AddEntry(h6,"With outliers","lep");
    legend6->AddEntry(h7,"Without outliers","lep");
    legend6->SetTextSize(0.028);
    legend6->Draw();
    
    h7->SetLineColor(kRed);
   
    
    
    h8->SetTitle("Six-track: Pairings; Invariant Mass1 [GeV] ; ");
    
    
}



