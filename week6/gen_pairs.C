#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <tuple>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <vector>
#include <utility>
#include <algorithm>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

TH1F *mass = new TH1F("mass","mass", 1000, 0.2,2);
TH1F *fourmass = new TH1F("4mass","4mass", 1000, 0.8,3.5);

Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

Float_t calc_E(const std::vector<float>& p)
{
    return TMath::Sqrt(mpi*mpi+p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

Double_t calc_invM(const std::vector<float> p1, const std::vector<float> p2)
{
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    
    return TMath::Sqrt(TMath::Power(E1+E2,2) - TMath::Power(p1[0]+p2[0],2)- TMath::Power(p1[1]+p2[1],2)- TMath::Power(p1[2]+p2[2],2));
}

   // Calculate invariant four mass
Double_t calc_fourmass(std::vector<float> p1, std::vector<float> p2,std::vector<float> p3,std::vector<float> p4)
{
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    Double_t E3 = calc_E(p3);
    Double_t E4 = calc_E(p4);
    
    return TMath::Sqrt(TMath::Power(E1+E2+E3+E4,2) - TMath::Power(p1[0]+p2[0]+p3[0]+p4[0],2)- TMath::Power(p1[1]+p2[1]+p3[1]+p4[1],2)- TMath::Power(p1[2]+p2[2]+p3[2]+p4[2],2));
}

//                     ----------------------------------------------------

void gen_pairs()
{
    
    // Create variables to store the data
    
    static Float_t trk_pt[420], trk_eta[420], trk_phi[420], trk_dedx[420], trk_p[420], trk_dxy[420], trk_dz[420];
    static Int_t trk_isK[420], trk_isPi[420], trk_isP[420];
    
    static Float_t ThxR, ThxL, ThyR, ThyL, zPV;
    static Int_t trk_q[420];
    static Int_t ntrk;
    
    // Create branches of the TTree we are loading in
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q); 
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
    
    // Load in a given n-track

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
	    
	    std::vector<float> dxys;
	    std::vector<float> dzs;
	    std::vector<float> dds;
	    std::vector<int> qs;
	    
	    // Define variables
	    
	    for (int i = 0; i < track; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		pt[i] = trk_pt[i];
		p[i] = trk_p[i];
		q[i] = trk_q[i];
		dxy[i] = trk_dxy[i];
		dz[i] = trk_dz[i];
		
		dxys.push_back(TMath::Abs(dxy[i]));
		dzs.push_back(TMath::Abs(dz[i]));
		qs.push_back(q[i]);
		dds.push_back(TMath::Sqrt(dz[i]*dz[i]+dxy[i]*dxy[i]));

	    }
	    
	    std::vector<int> data = {0,1,2,3,4,5};
	    
	    int n = data.size();
	    std::vector<int> indices(n);
            std::iota(indices.begin(), indices.end(), 0);
            std::vector<int> combination(4, 0.0);
	    
	    do 
            {
                int total_charge = 0;
                for (int i = 0; i < 4; ++i) 
	        {
                     combination[i] = data[indices[i]];
	             total_charge += q[indices[i]];
                }
	
	        if (total_charge != 0)
	        {
	            continue;
	        }
	
	    
	    std::vector<float> p1 = {px[combination[0]],py[combination[0]],pz[combination[0]]};
	    std::vector<float> p2 = {px[combination[1]],py[combination[1]],pz[combination[1]]};
	    std::vector<float> p3 = {px[combination[2]],py[combination[2]],pz[combination[2]]};
	    std::vector<float> p4 = {px[combination[3]],py[combination[3]],pz[combination[3]]};
	    
	    float mass1, mass2, mass3, mass4;
    
            if (q[combination[0]] + q[combination[1]] == 0)
	    {
		 mass1 = calc_invM(p1, p2);
                 mass2 = calc_invM(p3, p4);
	
	         if (q[combination[0]] + q[combination[2]] == 0)
                 {
                      mass3 = calc_invM(p1, p3);
                      mass4 = calc_invM(p2, p4);
                 }
		 
                 else
                 {
                      mass3 = calc_invM(p1, p4);
                      mass4 = calc_invM(p2, p3);
                 }
	
            }
            else
	    
           {
               mass1 = calc_invM(p1, p3);
               mass2 = calc_invM(p2, p4);
	
	       mass3 = calc_invM(p1, p4);
               mass4 = calc_invM(p2, p3);

           }
	   
	   if (TMath::Abs(mass1-mrho) < 0.15)
	   {
	       mass->Fill(mass2);
	       if (TMath::Abs(mass2-mrho) < 0.15)
	       {
	           fourmass->Fill(calc_fourmass(p1,p2,p3,p4));
	       }
	   }
	   
	   if (TMath::Abs(mass3-mrho) < 0.15)
	   {
	       mass->Fill(mass4);
	       if (TMath::Abs(mass4-mrho) < 0.15)
	       {
	           fourmass->Fill(calc_fourmass(p1,p2,p3,p4));
	       }
	   }

           }
	   while (std::next_permutation(indices.begin(), indices.end()));
	    

	    
	    
        }
     }
      

    TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
    mass->Draw("E1");
   
    
    TLatex t11;
    t11.SetTextSize(0.06);
    t11.DrawLatex(0.76*3,1, "#font[22]{CMS}");
    
    TLatex t22;
    t22.SetTextSize(0.034);
    t22.DrawLatex(0.86*3,1.2, "#sqrt{s} = 13 #font[22]{TeV}");
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800,600);
    fourmass->Draw("E1");
    

    
    
}



