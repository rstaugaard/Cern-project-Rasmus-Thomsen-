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

TFile *file1 = new TFile("../data_files/reduced_data_dxy.root","read");
TTree *my_tree = (TTree*)file1->Get("t1");

TH1F *looper = new TH1F("looper","looper",100,0,2);
TH1F *alltwo = new TH1F("alltwo", "alltwo", 200, 0.2, 1.6);
TH1F *htwo = new TH1F("htwo", "htwo", 200, 0.2, 1.6);
TH1F *hfour = new TH1F("hfour", "hfour", 200, 1.1, 3.5);

TH1F *htwo_cut = new TH1F("htwo_cut","htwo_cut",80,0.2,1.8);
TH1F *hfour_cut = new TH1F("hfour_cut","hfour_cut",50,1.2,3.5);

TH2F *corr1 = new TH2F("corr1", "corr1",300,0.2,1.6,300,0,5);
TH2F *corr2 = new TH2F("corr2", "corr2",300,0.2,1.6,300,0,5);
TH2F *corr3 = new TH2F("corr3", "corr3",300,0.2,1.6,300,0,2.5);
TH2F *corr4 = new TH2F("corr4", "corr4",300,0.2,1.6,300,0,2.5);

     //  Handy functions     ----------------------------------------------------
Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

// Struct to safe information of a track more conveniently

struct Particle
{
    float px;
    float py;
    float pz;
    int charge;
    float dxy;
    float dz; 
    int idx;
    
    Particle(const std::vector<float>& m, int c, float d1,float d2, int i) : px(m[0]), py(m[1]), pz(m[2]), charge(c), dxy(d1), dz(d2), idx(i) {}
};

  //Calculate energy assuming pion 
Float_t calc_E(Particle p1)
{
    return TMath::Sqrt(mpi*mpi+p1.px*p1.px+p1.py*p1.py+p1.pz*p1.pz);
}
 
 // Calculate invariant two-mass
Double_t calc_invM(Particle p1, Particle p2)
{
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    
    return TMath::Sqrt(TMath::Power(E1+E2,2) - TMath::Power(p1.px+p2.px,2)- TMath::Power(p1.py+p2.py,2)- TMath::Power(p1.pz+p2.pz,2));
}

   // Calculate invariant four mass
Double_t calc_fourmass(std::vector<Particle> pairs)
{
    Particle p1 = pairs[0];
    Particle p2 = pairs[1];
    Particle p3 = pairs[2];
    Particle p4 = pairs[3];
    
    Double_t E1 = calc_E(p1);
    Double_t E2 = calc_E(p2);
    Double_t E3 = calc_E(p3);
    Double_t E4 = calc_E(p4);
    
    return TMath::Sqrt(TMath::Power(E1+E2+E3+E4,2) - TMath::Power(p1.px+p2.px+p3.px+p4.px,2)- TMath::Power(p1.py+p2.py+p3.py+p4.py,2)- TMath::Power(p1.pz+p2.pz+p3.pz+p4.pz,2));
}

Float_t vector_sum(Particle p1, Particle p2)
{
    Float_t val = TMath::Sqrt(std::pow(p1.px+p2.px,2)+std::pow(p1.py+p2.py,2)+std::pow(p1.pz+p2.pz,2));
    return val;
}

// Function to determine the two different ways to pair up (construct two-mass) of a net-neutral four track
// It will return the pair of two (of the two) two-masses, which have the smallest difference in
// some variable (usually dz)

std::vector<float> pair_up(std::vector <Particle> pairs)
{
    Particle p1 = pairs[0];
    Particle p2 = pairs[1];
    Particle p3 = pairs[2];
    Particle p4 = pairs[3];
    
    float mass1, mass2, mass3, mass4;
    float DeltaDxy1, DeltaDxy2, DeltaDxy3, DeltaDxy4;
    
    if (p1.charge + p2.charge == 0)
    {
        mass1 = calc_invM(p1, p2);
        mass2 = calc_invM(p3, p4);
        DeltaDxy1 = TMath::Abs(p1.dz - p2.dz);
        DeltaDxy2 = TMath::Abs(p2.dz - p3.dz);
	
        looper->Fill(vector_sum(p1,p2)/mpi);
	looper->Fill(vector_sum(p3,p4)/mpi);
	
	if (p1.charge + p3.charge == 0)
        {
            mass3 = calc_invM(p1, p3);
            mass4 = calc_invM(p2, p4);
            DeltaDxy3 = TMath::Abs(p1.dz - p3.dz);
            DeltaDxy4 = TMath::Abs(p2.dz - p4.dz);
	    	
            looper->Fill(vector_sum(p1,p3)/mpi);
	    looper->Fill(vector_sum(p2,p4)/mpi);
	    
	    

        }
        else
        {
            mass3 = calc_invM(p1, p4);
            mass4 = calc_invM(p2, p3);
            DeltaDxy3 = TMath::Abs(p1.dz - p4.dz);
            DeltaDxy4 = TMath::Abs(p2.dz - p3.dz);
	    
            looper->Fill(vector_sum(p1,p4)/mpi);
	    looper->Fill(vector_sum(p2,p3)/mpi);
        }
	
    }
    else
    {
        mass1 = calc_invM(p1, p3);
        mass2 = calc_invM(p2, p4);
        DeltaDxy1 = TMath::Abs(p1.dz - p2.dz);
        DeltaDxy2 = TMath::Abs(p2.dz - p3.dz);
	
	looper->Fill(vector_sum(p1,p3)/mpi);
	looper->Fill(vector_sum(p2,p4)/mpi);
	
	mass3 = calc_invM(p1, p4);
        mass4 = calc_invM(p2, p3);
        DeltaDxy3 = TMath::Abs(p1.dz - p4.dz);
        DeltaDxy4 = TMath::Abs(p2.dz - p3.dz);
		
        looper->Fill(vector_sum(p1,p4)/mpi);
	looper->Fill(vector_sum(p2,p3)/mpi);
    }
    
    std::vector<float> masses;
    
    if (DeltaDxy1 + DeltaDxy2 > DeltaDxy3 + DeltaDxy4)
    {
	masses.push_back(mass3);
        masses.push_back(mass4);
	      
	if (TMath::Abs(mass3-mrho) < 0.15)
        {
	    alltwo->Fill(mass4);
        }
	
	if (TMath::Abs(mass1-mrho) < 0.15)
        {
	    alltwo->Fill(mass2);
        }
    }
    else
    {
        masses.push_back(mass1);
        masses.push_back(mass2);
	
	if (TMath::Abs(mass1-mrho) < 0.15)
        {
	    alltwo->Fill(mass2);
        }
	
	if (TMath::Abs(mass3-mrho) < 0.15)
        {
	    alltwo->Fill(mass4);
        }
    }
  
    return masses; 
}

//bool IsDivisible(std::vector<Particle> particles)
//{
//    std::vector<float> dzs = {particles[0].dz,particles[1].dz,particles[2].dz,particles[3].dz}; 
//    std::sort(dzs.begin(), dzs.end());
    
//    if 

//}



void unfiltered_analysis()
{
    static Float_t pT[8], eta[8], px[8], py[8], pz[8], dxy[8], dz[8], q[8];
    
    static Float_t ThX, ThY;

    static Float_t red_px[6], red_py[6], red_pz[6], red_dxy[6], red_dz[6], red_q[6];
    static Int_t indx[6];
    
    // Create branches of the TTree we are loading in

    my_tree->SetBranchAddress("pT",&pT);
    my_tree->SetBranchAddress("eta",&eta);
    my_tree->SetBranchAddress("px",&px);
    my_tree->SetBranchAddress("py",&py); 
    my_tree->SetBranchAddress("pz",&pz);
    my_tree->SetBranchAddress("dxy",&dxy); 
    my_tree->SetBranchAddress("dz",&dz);
    my_tree->SetBranchAddress("q",&q); 
    my_tree->SetBranchAddress("ThX",&ThX);
    my_tree->SetBranchAddress("ThY",&ThY);
    
    
    my_tree->SetBranchAddress("indx",&indx);
    my_tree->SetBranchAddress("red_dxy",&red_dxy);
    my_tree->SetBranchAddress("red_dz",&red_dz);
    my_tree->SetBranchAddress("red_px",&red_px);
    my_tree->SetBranchAddress("red_py",&red_py);
    my_tree->SetBranchAddress("red_pz",&red_pz);
    my_tree->SetBranchAddress("red_q",&red_q);
   
    int numEnt = my_tree->GetEntries();
    
    int count = 0;

    for (int irow = 0; irow< numEnt; irow++)
    {    
        my_tree->GetEntry(irow);
	
        std::vector<float> p1 = {px[0],py[0],pz[0]};
        std::vector<float> p2 = {px[1],py[1],pz[1]};
        std::vector<float> p3 = {px[2],py[2],pz[2]}; 
        std::vector<float> p4 = {px[3],py[3],pz[3]}; 
        std::vector<float> p5 = {px[4],py[4],pz[4]};
        std::vector<float> p6 = {px[5],py[5],pz[5]};
    
        Particle particle1(p1, q[0], TMath::Abs(dxy[0]), TMath::Abs(dz[0]), 0);
        Particle particle2(p2, q[1], TMath::Abs(dxy[1]), TMath::Abs(dz[1]), 1);
        Particle particle3(p3, q[2], TMath::Abs(dxy[2]), TMath::Abs(dz[2]), 2);                         	    
        Particle particle4(p4, q[3], TMath::Abs(dxy[3]), TMath::Abs(dz[3]), 3);
        Particle particle5(p5, q[4], TMath::Abs(dxy[4]), TMath::Abs(dz[4]), 4);
        Particle particle6(p6, q[5], TMath::Abs(dxy[5]), TMath::Abs(dz[5]), 5);
	
	float rpx[4], rpy[4], rpz[4], rdxy[4], rdz[4], rq[4];
		
	for (int i = 0; i < 4; i++)
	{
	    rpx[i] = red_px[i]; 
	    rpy[i] = red_py[i];
	    rpz[i] = red_pz[i];
	    rdxy[i] = red_dxy[i];
	    rdz[i] = red_dz[i];
	    rq[i] = red_q[i];
	}
	
        if (count < 20)
	{
	    std::cout << rpx[0] << " , " << rpy[0] << " , " << rpz[0] << " | " << px[indx[0]] << " , "<< py[indx[0]] << " , "<< pz[indx[0]] << std::endl;
	}
	
	count += 1;
	
	
	
	std::vector<float> red_p1 = {rpx[0],rpy[0],rpz[0]};
	std::vector<float> red_p2 = {rpx[1],rpy[1],rpz[1]}; 
	std::vector<float> red_p3 = {rpx[2],rpy[2],rpz[2]}; 
	std::vector<float> red_p4 = {rpx[3],rpy[3],rpz[3]};  
	
	Particle red_particle1(red_p1, rq[0],TMath::Abs(rdxy[0]), TMath::Abs(rdz[0]), 0);
	Particle red_particle2(red_p2, rq[1],TMath::Abs(rdxy[1]), TMath::Abs(rdz[1]), 1);
	Particle red_particle3(red_p3, rq[2],TMath::Abs(rdxy[2]), TMath::Abs(rdz[2]), 2);
	Particle red_particle4(red_p4, rq[3],TMath::Abs(rdxy[3]), TMath::Abs(rdz[3]), 3);
	
        //if (TMath::Abs(red_particle1.px+red_particle2.px+red_particle3.px+red_particle4.px-6500*(ThX)) > 0.25 || TMath::Abs(red_particle1.py+red_particle2.py+red_particle3.py+red_particle4.py+6500*(ThY)) > 0.25)
        //{
	//     continue;
        //}
    
        std::vector<Particle> pairs = {red_particle1, red_particle2, red_particle3, red_particle4};
    
        std::vector<Particle> allpairs = {particle1, particle2, particle3, particle4, particle5, particle6};
	
        std::vector<float> mass = pair_up(pairs);
	
        if (mass.size() == 0)
        {
            continue;
        }
	
	corr1->Fill(mass[0],red_particle1.dz+red_particle2.dz+red_particle3.dz+red_particle4.dz);
	corr1->Fill(mass[1],red_particle1.dz+red_particle2.dz+red_particle3.dz+red_particle4.dz);
	
	corr2->Fill(mass[0],red_particle1.dxy+red_particle2.dxy+red_particle3.dxy+red_particle4.dxy);
	corr2->Fill(mass[1],red_particle1.dxy+red_particle2.dxy+red_particle3.dxy+red_particle4.dxy);
	
	corr3->Fill(mass[0],TMath::Abs(red_particle1.px+red_particle2.px+red_particle3.px+red_particle4.px-6500*(ThX)));
	corr3->Fill(mass[1],TMath::Abs(red_particle1.px+red_particle2.px+red_particle3.px+red_particle4.px-6500*(ThX)));
	
	corr4->Fill(mass[0],TMath::Abs(red_particle1.py+red_particle2.py+red_particle3.py+red_particle4.py+6500*(ThY)));
	corr4->Fill(mass[1],TMath::Abs(red_particle1.py+red_particle2.py+red_particle3.py+red_particle4.py+6500*(ThY)));
	
    
        if (TMath::Abs(mass[0]-mrho) < 0.15)
        {
            htwo->Fill(mass[1]);
	    
	    
            if (TMath::Abs(mass[1]-mrho) < 0.15)
            { 
                hfour->Fill(calc_fourmass(pairs));
            }
    
        }

        if (TMath::Abs(red_particle1.py+red_particle2.py+red_particle3.py+red_particle4.py+6500*(ThY)) > 0.2 ||
	    TMath::Abs(red_particle1.px+red_particle2.px+red_particle3.px+red_particle4.px-6500*(ThX)) > 0.2)
        {
            continue;
        }
	
	if (pT[indx[0]]+pT[indx[1]]+pT[indx[2]]+pT[indx[3]] > 1.7 || pT[indx[0]]+pT[indx[1]]+pT[indx[2]]+pT[indx[3]] < 1)
	{
	    continue;
	}
	
	if (red_particle1.dz+red_particle2.dz+red_particle3.dz+red_particle4.dz > 2 || red_particle1.dz+red_particle2.dz+red_particle3.dz+red_particle4.dz < 1)
	{
	    continue;
	}
	
	if (red_particle1.dxy+red_particle2.dxy+red_particle3.dxy+red_particle4.dxy > 0.6)
	{
	    continue;
	}	
	
	
	if (TMath::Abs(mass[0]-mrho) < 0.1)
        {
            htwo_cut->Fill(mass[1]);
	    
	    
            if (TMath::Abs(mass[1]-mrho) < 0.1)
            { 
                hfour_cut->Fill(calc_fourmass(pairs));
            }
    
        }
	
	
	
	
    
    }



    TCanvas *c1 = new TCanvas("c1", "alltwo", 800,600);
    alltwo->Draw("E1");

    TCanvas *c2 = new TCanvas("c2", "htwo", 800,600);
    //htwo->Draw("E1");
    htwo_cut->Draw("SAMEE1");
    
    htwo_cut->SetLineColor(kRed);

    TCanvas *c3 = new TCanvas("c3", "hfour", 800,600);
    //hfour->Draw("E1");
    hfour_cut->Draw("SAMEE1");
    
    hfour_cut->SetLineColor(kRed);
    
    TCanvas *c4 = new TCanvas("c4","looper",800,600);
    looper->Draw("E1");
    
    TCanvas *c5 = new TCanvas("c5","corr1",800,600);
    corr1->Draw("Colz");

    TCanvas *c6 = new TCanvas("c6","corr2",800,600);
    corr2->Draw("Colz");
    
    TCanvas *c7 = new TCanvas("c7","corr3",800,600);
    corr3->Draw("Colz");
    
    TCanvas *c8 = new TCanvas("c8","corr4",800,600);
    corr4->Draw("Colz");

}

