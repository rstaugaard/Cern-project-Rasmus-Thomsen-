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

TFile *file1 = new TFile("/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/ntuples/data/TOTEM43/110000_nocut.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

TH1F *hallmass = new TH1F("hallmass","hallmass", 100,0.2,1.6);
TH1F *bothmass = new TH1F("bothmass","bothmass", 140,0.2,1.4);
TH1F *offdiagmass = new TH1F("offdiagmass","offdiagmass", 90,0.2,1.4);
TH1F *diagmass = new TH1F("diagmass","diagmass", 90,0.2,1.4);

TH1F *h6 = new TH1F("h6","h6", 90,0.2,1.4);
TH1F *h7 = new TH1F("h7","four-mass", 40,1.2,3.5);

TH1F *h6_2 = new TH1F("h6_2","h6_2", 55,0.2,1.4);
TH1F *h7_2 = new TH1F("h7_2","four-mass", 30,1.2,3.5);

TH1F *h8 = new TH1F("h8","h8", 140,0.2,1.4);
TH1F *h9 = new TH1F("h9","h9", 200,-4,4);
TH2F *d1 = new TH2F("d1","d1",200,0.2,2,200,0.2,2);
TH2F *d2 = new TH2F("d2","d2",100,0.2,2,100,0,2);
TH2F *d3 = new TH2F("d3","d3",100,0.2,2,100,0,2);

TH2F *d4 = new TH2F("d4","d4",100,0.2,2,100,0,3);
TH2F *d5 = new TH2F("d5","d5",100,0,3,100,0,3);

TH1F *h11 = new TH1F("h11","dx", 100, -4,4);
TH1F *h12 = new TH1F("h12","dy", 100, -4,4);
TH1F *h13 = new TH1F("h13","dz", 100, -4,4);

// Struct to safe information of a track more conveniently

struct Particle
{
    float px;
    float py;
    float pz;
    int charge;
    float dz; 
    int idx;
    
    Particle(const std::vector<float>& m, int c, float d, int i) : px(m[0]), py(m[1]), pz(m[2]), charge(c), dz(d), idx(i) {}
};

     //  Handy functions     ----------------------------------------------------
Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

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

     //  fit functions         ----------------------------------------------------

Double_t gauss(float x,Double_t *par)
{
    Double_t value = par[0]*TMath::Gaus(x,par[1],par[2]);
    return value;
}
 
 
Double_t background(Float_t x,Double_t *par)
{
    Double_t value = par[0]*TMath::Exp(-x*par[1]);
    return value;
}

Double_t totalfit(Float_t x,Double_t *par)
{
   Double_t bp[4] = { 1.23189e+04,2.47533e-01,8.05796e-01,3.77877e+00};
   Double_t pars1[4] = {par[0],par[1]};
   Double_t pars2[3] = {par[2],par[3],par[4]};
   //Double_t pars3[3] = {par[7],par[8],par[9]};
   
   Double_t value = background(x,pars1)+gauss(x,pars2); //+gauss(x,pars3);
   return value;
}


Int_t numbins;
 

   // Custom fcn function (function to be minimized)
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
       xval = h6->GetBinCenter(i);
       yval = h6->GetBinContent(i);
       err = h6->GetBinError(i);
       if (yval > 0 && xval > 0.58)
       {
          dx = (totalfit(xval, par) - yval)/err;
          chi2 += (dx*dx);
	  numbins++;
       }
   }
       //std::cout << chi2 << std::endl;     
       f = chi2 ;
}

 // Number of parameters in model
const Int_t numpar = 5;
Int_t num_entries = 0;

   // Minimzer function
void myminimizer(Double_t *par, Double_t *err)
{
    Double_t arglist[numpar] = {1.49602e+03,3.07174e+00,5.30886e+01,7.52843e-01,5.19705e-02};
    
    TMinuit *gMinuit2 = new TMinuit(10);
    gMinuit2->SetFCN(fcn);
    
    Int_t ierflg = 0 ;
    gMinuit2->mnexcm("SET ERR", arglist ,0,ierflg);
    gMinuit2->mnparm(0, "N", arglist[0], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(1, "tau", arglist[1], 0.0001,0, 0, ierflg);
    gMinuit2->mnparm(2, "N",arglist[2], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(3, "mu",arglist[3], 0.0001, 0, 0, ierflg);
    gMinuit2->mnparm(4, "sigma", arglist[4], 0.0001, 0, 0, ierflg);

    gMinuit2->mnexcm("MINOS", arglist , 2,ierflg);
    
    for (int j = 0; j < numpar; j++)
    {
       gMinuit2->GetParameter(j,par[j],err[j]);
    }
}

 // Function to determine if three entries are closer together (in dxy) than 
 // they are to the fourth entry
 
 // Returns true if the sum of distance of the three closest entries is less than 
 // the distance of the fourth closest entry
 
 
bool isThreeEntriesCloser(const std::vector<Particle>& input) 
{
    if (input.size() != 4) 
    {
        std::cerr << "Input vector must have exactly four entries.\n";
        return false;
    }
    

    // Calculate pairwise distances
    std::vector<double> distances;
    for (size_t i = 0; i < input.size(); ++i) 
    {
        for (size_t j = i + 1; j < input.size(); ++j) 
	{
            double distance = std::abs(input[i].dz - input[j].dz);
            distances.push_back(distance);
        }

    }
    
    // Sort distances in ascending order
    std::sort(distances.begin(), distances.end());

    // Check if three smallest distances are closer together than the largest distance
    double sumOfSmallestDistances = distances[0] + distances[1] + distances[2];
    return sumOfSmallestDistances < distances[3];
}

// Function to reduce a combination of n > 4 down to a combination of 4.
// It calculates all given combinations of 4 elements. If a given combination
// Is net-neutral and does not fail the isThreeEntriesCloser test
// Then takes the standard deviation. The returned element will be the combination
// With the lowest standard deviation

std::vector<Particle> reduction(const std::vector<Particle>& data1, const std::vector<Particle>& data2, float ThX, float ThY) 
{
    int n1 = data1.size();
    int n2 = data2.size();
    std::vector<Particle> result;
    
    if (n1 < 2 || n2 < 2) 
    {
        return result;
    }
    
    double minStdDev = 0.2;
    
    for (int i1 = 0; i1 < n1; ++i1) 
    {
        for (int j1 = i1 + 1; j1 < n1; ++j1) 
	{
	    for (int i2 = 0; i2 < n2; ++i2)
	    {
		for (int j2 = i2 + 1; j2 < n2; ++j2)
		{
		    
		    std::vector<Particle> combination = {data1[i1],data1[j1],data2[i2],data2[j2]};
		    
		    // Check if three is closer
		    if (isThreeEntriesCloser(combination) == true)
	            {
	                continue;
	            }
		    
		    if (TMath::Abs(combination[0].px+combination[1].px+combination[2].px+combination[3].px-6500*(ThX)) > 0.1 &&
	    TMath::Abs(combination[0].py+combination[1].py+combination[2].py+combination[3].py+6500*(ThY) > 0.1))
		    {
		        continue;
		    }
		    		    	    
	            // Calculate standard deviation
		    double mean = (combination[0].dz+combination[1].dz+combination[2].dz+combination[3].dz) / 4.0;
		    double sumSquaredDiff = 0;
		    for (int j = 0; j < 4; j++)
		    {
			sumSquaredDiff += std::pow(combination[j].dz - mean, 2);
		    };

		    double stdDev = std::sqrt(sumSquaredDiff / 4.0);

		    // Update result if standard deviation is smaller
		    if (stdDev < minStdDev) 
		    {
			minStdDev = stdDev;
			result = combination;
		    }
		}
	    }
	}
    }

    return result;
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
	
	if (p1.charge + p3.charge == 0)
        {
            mass3 = calc_invM(p1, p3);
            mass4 = calc_invM(p2, p4);
            DeltaDxy3 = TMath::Abs(p1.dz - p3.dz);
            DeltaDxy4 = TMath::Abs(p2.dz - p4.dz);
        }
        else
        {
            mass3 = calc_invM(p1, p4);
            mass4 = calc_invM(p2, p3);
            DeltaDxy3 = TMath::Abs(p1.dz - p4.dz);
            DeltaDxy4 = TMath::Abs(p2.dz - p3.dz);
        }
	
    }
    else
    {
        mass1 = calc_invM(p1, p3);
        mass2 = calc_invM(p2, p4);
        DeltaDxy1 = TMath::Abs(p1.dz - p2.dz);
        DeltaDxy2 = TMath::Abs(p2.dz - p3.dz);
	
	mass3 = calc_invM(p1, p4);
        mass4 = calc_invM(p2, p3);
        DeltaDxy3 = TMath::Abs(p1.dz - p4.dz);
        DeltaDxy4 = TMath::Abs(p2.dz - p3.dz);
    }
    
    d1->Fill(mass1,mass2);
    d1->Fill(mass3,mass4);
    
    std::vector<float> masses;
    
    if (DeltaDxy1 + DeltaDxy2 > DeltaDxy3 + DeltaDxy4)
    {
	masses.push_back(mass3);
        masses.push_back(mass4);
	d2->Fill(mass1,DeltaDxy1);
        d2->Fill(mass2,DeltaDxy2);
        d3->Fill(mass3,DeltaDxy3);
        d3->Fill(mass4,DeltaDxy4);
	      
	if (TMath::Abs(mass3-mrho) < 0.15)
        {
            diagmass->Fill(mass4);
	    h8->Fill(mass4);
        }
	
	if (TMath::Abs(mass1-mrho) < 0.15)
        {
            offdiagmass->Fill(mass2);
	    h8->Fill(mass2);
        }
    }
    else
    {
        masses.push_back(mass1);
        masses.push_back(mass2);
	d2->Fill(mass3,DeltaDxy3);
        d2->Fill(mass4,DeltaDxy4);
	d3->Fill(mass1,DeltaDxy1);
        d3->Fill(mass2,DeltaDxy2);
	
	if (TMath::Abs(mass1-mrho) < 0.15)
        {
            diagmass->Fill(mass2);
	    h8->Fill(mass2);
        }
	
	if (TMath::Abs(mass3-mrho) < 0.15)
        {
            offdiagmass->Fill(mass4);
	    h8->Fill(mass4);
        }
    }
  
    d4->Fill(mass1,DeltaDxy1);
    d4->Fill(mass2,DeltaDxy2);
    d4->Fill(mass3,DeltaDxy3);
    d4->Fill(mass4,DeltaDxy4);
  
    d5->Fill((DeltaDxy1+DeltaDxy2)/2,(DeltaDxy3+DeltaDxy4)/2);
  
    return masses; 
}

struct tree_variables
{
    Float_t px[6];
    Float_t py[6];
    Float_t pz[6];
    Float_t pT[6];
    Int_t q[6];
    Float_t dz[6];
    Float_t eta[6];
    Float_t ThX[1];
    Float_t ThY[1];
    
    Int_t indx[4];
    Float_t red_dz[4];
    Float_t red_px[4];
    Float_t red_py[4];
    Float_t red_pz[4];
};


//                     ----------------------------------------------------

void unfiltered_reduc()
{
    TFile file1("../data_files/reduced_data_dz_p_restrict_.root","RECREATE");

    TTree t1("t1","tree1");


    tree_variables tvars;
    t1.Branch("px",tvars.px,"px[6]/F");
    t1.Branch("py",tvars.py,"py[6]/F");
    t1.Branch("pz",tvars.pz,"pz[6]/F");
    t1.Branch("pT",tvars.pT,"pT[6]/F");
    t1.Branch("q",tvars.q,"q[6]/I");
    t1.Branch("dz",tvars.dz,"dz[6]/F");
    t1.Branch("eta",tvars.eta,"eta[6]/F");
    t1.Branch("ThX",tvars.ThX,"ThX[1]/F");
    t1.Branch("ThY",tvars.ThY,"ThY[1]/F");

    t1.Branch("indx",tvars.indx,"indx[4]/I");
    t1.Branch("red_dz",tvars.red_dz,"dx[4]/F");
    t1.Branch("red_px",tvars.red_px,"red_px[4]/F");
    t1.Branch("red_py",tvars.red_py,"red_py[4]/F");
    t1.Branch("red_pz",tvars.red_pz,"red_pz[4]/F");

    
    
    int nthreads = 4;
    ROOT::EnableImplicitMT(nthreads);
    
    TH1F *h1 = new TH1F("h1","h1", 200,-4,4);
    TH1F *h2 = new TH1F("h2","h2", 200,-4,4);
    TH1F *h3 = new TH1F("h3","h3", 200,-4,4);
    
    TH2F *d1 = new TH2F("d1","d1",100,0.2,1.6,100,0.2,1.6);
    
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
    
    //my_tree->SetBranchAddress("trk_dxy",&trk_dxy);
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
    
    int count = 0;

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
		
		dxy[i] = trk_dz[i];
		dz[i] = trk_dz[i];
		//dxys.push_back(TMath::Abs(dxy[i]));
		dzs.push_back(TMath::Abs(dz[i]));
		qs.push_back(q[i]);
		//dds.push_back(TMath::Sqrt(dz[i]*dz[i]+dxy[i]*dxy[i]));
		
		
		 tvars.px[i] = px[i];
	         tvars.py[i] = py[i];
	         tvars.pz[i] = pz[i];
	         tvars.pT[i] = pt[i];
	         tvars.q[i] = q[i];
	         tvars.dz[i] = dz[i];
	         tvars.eta[i] = eta[i];
	    }
	   
	    std::vector<float> p1 = {px[0],py[0],pz[0]};
	    std::vector<float> p2 = {px[1],py[1],pz[1]};
	    std::vector<float> p3 = {px[2],py[2],pz[2]};
	    std::vector<float> p4 = {px[3],py[3],pz[3]};
	    std::vector<float> p5 = {px[4],py[4],pz[4]};
	    std::vector<float> p6 = {px[5],py[5],pz[5]};
	    
	    Particle particle1(p1,qs[0], TMath::Abs(dzs[0]), 0);
	    Particle particle2(p2,qs[1], TMath::Abs(dzs[1]), 1);
	    Particle particle3(p3,qs[2], TMath::Abs(dzs[2]), 2);
	    Particle particle4(p4,qs[3], TMath::Abs(dzs[3]), 3);
	    Particle particle5(p5,qs[4], TMath::Abs(dzs[4]), 4);
	    Particle particle6(p6,qs[5], TMath::Abs(dzs[5]), 5);
	    
	    std::vector<Particle> allparticles= {particle1,particle2,particle3,particle4,particle5,particle6};
	    std::vector<Particle> dataplus;
	    std::vector<Particle> dataminus;
	    
	    for (int i = 0; i < track; i++)
	    {
	        if(q[i] > 0)
		{
		    dataplus.push_back(allparticles[i]);
		}
		
		else
		{
		    dataminus.push_back(allparticles[i]);
		}
	    }
	    
	    // Reduce track to 4-track	
	    std::vector<Particle> pairs = reduction(dataplus,dataminus, ThxL+ThxR, ThyL+ThyR);
	    
	    // If no 4-track is found, then continue to next data point
	    if (pairs.size() == 0)
	    {
	        continue;
	    }

	   
	    if (count < 20)
	    {
	    // Check if things work as intended 
	        std::cout<< dz[0] << " , " << dz[1] << " , " << dz[2] << " , " << dz[3] << " , " << dz[4] << " , " << dz[5] << std::endl;
	        std::cout<< qs[0] << " , " << qs[1] << " , " << qs[2] << " , " << qs[3] << " , " << qs[4] << "  ,  " << qs[5]<< std::endl;
	        std::cout<< "Reduction;" << " " << pairs[0].dz << " , " << pairs[1].dz << " , " << pairs[2].dz << " , " << pairs[3].dz << std::endl;
		 std::cout<< pairs[0].idx << " , " << pairs[1].idx << " , " << pairs[2].idx << " , " << pairs[3].idx << std::endl;
	    }
	    
	    count += 1;
	   
	    //
	  
	    tvars.ThX[0] = ThxL+ThxR;
            tvars.ThY[0] = ThyL+ThyR;
	    
	    for (int m = 0; m < 4; m++)
	    {
	        tvars.indx[m] = pairs[m].idx;
	        tvars.red_dz[m] = pairs[m].dz;
	        tvars.red_px[m] = pairs[m].px;
	        tvars.red_py[m] = pairs[m].py;
	        tvars.red_pz[m] = pairs[m].pz;
	    }
	    

	    
	    t1.Fill();
	   

	    // Finally pair up the 4-track into two 2-tracks
	    std::vector<float> mass = pair_up(pairs);
	    
	    if (mass.size() == 0)
	    {
	        continue;
	    }

	    
	    // Cut on the 2-mass, so we only consider interesting pairs
	    if (TMath::Abs(mass[0]-mrho) < 0.15)
	    {
	        h6->Fill(mass[1]);
		if (TMath::Abs(mass[1]-mrho) < 0.15)
		{
		    h7->Fill(calc_fourmass(pairs));
		    
		}
		
		
	    }
	    
	    if (TMath::Abs(pairs[0].px+pairs[1].px+pairs[2].px+pairs[3].px-6500*(ThxL+ThxR))< 0.1 && TMath::Abs(pairs[0].py+pairs[1].py+pairs[2].py+pairs[3].py+6500*(ThyL+ThyR)) < 0.1)
	    {
	        if (TMath::Abs(mass[0]-mrho) < 0.15)
	        {
	            h6_2->Fill(mass[1]);
		    if (TMath::Abs(mass[1]-mrho) < 0.15)
		    {
		        h7_2->Fill(calc_fourmass(pairs));
		    
		    } 
		
		
	        }
	    }
	    
        }
     }
      
      
    TCanvas *c1 = new TCanvas("c1","c1",1800,800);
    c1->Divide(3,1);

    gStyle->SetPalette(kCividis);
    gStyle->SetOptStat(false);
    c1->cd(1); h1->Draw();
    c1->cd(2); h2->Draw();
    c1->cd(3); h3->Draw();
    

    
    
    
    TCanvas *c6 = new TCanvas("c6","c6",800,600);
    h6->Draw("E1");
    
    

    Double_t params[7];
    Double_t errs[7];
    myminimizer(params,errs);
    Double_t x[1400], y[1400];
    
    for (Int_t i=0;i< 1400;i++) 
    {
        x[i] = 0.55+i*0.005;
        y[i] = totalfit(x[i],params);
     }

    auto g2 = new TGraph(1400,x,y);
    g2->SetLineWidth(2);
    g2->SetLineColor(kRed);
    g2->Draw("SAME");
    
    Double_t x2[500], y2[500];
    for (Int_t i=0;i< 400;i++) 
    {
        x2[i] = 0.55+i*0.001;
	Double_t params2[3] = {params[2],params[3],params[4]};
        y2[i] = gauss(x2[i],params2);
     }

    auto g4 = new TGraph(400,x2,y2);
    g4->SetLineWidth(2);
    g4->SetLineColor(kGreen);
    g4->Draw("SAME");
    
    for (Int_t i=0;i< 300;i++) 
    {
        Double_t back_params[2] = {params[0],params[1]};
        y[i] = background(x[i],back_params);
     }

    auto g5 = new TGraph(1400,x,y);
    g5->SetLineWidth(2);
    g5->SetLineColor(kBlue);
    g5->Draw("SAME");
    
     auto legend6 = new TLegend(0.64,0.46,0.99,0.8);
    //legend->SetHeader("Fit","C");
    legend6->AddEntry(h6,"Data","lep");
    legend6->AddEntry(g5,"Background: Nexp(-tx)","l");
    //legend6->AddEntry(g3,"#splitline{Kaon peak:}{#mu = 0.50 #pm 0.02 GeV}","l");
    legend6->AddEntry(g4,"#splitline{Rho peak:}{#splitline{ #mu = 0.755 #pm 0.009 GeV}{#sigma = 0.05 #pm 0.01 GeV}} ","l");
    legend6->AddEntry(g2,"Total fit: Chi2 / NDof : 73 / 55","l");
    legend6->SetTextSize(0.028);
    legend6->Draw();
    
    h6->SetTitle(" ; Inv. Mass [GeV]; ");
    
    TLatex s1;
    s1.SetTextSize(0.06);
    s1.DrawLatex(0.82*1.4,218, "#font[22]{CMS}");
    
    TLatex s2;
    s2.SetTextSize(0.034);
    s2.DrawLatex(0.87*1.4,241, "#sqrt{s} = 13 #font[22]{TeV}");
   
    TCanvas *c7 = new TCanvas("c7","c7",1600,600);
    c7->Divide(2,1);
    c7->cd(1); d2->Draw("Colz");
    
    TLatex s27;
    s27.SetTextSize(0.06);
    s27.DrawLatex(0.3,2.02, "#font[22]{CMS}");
    
    TLatex s28;
    s28.SetTextSize(0.034);
    s28.DrawLatex(0.85*2,2.02, "#sqrt{s} = 13 #font[22]{TeV}");
    
    
    c7->cd(2); d3->Draw("Colz");
    d2->SetTitle("#Deltadz large; Inv. Mass; #Delta|dz|");
    d3->SetTitle("#Deltadz small; Inv. Mass; #Delta|dz|");
    
    TLatex s7;
    s7.SetTextSize(0.06);
    s7.DrawLatex(0.3,2.02, "#font[22]{CMS}");
    
    TLatex s8;
    s8.SetTextSize(0.034);
    s8.DrawLatex(0.85*2,2.02, "#sqrt{s} = 13 #font[22]{TeV}");
    
  
  
    std::cout << "Number of bins: " << "  "  << numbins << std::endl;
    std::cout << "mu: " << "  "  << params[3] << "+-" << errs[3] << std::endl;
    std::cout << "sigma: " << "  "  << params[4] << "+-" << errs[4] << std::endl;
    std::cout << "Significans: " << "  "  << params[2]/errs[2]  << std::endl;
    
    TCanvas *c8 = new TCanvas("c8","c8",800,600);
    h7->Draw("E1");
    
    TLatex s3;
    s3.SetTextSize(0.06);
    s3.DrawLatex(0.82*3,222, "#font[22]{CMS}");
    
    TLatex s4;
    s4.SetTextSize(0.034);
    s4.DrawLatex(0.87*3,260, "#sqrt{s} = 13 #font[22]{TeV}");
    
    h7->SetTitle(" ; Inv. 4-Mass [GeV] ; ");
    
    //h9->SetLineColor(kRed);
    
    TCanvas *c9 = new TCanvas("c9","c9",800,600);
    h8->Draw("E1");
    
    TLatex s88;
    s88.SetTextSize(0.06);
    s88.DrawLatex(0.82*3,222, "#font[22]{CMS}");
    
    TLatex s89;
    s89.SetTextSize(0.034);
    s89.DrawLatex(0.87*3,260, "#sqrt{s} = 13 #font[22]{TeV}");
    
    
    h8->SetTitle(" ; Inv. Mass [GeV];");
    
    TCanvas *c10 = new TCanvas("c10","c10",1600,600);
    c10->Divide(2,1);
    c10->cd(1); d4->Draw("Colz");
    
    TLatex s25;
    s25.SetTextSize(0.06);
    s25.DrawLatex(0.3,3.02, "#font[22]{CMS}");
    
    TLatex s26;
    s26.SetTextSize(0.034);
    s26.DrawLatex(0.84*3,3.02, "#sqrt{s} = 13 #font[22]{TeV}");
    
    c10->cd(2); d5->Draw("Colz");
    
    d4->SetTitle(" ; Inv. Mass ; #Delta|dz|");
    d5->SetTitle(" ;#Delta|dz_{1}| ;#Delta|dz_{2}|");
    
    TLatex s5;
    s5.SetTextSize(0.06);
    s5.DrawLatex(0.3,3.02, "#font[22]{CMS}");
    
    TLatex s6;
    s6.SetTextSize(0.034);
    s6.DrawLatex(0.84*3,3.02, "#sqrt{s} = 13 #font[22]{TeV}");
    
    
    TCanvas *c11 = new TCanvas("c11", "c11", 1600,600);
    c11->Divide(3,1);
    c11->cd(1); h11->Draw();
    c11->cd(2); h12->Draw();
    c11->cd(3); h13->Draw();
    
    TCanvas *c12 = new TCanvas("c12", "c12", 1600,600);
    c12->Divide(2,1);
    c12->cd(1); h6_2->Draw("E1");
    c12->cd(2); h7_2->Draw("E1");
    
    h6_2->SetTitle(" Two-mass; Inv. Mass; ");
    h7_2->SetTitle("Four-mass; Inv. 4-mass; ");
        
	
	
    TCanvas *allmass = new TCanvas("allmass","allmass",800,600);
    bothmass->Draw("E1");
    h6->Draw("SAMEE1");
    
    bothmass->SetLineColor(kRed);
    
    TCanvas *diag = new TCanvas("diag","diag",1600,600);
    diag->Divide(2,1);
    diagmass->Scale(1./diagmass->Integral(), "width");
    diag->cd(1); diagmass->Draw();
    diagmass->SetLineColor(kBlue);
    diagmass->SetMarkerStyle(25);
    diagmass->SetMarkerColorAlpha(kBlue,0.5);
    diagmass->SetTitle("Selected ; Inv. Mass [GeV] ; ");
    
    
    TLatex t2;
    t2.SetTextSize(0.034);
    t2.DrawLatex(0.86*1.4,2.32, "#sqrt{s} = 13 #font[22]{TeV}");
    
    
    diagmass->Sumw2();
    offdiagmass->Scale(1./offdiagmass->Integral(), "width");
    diag->cd(2); offdiagmass->Draw("SAME");
    offdiagmass->SetMarkerStyle(24);
    offdiagmass->SetMarkerColorAlpha(kRed,0.5);
    offdiagmass->Sumw2();
    offdiagmass->SetTitle("Not selected ; Inv. Mass [GeV]; ");
    
    offdiagmass->SetLineColor(kRed);
    
    TLatex t11;
    t11.SetTextSize(0.06);
    t11.DrawLatex(0.76*1.4,1.76, "#font[22]{CMS}");
    
    TLatex t22;
    t22.SetTextSize(0.034);
    t22.DrawLatex(0.86*1.4,2.1, "#sqrt{s} = 13 #font[22]{TeV}");
    
    
    
    //offdiagmass->Draw("E1");
    
    
}



