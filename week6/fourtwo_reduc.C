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

TH1F *h6 = new TH1F("h6","h6", 70,0.2,1.4);
TH1F *h7 = new TH1F("h7","four-mass", 40,0.8,3);

TH1F *h6_2 = new TH1F("h6_2","h6_2", 55,0.2,1.4);
TH1F *h7_2 = new TH1F("h7_2","four-mass", 30,0.8,3);

TH1F *h8 = new TH1F("h8","h8", 200,0.2,1.4);
TH1F *h9 = new TH1F("h9","h9", 200,-4,4);
TH2F *d1 = new TH2F("d1","d1",200,0.2,2,200,0.2,2);
TH2F *d2 = new TH2F("d2","d2",100,0.2,2,100,0,2);
TH2F *d3 = new TH2F("d3","d3",100,0.2,2,100,0,2);

TH2F *d4 = new TH2F("d4","d4",100,0.2,2,100,0,3);
TH2F *d5 = new TH2F("d5","d5",100,0,3,100,0,3);

TH1F *h11 = new TH1F("h11","dx", 100, -2,2);
TH1F *h12 = new TH1F("h12","dy", 100, -2,2);
TH1F *h13 = new TH1F("h13","dz", 100, -3,3);

     //  Handy functions     ----------------------------------------------------
Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

  //Calculate energy assuming pion 
Float_t calc_E(const std::vector<float>& p)
{
    return TMath::Sqrt(mpi*mpi+p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}
 
 // Calculate invariant two-mass
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
 
 
bool isThreeEntriesCloser(const std::vector<float>& input) 
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
            double distance = std::abs(input[i] - input[j]);
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

std::vector<float> reduction(const std::vector<float>& data, const std::vector<int>& q) 
{
    int n = data.size();
    if (n < 4) 
    {
        throw std::invalid_argument("Input vector must have at least four elements.");
    }

    std::vector<float> result;
    double minStdDev = std::numeric_limits<double>::max();

    // Generate all combinations of four elements from the vector
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<float> combination(4, 0.0);

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
	
	if (isThreeEntriesCloser(combination) == true)
	{
	    continue;
	}
	

        // Calculate standard deviation
        double mean = static_cast<double>(std::accumulate(combination.begin(), combination.end(), 0.0)) / 4.0;
        double sumSquaredDiff = std::accumulate(combination.begin(), combination.end(), 0.0, [mean](double sum, float num) 
	{
            return sum + std::pow(num - mean, 2);
        });
        
	double stdDev = std::sqrt(sumSquaredDiff / 4.0);

        // Update result if standard deviation is smaller
        if (stdDev < minStdDev) 
	{
            minStdDev = stdDev;
            result = combination;
        }
    } 
    while (std::next_permutation(indices.begin(), indices.end()));

    return result;
}

// Function to check whether elements are close w.r.t some given variable (usually dxy/dz)
// If the elements are closer than some given threshold, then it returns true

bool areElementsClose(const std::vector<float>& data, float threshold) 
{
    for (int i = 0; i < data.size(); i++) 
    {
        for (int j = i; j < data.size(); j++)
	{
	    float diff = std::abs(data[i] - data[j]);
            if (diff > threshold) 
	    {
                return false; // If the difference exceeds the threshold, elements are not close.
            }
	}
    }
    
    return true; // All elements are close enough.
}

// Find the original indices of a given dxy value

std::vector<int> findIndices(const std::vector<float>& dxy, const std::vector<float>& data)
{
    std::vector<int> indxs;
    for (int i = 0; i < data.size(); i++) // Loop through dxy
    {
        for(int j = 0; j < dxy.size(); j++) // Loop through 4-pair
	{
	    if (data[i] == dxy[j])  //  If the dxy equals reduced pair then give index
	    {
	        indxs.push_back(j);
		break;
	    }
	}
    }
    
    return indxs;
}

// Working idea 
bool isSplitThreeTrack(std::vector<float> picked_dxys,std::vector<float> remaining_dxys)
{
    int count = 0;
    std::vector<float> distance1;
    std::vector<float> distance2;
    for(int i = 0; i < picked_dxys.size(); i++)
    {
	for (int j = 0; j < picked_dxys.size(); j++)
	{
	    distance1.push_back(TMath::Abs(picked_dxys[i]-picked_dxys[j]));
	}
	
        distance2.push_back(TMath::Abs(picked_dxys[i]-remaining_dxys[0]));

	

    }
    
    std::sort(distance2.begin(), distance2.end());
    
    float sum = 0;
    for (int i = 0; i< distance1.size(); i++)
    {
        sum += distance1[i];
    } 
    
    if (2*sum/distance1.size() > distance1[0]+distance2[0])
    {
        return true;
    }
    
    else
    {
        return false;
    }
}





// Struct to safe information of a track more conveniently

struct Particle
{
    std::vector<float> momenta;
    int charge;
    float dz; 
    
    Particle(const std::vector<float>& m, int c, float d) : momenta(m), charge(c), dz(d) {}
};

// Function to determine the two different ways to pair up (construct two-mass) of a net-neutral four track
// It will return the pair of two (of the two) two-masses, which have the smallest difference in
// some variable (usually dz)

std::vector<float> pair_up(Particle p1, Particle p2, Particle p3, Particle p4)
{
    float mass1, mass2, mass3, mass4;
    float DeltaDxy1, DeltaDxy2, DeltaDxy3, DeltaDxy4;
    
    if (p1.charge + p2.charge == 0)
    {
        mass1 = calc_invM(p1.momenta, p2.momenta);
        mass2 = calc_invM(p3.momenta, p4.momenta);
        DeltaDxy1 = TMath::Abs(p1.dz - p2.dz);
        DeltaDxy2 = TMath::Abs(p2.dz - p3.dz);
	
	if (p1.charge + p3.charge == 0)
        {
            mass3 = calc_invM(p1.momenta, p3.momenta);
            mass4 = calc_invM(p2.momenta, p4.momenta);
            DeltaDxy3 = TMath::Abs(p1.dz - p3.dz);
            DeltaDxy4 = TMath::Abs(p2.dz - p4.dz);
        }
        else
        {
            mass3 = calc_invM(p1.momenta, p4.momenta);
            mass4 = calc_invM(p2.momenta, p3.momenta);
            DeltaDxy3 = TMath::Abs(p1.dz - p4.dz);
            DeltaDxy4 = TMath::Abs(p2.dz - p3.dz);
        }
	
    }
    else
    {
        mass1 = calc_invM(p1.momenta, p3.momenta);
        mass2 = calc_invM(p2.momenta, p4.momenta);
        DeltaDxy1 = TMath::Abs(p1.dz - p2.dz);
        DeltaDxy2 = TMath::Abs(p2.dz - p3.dz);
	
	mass3 = calc_invM(p1.momenta, p4.momenta);
        mass4 = calc_invM(p2.momenta, p3.momenta);
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
    }
    else
    {
        masses.push_back(mass1);
        masses.push_back(mass2);
	d2->Fill(mass3,DeltaDxy3);
        d2->Fill(mass4,DeltaDxy4);
	d3->Fill(mass1,DeltaDxy1);
        d3->Fill(mass2,DeltaDxy2);
    }
  
  d4->Fill(mass1,DeltaDxy1);
  d4->Fill(mass2,DeltaDxy2);
  d4->Fill(mass3,DeltaDxy3);
  d4->Fill(mass4,DeltaDxy4);
  
  d5->Fill((DeltaDxy1+DeltaDxy2)/2,(DeltaDxy3+DeltaDxy4)/2);
  
  return masses; 
}






//                     ----------------------------------------------------

void fourtwo_reduc()
{
    TH1F *h1 = new TH1F("h1","h1", 200,-2,2);
    TH1F *h2 = new TH1F("h2","h2", 200,-2,2);
    TH1F *h3 = new TH1F("h3","h3", 200,-2,2);
    
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
	    
	    // Check momentum conservation in histograms
	    
	    h1->Fill(px[0]+px[1]+px[2]+px[3]+px[4]+px[5]-6500*(ThxL+ThxR));
	    h2->Fill(py[0]+py[1]+py[2]+py[3]+py[4]+py[5]+6500*(ThyL+ThyR));
	    h3->Fill(pz[0]+pz[1]+pz[2]+pz[3]+pz[4]+pz[5]);
	    
	    
	    // Reduce track to 4-track	
	    std::vector<float> pairs = reduction(dxys,qs);
	    
	    // If no 4-track is found, then continue to next data point
	    if (pairs.size() == 0)
	    {
	        continue;
	    }

	    
	    // Check if things work as intended 
	    std::cout<< dxy[0] << " , " << dxy[1] << " , " << dxy[2] << " , " << dxy[3] << " , " << dxy[4] << " , " << dxy[5] << std::endl;
	    std::cout<< qs[0] << " , " << qs[1] << " , " << qs[2] << " , " << qs[3] << " , " << qs[4] << std::endl;
	    std::cout<< "Reduction;" << " " << pairs[0] << " , " << pairs[1] << " , " << pairs[2] << " , " << pairs[3] << std::endl;
	    
	    
            // Find indices to store data in particle-struct
	    
	    std::vector<int> indxs = findIndices(dxys,pairs);
	    
	    
	    //std::vector<float> curr_dxy = {dxy[indxs[0]],dxy[indxs[1]],dxy[indxs[2]],dxy[indxs[3]]};
	    //std::vector<float> rem_dxy;
	    
	    //for (int i = 0; i< track; i++)
	    //{
	    //    if(i != indxs[0] && i != indxs[1] && i != indxs[2] && i != indxs[3])
	    //  {
	    //      rem_dxy.push_back(dxy[i]);
	    //  }
	    //} 
	    
	    //if (isSplitThreeTrack(curr_dxy,rem_dxy) == true)
	    //{
	    //    continue;
	    //}
	    
	    if ((dxy[indxs[0]]+dxy[indxs[1]]+dxy[indxs[2]]+dxy[indxs[3]])/4 > 0.2)
	    {
	        continue;
	    }
	    
	    
	    
	    // Check momentum conservation in four-track
	    h11->Fill(px[indxs[0]]+px[indxs[1]]+px[indxs[2]]+px[indxs[3]]-6500*(ThxL+ThxR));
	    h12->Fill(py[indxs[0]]+py[indxs[1]]+py[indxs[2]]+py[indxs[3]]+6500*(ThyL+ThyR));
	    h13->Fill(pz[indxs[0]]+pz[indxs[1]]+pz[indxs[2]]+pz[indxs[3]]);
	    
	    std::cout<< indxs[0] << " , " << indxs[1] << " , " << indxs[2] << " , " << indxs[3] << std::endl;
	    
	    std::vector<float> p1 = {px[indxs[0]],py[indxs[0]],pz[indxs[0]]};
	    std::vector<float> p2 = {px[indxs[1]],py[indxs[1]],pz[indxs[1]]};
	    std::vector<float> p3 = {px[indxs[2]],py[indxs[2]],pz[indxs[2]]};
	    std::vector<float> p4 = {px[indxs[3]],py[indxs[3]],pz[indxs[3]]};
	    
	    Particle particle1(p1,qs[indxs[0]], TMath::Abs(dzs[indxs[0]]));
	    Particle particle2(p2,qs[indxs[1]], TMath::Abs(dzs[indxs[1]]));
	    Particle particle3(p3,qs[indxs[2]], TMath::Abs(dzs[indxs[2]]));
	    Particle particle4(p4,qs[indxs[3]], TMath::Abs(dzs[indxs[3]]));
	    
	    
	    // Finally pair up the 4-track into two 2-tracks
	    std::vector<float> mass = pair_up(particle1,particle2,particle3,particle4);

	    h8->Fill(mass[0]);
	    h8->Fill(mass[1]);
	    
	    // Cut on the 2-mass, so we only consider interesting pairs
	    if (TMath::Abs(mass[0]-mrho) < 0.15)
	    {
	        h6->Fill(mass[1]);
		if (TMath::Abs(mass[1]-mrho) < 0.15)
		{
		    h7->Fill(calc_fourmass(p1,p2,p3,p4));
		    
		}
		
		
	    }
	    
	    if (TMath::Abs(px[indxs[0]]+px[indxs[1]]+px[indxs[2]]+px[indxs[3]]-6500*(ThxL+ThxR))< 0.25 && TMath::Abs(py[indxs[0]]+py[indxs[1]]+py[indxs[2]]+py[indxs[3]]+6500*(ThyL+ThyR)) < 0.25)
	    {
	        if (TMath::Abs(mass[0]-mrho) < 0.15)
	        {
	            h6_2->Fill(mass[1]);
		    if (TMath::Abs(mass[1]-mrho) < 0.15)
		    {
		        h7_2->Fill(calc_fourmass(p1,p2,p3,p4));
		    
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
    legend6->AddEntry(g4,"#splitline{Rho peak:}{#splitline{ #mu = 0.753 #pm 0.009 GeV}{#sigma = 0.05 #pm 0.01 GeV}} ","l");
    legend6->AddEntry(g2,"Total fit: Chi2 / NDof : 56 / 41","l");
    legend6->SetTextSize(0.028);
    legend6->Draw();
    
    h6->SetTitle(" ; Inv. Mass [GeV]; ");
    
    TLatex s1;
    s1.SetTextSize(0.06);
    s1.DrawLatex(0.82*1.4,268, "#font[22]{CMS}");
    
    TLatex s2;
    s2.SetTextSize(0.034);
    s2.DrawLatex(0.87*1.4,301, "#sqrt{s} = 13 #font[22]{TeV}");
   
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
    
    TCanvas *c10 = new TCanvas("c10","c10",1600,600);
    d4->SetTitle(" ; Inv. Mass ; #Delta|dz|");
    d5->SetTitle(" ;#Delta|dz_{1}| ;#Delta|dz_{2}|");
    
    c10->Divide(2,1);
    c10->cd(1); d4->Draw("Colz");
    
    TLatex s15;
    s15.SetTextSize(0.06);
    s15.DrawLatex(0.3,3.02, "#font[22]{CMS}");
    
    TLatex s16;
    s16.SetTextSize(0.034);
    s16.DrawLatex(0.84*2,3.02, "#sqrt{s} = 13 #font[22]{TeV}");
    
    c10->cd(2); d5->Draw("Colz");
   
    
    TLatex s35;
    s35.SetTextSize(0.06);
    s35.DrawLatex(0.3,3.02, "#font[22]{CMS}");
    
    TLatex s36;
    s36.SetTextSize(0.034);
    s36.DrawLatex(0.84*3,3.02, "#sqrt{s} = 13 #font[22]{TeV}");
    
    
    
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
        
    
    
}



