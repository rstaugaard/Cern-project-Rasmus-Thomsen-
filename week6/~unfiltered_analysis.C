
TFile *file1 = new TFile("../data_files/reduced_data_dz.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

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



void unfiltered_analysis()
