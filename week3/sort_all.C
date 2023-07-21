#include <tuple>

     //  Handy functions     ----------------------------------------------------
Float_t mkaon = 0.493;
Float_t mpi = 0.1396;
Float_t mrho = 0.770;

Float_t calc_E(Float_t px,Float_t py,Float_t pz)
{
  return TMath::Sqrt(mpi*mpi+px*px+py*py+pz*pz);

}

Float_t calc_InvM(Float_t px1,Float_t py1,Float_t pz1,Float_t px2,Float_t py2,Float_t pz2)
{
    Float_t E1 = calc_E(px1,py1,pz1);
    Float_t E2 = calc_E(px2,py2,pz2);
    
    return TMath::Sqrt(TMath::Power(E1+E2,2)-TMath::Power(px1+px2,2)-TMath::Power(py1+py2,2)-TMath::Power(pz1+pz2,2));
}

TH2F *d1 = new TH2F("d1","d1", 200, 0.2,1.2,200,0.2,1.2);
TH2F *d2 = new TH2F("d1","d1", 200, 0.2,1.2,200,0.2,1.2);

Int_t count1 = 0;
Int_t count2 = 0; 

Int_t is_rho(Float_t invM1, Float_t invM2)
{
   if (TMath::Abs(invM1 - (mrho-0.03)) < 0.1 && TMath::Abs(invM2 - (mrho-0.03)) < 0.1)
   {  
      d1->Fill(invM1,invM2);
      count1 += 1;
      return 1;
   }
   
   else
   {
      return 0;
   }
}

Int_t is_kaon(Float_t invM1, Float_t invM2)
{
   if (TMath::Abs(invM1 - (mkaon)) < 0.1 && TMath::Abs(invM2 - (mkaon)) < 0.1)
   {
      d2->Fill(invM1,invM2);
      count2 += 1;
      return 2;
   }
   
   else
   {
      return 0;
   }
}


     //                      ----------------------------------------------------

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TFile *file2 = new TFile("sorted_all.root","RECREATE");

void sort_all()
{ 

    TTree *my_tree = (TTree*)file1->Get("tree");
    TNtuple *tuple1 = new TNtuple("tuple1","tuple1","events");


    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_dxy[200];
    static Int_t trk_q[200];
    static Int_t ntrk;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow< numEnt; irow++)
    {
    my_tree->GetEntry(irow);
        if (ntrk == 4)
	{
	    Float_t px[5];
            Float_t py[5];
	    Float_t pz[5];
	    Int_t event;
	
	    for (int i = 0; i < 5; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		
	    }
	    
	    Float_t invM1, invM2,invM3, invM4; 
	    
	    if (trk_q[0]+trk_q[1] == 0)
	    {   
	        if (trk_q[0] + trk_q[2] == 0)
		{
		    invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	            invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		    invM3 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	            invM4 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		}
		
		else
		{
	            invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	            invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		    invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	            invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		}
		
	    }
	     
	    else
	    {
		invM1 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	        invM2 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		
		invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	        invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
	    }
	    
	    event = is_rho(invM1,invM2)+is_kaon(invM1,invM2);
	    tuple1->Fill(event); 
            event = is_rho(invM3,invM4)+is_kaon(invM3,invM4);
	    tuple1->Fill(event); 
	 
	 
	 }    

     }
      
    file2->Write();
    file2->Close();
    
    TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
    c1->Divide(2,2);
    c1->cd(1); d1->Draw("Colz");
    c1->cd(2); d2->Draw("Colz"); 
    std::cout << count1 << "  " << count2 << std::endl;
}



