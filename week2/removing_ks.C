#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");


     //  Handy functions     ----------------------------------------------------
Float_t mpi = 0.1396;
Float_t mrho = 0.770;
Float_t mkaon = 0.493;

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

tuple<Float_t,Float_t> best_pair(Float_t invM1,Float_t invM2,Float_t invM3,Float_t invM4, Float_t m)
{

    if (TMath::Abs(invM1-m)+ TMath::Abs(invM2-m) < TMath::Abs(invM3-m) + TMath::Abs(invM4-m))
    {
        return make_tuple(invM1, invM2);
    }
    
    else
    {
        return make_tuple(invM3, invM4);
    }
    
}



     //                      ----------------------------------------------------


void removing_ks()
{

    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_dxy[200];
    static Int_t trk_q[200];
    static Int_t ntrk;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q);
    my_tree->SetBranchAddress("trk_dxy",&trk_dxy);

    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,0.25,1.2,200,0.25,1.2);
    
    TH2F *d2 = new TH2F("d2", "dxy /InvM1",200,-15,15,200,0.25,1.2);
    TH2F *d3 = new TH2F("d3", "dxy /InvM2",200,-15,15,200,0.25,1.2);
    
    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow< numEnt; irow++)
    {
    my_tree->GetEntry(irow);
        if (ntrk == 4)
	{
	    Float_t px[5];
            Float_t py[5];
	    Float_t pz[5];
	    Float_t dxy[5];
	
	    for (int i = 0; i < 5; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		dxy[i] = trk_dxy[i];
		
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
		
		auto [value1, value2] = best_pair(invM1,invM2,invM3,invM4,mrho);
		auto [kval1,kval2] = best_pair(invM1,invM2,invM3,invM4,mkaon);
		
		if (TMath::Abs(kval1 - mkaon) > 0.1 & TMath::Abs(kval2 - mkaon) > 0.1)
		{
		    d1->Fill(value1,value2);
		}
		
		
		}
		
		else
		{
	        invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	        invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	        invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		
		auto [value1, value2] = best_pair(invM1,invM2,invM3,invM4,mrho);
		auto [kval1,kval2] = best_pair(invM1,invM2,invM3,invM4,mkaon);
		
		if (TMath::Abs(kval1 - mkaon) > 0.1 & TMath::Abs(kval2 - mkaon) > 0.1)
		{
		    d1->Fill(value1,value2);
		}
		
		}
		
	    }
	     
	    else
	    {
		invM1 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	        invM2 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		
		invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	        invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		
		auto [value1, value2] = best_pair(invM1,invM2,invM3,invM4,mrho);
		auto [kval1,kval2] = best_pair(invM1,invM2,invM3,invM4,mkaon);
		
		if (TMath::Abs(kval1 - mkaon) > 0.1 && TMath::Abs(kval2 - mkaon) > 0.1)
		{
		    d1->Fill(value1,value2);
		}
		
	    }
	 
	 
	 }    


     }
      
      
    TCanvas *c1 = new TCanvas("c1","c1",600,600);

    d1->SetTitle("InvM1 / InvM2 ; InvM1 [GeV] ; InvM2 [GeV]");
    
    d1->Draw("Colz");
    
    const char* name1 = "X-projection";
    const char* name2 = "Y-projection";
    
    auto Xhist = d1->ProjectionX(name1);
    auto Yhist = d1->ProjectionY(name2);
    Yhist->Sumw2();
    Xhist->Sumw2();
     
    Xhist->SetTitle("X-projection; Energy [GeV] ; #Counts");
    Yhist->SetTitle("Y-projection; Energy [GeV] ; #Counts");
    
    
    TCanvas *c2 = new TCanvas("c2","c2",1600,600);
    c2->Divide(2,1);
    c2->cd(1); Xhist->Draw();
    c2->cd(2); Yhist->Draw();
    
  //  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  //  c3->Divide(2,1);
  //  c3->cd(1); d2->Draw("Colz");
  //  c3->cd(2); d3->Draw("Colz");
    
    // Save Ntuples
    
    TFile new_file("Projections_noK.root","RECREATE");

    TNtuple *tuple1 = new TNtuple("tuple1","tuple1","BinX:NX:errX:BinY:NY:errY");

    Float_t BinX, NX, errX, BinY, NY, errY;

    Float_t ent;
    
    if (Xhist->GetNbinsX() > Yhist->GetNbinsX())
    {
        ent = Xhist->GetNbinsX();
    } 
    
    else
    {
        ent = Yhist->GetNbinsX();
    }
    
    for (int i=0; i < ent ; ++i)
    {
    BinX = Xhist->GetBinCenter(i);
    NX = Xhist->GetBinContent(i);
    errX = Xhist->GetBinError(i);
    
    BinY = Yhist->GetBinCenter(i);
    NY = Yhist->GetBinContent(i);
    errY = Yhist->GetBinError(i);
    
    tuple1->Fill(BinX,NX,errX,BinY,NY,errY);

    
    }

    new_file.Write();
    new_file.Close();

    
    
     


}



