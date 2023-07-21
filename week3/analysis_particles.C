#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/old/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");


     //  Handy functions     ----------------------------------------------------
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

Int_t best_pair(Float_t invM1,Float_t invM2,Float_t invM3,Float_t invM4)
{
    if (TMath::Abs(invM1-invM2) < TMath::Abs(invM3-invM4))
    {
        return 1;
    }
    
    else
    {
        return 2;
    }
}

Int_t pick(Float_t invM1,Float_t invM2,Float_t invM3,Float_t invM4)
{
    if (TMath::Abs((invM1-invM2)/2) > 0.2 && TMath::Abs((invM3-invM4)/2) > 0.2)
    {
        return best_pair(invM1,invM2,invM3,invM4);
    }
    
    else
    {
       if (TMath::Abs(invM1-invM3) > TMath::Abs(invM2-invM4))
       {
          return 1;
       }
       
       else
       {
          return 2;
       }
   
    }
}

Int_t count = 0;

     //                      ----------------------------------------------------


void analysis_particles()
{

    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_dxy[200], trk_dz[200];
    static Int_t trk_q[200];
    static Int_t ntrk;
    static Float_t ThxR, ThxL, ThyR, ThyL;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q);
    my_tree->SetBranchAddress("trk_dxy",&trk_dxy);
    my_tree->SetBranchAddress("trk_dz",&trk_dz);    
    my_tree->SetBranchAddress("ThxR",&ThxR);
    my_tree->SetBranchAddress("ThxL",&ThxL);
    my_tree->SetBranchAddress("ThyR",&ThyR);
    my_tree->SetBranchAddress("ThyL",&ThyL);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,-1.2,1.2,200,0.2,1.2);
    TH2F *d2 = new TH2F("d2","InvM1 / InvM2",200,-1.2,1.2,200,0.2,1.2);
    TH1F *h1 = new TH1F("h1","h1",200,0.25,1.2);
    TH1F *h2 = new TH1F("h2","h2",200,0.25,1.2);
    
    TH2F *d3 = new TH2F("d3","InvM / pz",200,-3,3,200,0.25,1.2);
    TH2F *d4 = new TH2F("d4","InvM / pT",200,0,3,200,0.25,1.2);
    
    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow< numEnt; irow++)
    {
    my_tree->GetEntry(irow);
        if (ntrk == 4)
	{
	    Float_t px[5];
            Float_t py[5];
	    Float_t pz[5];
	    Float_t pt[5];
	    Float_t eta[5];
	    Float_t dxy[5];
	    Float_t q[5];
	    Float_t dz[5];
	    Float_t dZ[5];
	
	    for (int i = 0; i < 5; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		pt[i] = trk_pt[i];
		eta[i] = trk_eta[i];
		dxy[i] = trk_dxy[i];
		q[i] = trk_q[i];
		dz[i] = trk_dz[i];
		dZ[i] = TMath::Sqrt(TMath::Power(dxy[i],2)+TMath::Power(dz[i],2));
		
	    }
	    
	    if (q[0]+q[1]+q[2]+q[3] != 0) 
	    {
	        continue;
	    }
	    
	    Float_t invM1, invM2,invM3, invM4;
	    Int_t val;
	    
	    if (q[0]+q[1] == 0)
	    {   
	        if (q[0] + q[2] == 0)
		{
		    invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	            invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		    invM3 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	            invM4 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		    
		    d3->Fill(TMath::Abs(pz[0])+TMath::Abs(pz[1]),invM1);
		    d3->Fill(TMath::Abs(pz[2])+TMath::Abs(pz[3]),invM2);
		    d3->Fill(TMath::Abs(pz[0])+TMath::Abs(pz[2]),invM3);
		    d3->Fill(TMath::Abs(pz[1])+TMath::Abs(pz[3]),invM4);
		    
		    d4->Fill(pt[0]+pt[1],invM1);
		    d4->Fill(pt[2]+pt[3],invM2);
		    d4->Fill(pt[0]+pt[2],invM3);
		    d4->Fill(pt[1]+pt[3],invM4);
		    

		    
		}
		
		else
		{
	            invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	            invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		    invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	            invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		    
		    
		    d3->Fill(TMath::Abs(pz[0])+TMath::Abs(pz[1]),invM1);
		    d3->Fill(TMath::Abs(pz[2])+TMath::Abs(pz[3]),invM2);
		    d3->Fill(TMath::Abs(pz[0])+TMath::Abs(pz[2]),invM3);
		    d3->Fill(TMath::Abs(pz[1])+TMath::Abs(pz[3]),invM4);
		    
		    d4->Fill(pt[0]+pt[1],invM1);
		    d4->Fill(pt[2]+pt[3],invM2);
		    d4->Fill(pt[0]+pt[3],invM3);
		    d4->Fill(pt[1]+pt[2],invM4);
		}
		
	    }
	     
	    else
	    {
		invM1 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	        invM2 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		
		invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	        invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		
		d3->Fill(TMath::Abs(pz[0])+TMath::Abs(pz[1]),invM1);
	        d3->Fill(TMath::Abs(pz[2])+TMath::Abs(pz[3]),invM2);
	        d3->Fill(TMath::Abs(pz[0])+TMath::Abs(pz[2]),invM3);
	        d3->Fill(TMath::Abs(pz[1])+TMath::Abs(pz[3]),invM4);
		    
	        d4->Fill(pt[0]+pt[2],invM1);
	        d4->Fill(pt[1]+pt[3],invM2);
	        d4->Fill(pt[0]+pt[3],invM3);
	        d4->Fill(pt[1]+pt[2],invM4);
		
	    }
	    
	   d1->Fill(invM1,invM2);
	   d1->Fill(invM1,invM2);
		     
	    if (pt[0]+pt[1]+pt[2]+pt[3] > 0)
	    
	    {
	         if (TMath::Abs(invM1-invM2) < 0.2)
		 {   
		     d2->Fill(invM1,invM2);
		     h1->Fill(invM1);
		     h2->Fill(invM2);
		 }
	         
		 if (TMath::Abs(invM3-invM4) < 0.2)
		 {   
		     d2->Fill(invM3,invM4);
		     h1->Fill(invM3);
		     h2->Fill(invM4);
		 }
	
            }
	    
	    
	    else
	    {   
	        continue;
	    }
		
	 }    


     }
      
      
    TCanvas *c1 = new TCanvas("c1","c1",1600,800);
    c1->Divide(2,1);
    d1->SetTitle("InvM1 / InvM2 ; InvM1 [GeV] ; InvM2 [GeV]");
    
    // gStyle->SetPalette("kCividis");
    c1->cd(1); d1->Draw("Colz");
    c1->cd(2); d2->Draw("Colz");
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
    c2->cd(1); Xhist->Draw("SAME"); h1->Draw("SAMEE1");
    c2->cd(2); Yhist->Draw("SAME"); h2->Draw("SAMEE1");
    Xhist->SetLineColor(kRed);
    h1->SetLineColor(kBlue);
    Yhist->SetLineColor(kRed);
   
    TCanvas *c4 = new TCanvas("c4","c4",1600,600); 
    c4->Divide(2,1);
    c4->cd(1); d3->Draw("Colz");
    c4->cd(2); d4->Draw("Colz");
    
    

}



