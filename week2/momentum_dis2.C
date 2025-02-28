#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
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
    if (TMath::Abs(invM1-mrho)+ TMath::Abs(invM2-mrho) < TMath::Abs(invM3-mrho) + TMath::Abs(invM4-mrho))
    {
        return 1;
    }
    
    else
    {
        return 2;
    }
}


     //                      ----------------------------------------------------


void momentum_dis2()
{

    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_dxy[200];
    static Int_t trk_q[200];
    static Int_t ntrk;
    static Float_t ThxR, ThxL, ThyR, ThyL;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q);
    my_tree->SetBranchAddress("trk_dxy",&trk_dxy);
    my_tree->SetBranchAddress("ThxR",&ThxR);
    my_tree->SetBranchAddress("ThxL",&ThxL);
    my_tree->SetBranchAddress("ThyR",&ThyR);
    my_tree->SetBranchAddress("ThyL",&ThyL);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    TH1F *h1 = new TH1F("h1","px",200,-2,2);
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-2,2);
    
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,0.25,1.2,200,0.25,1.2);
    
    TH2F *d2 = new TH2F("d2", "Pt /dxy",200,0,2.5,200,-3,3);
    TH2F *d3 = new TH2F("d3", "Pt /dxy",200,0,2.5,200,-3,3);
    TH2F *d4 = new TH2F("d4", "Pt /dxy",200,-4,4,200,-4,4);
    TH2F *d5 = new TH2F("d5", "Pt /dxy",200,0,4,200,0,4);  
  
    
    TPrincipal *PCA = new TPrincipal(2); 
    
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
	
	    for (int i = 0; i < 5; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		pt[i] = trk_pt[i];
		eta[i] = trk_eta[i];
		dxy[i] = trk_dxy[i];
		
		
		h1->Fill(px[i]);
	        h2->Fill(py[i]); 
		h3->Fill(pz[i]);
		
		
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
		
		    auto val = best_pair(invM1,invM2,invM3,invM4);
		
		    if (val == 1)
		    {
		        d1->Fill(invM1,invM2);
		        d2->Fill(pt[0]+pt[1],dxy[0]+dxy[1]);
		        d3->Fill(pt[2]+pt[3],dxy[2]+dxy[3]);
		        d4->Fill(dxy[0]+dxy[1],dxy[2]+dxy[3]);
		        d5->Fill(pt[0]+pt[1],pt[2]+pt[3]);
		    }
		
		    else
		    {
		        d1->Fill(invM3,invM4);
		        d2->Fill(pt[0]+pt[2],dxy[0]+dxy[2]);
		        d3->Fill(pt[1]+pt[3],dxy[1]+dxy[3]);
		        d4->Fill(dxy[0]+dxy[2],dxy[1]+dxy[3]);
		        d5->Fill(pt[0]+pt[2],pt[1]+pt[3]);
		     }
		    		
		}
		
		else
		{
	            invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	            invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		    invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	            invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		
		    auto val = best_pair(invM1,invM2,invM3,invM4);
		    if (val == 1)
		    {
		        d1->Fill(invM1,invM2);
		        d2->Fill(pt[0]+pt[1],dxy[0]+dxy[1]);
		        d3->Fill(pt[2]+pt[3],dxy[2]+dxy[3]);
		        d4->Fill(dxy[0]+dxy[1],dxy[2]+dxy[3]);
		        d5->Fill(pt[0]+pt[1],pt[2]+pt[3]);
		    }
		
	            else
		    {
		        d1->Fill(invM3,invM4);
		        d2->Fill(pt[0]+pt[3],dxy[0]+dxy[3]);
		        d3->Fill(pt[1]+pt[2],dxy[1]+dxy[2]);
		        d4->Fill(dxy[0]+dxy[3],dxy[1]+dxy[2]);
		        d5->Fill(pt[1]+pt[2],pt[0]+pt[3]);
		     }
	
		
		}
		
	    }
	     
	    else
	    {
		invM1 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	        invM2 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		
		invM3 = calc_InvM(px[0],py[0],pz[0],px[3],py[3],pz[3]);
	        invM4 = calc_InvM(px[1],py[1],pz[1],px[2],py[2],pz[2]);
		
		auto val = best_pair(invM1,invM2,invM3,invM4);		
		if (val == 1)
		{
		    d1->Fill(invM1,invM2);
		    d2->Fill(pt[0]+pt[2],dxy[0]+dxy[2]);
		    d3->Fill(pt[1]+pt[3],dxy[1]+dxy[3]);
		    d4->Fill(dxy[0]+dxy[2],dxy[1]+dxy[3]);
		    d5->Fill(pt[0]+pt[2],pt[1]+pt[3]);
		}
		
		else
		{
		    d1->Fill(invM3,invM4);
		    d2->Fill(pt[0]+pt[3],dxy[0]+dxy[3]);
		    d3->Fill(pt[1]+pt[2],dxy[1]+dxy[2]);
		    d4->Fill(dxy[0]+dxy[3],dxy[1]+dxy[2]);
		    d5->Fill(pt[0]+pt[3],pt[1]+pt[2]);
		}
		
		
	    }
	 
	 
	 }    


     }
      
      
    TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
    h1->SetTitle("px; px [GeV]; # Counts");
    h2->SetTitle("py ; py [GeV] ; # Counts");
    h3->SetTitle("pz ; pz [GeV] ; # Counts");
    d1->SetTitle("InvM1 / InvM2 ; InvM1 [GeV] ; InvM2 [GeV]");
    
    
    
    c1->Divide(2,2);
    c1->cd(1);h1->Draw();
    c1->cd(2);h2->Draw();
    c1->cd(3);h3->Draw();
    c1->cd(4);d1->Draw("Colz");
    
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
    
    TCanvas *c3 = new TCanvas("c3","c3",1600,1600);
    c3->Divide(2,2);
    c3->cd(1); d2->Draw("Colz");
    c3->cd(2); d3->Draw("Colz");
    c3->cd(3); d4->Draw("Colz");
    c3->cd(4); d5->Draw("Colz");
    
    Double_t cov1 = d2->GetCovariance();
    Double_t means1[2] = {d2->GetMean(1),d2->GetMean(2)};
    Double_t std1[2] = {d2->GetStdDev(1),d2->GetStdDev(2)};
    
    
    Double_t cov2 = d3->GetCovariance();
    Double_t means2[2] = {d3->GetMean(1),d3->GetMean(2)};
    Double_t std2[2] = {d3->GetStdDev(1),d3->GetStdDev(2)};
    
    // Save Ntuples
    
    TFile new_file("Projections.root","RECREATE");

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



