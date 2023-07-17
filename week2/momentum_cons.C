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

tuple<Float_t,Float_t> best_pair(Float_t invM1,Float_t invM2,Float_t invM3,Float_t invM4)
{
    if (TMath::Abs(invM1-mrho)+ TMath::Abs(invM2-mrho) < TMath::Abs(invM3-mrho) + TMath::Abs(invM4-mrho))
    {
        return make_tuple(invM1, invM2);
    }
    
    else
    {
        return make_tuple(invM3, invM4);
    }
}


     //                      ----------------------------------------------------


void momentum_cons()
{

    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_p[200];
    static Int_t trk_q[200];
    static Float_t ThyR, ThyL,ThxR, ThxL;
    static Int_t ntrk;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q);
    
    my_tree->SetBranchAddress("ThxR",&ThxR);
    my_tree->SetBranchAddress("ThxL",&ThxL);
    my_tree->SetBranchAddress("ThyR",&ThyR);
    my_tree->SetBranchAddress("ThyL",&ThyL);
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);

    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    TH2F *d1 = new TH2F("d1","ThxL / ThxR",200,-1,1,200,-1,1);
    TH2F *d2 = new TH2F("d2","ThyL / ThyR",200,-1,1,200,-1,1);
    TH2F *d3 = new TH2F("d3","Rx / px",200,-2,2,200,-3,3);
    TH2F *d4 = new TH2F("d4","Ry / py",200,-2,2,200,-3,3);
    
    TH1F *h1 = new TH1F("h1","dx", 200,-1,1);
    TH1F *h2 = new TH1F("h2","dy",200,-1,1);
    TH2F *d5 = new TH2F("d5","dx/dy",100,-0.1,0.1,100,-0.1,0.1);
    
    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow<numEnt; irow++)
    {
    my_tree->GetEntry(irow);
        if (ntrk == 4)
	{   
	    Float_t pt[5];
	    Float_t px[5];
            Float_t py[5];
	    Float_t pz[5];
	    Float_t q[5];
	
	    for (int i = 0; i < 4; i++)
	    {   
	        px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	        py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	        pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		pt[i] = trk_pt[i];
		q[i] = trk_q[i];
	    }
	    
	    d1->Fill(6500*ThxL,6500*ThyL);
	    d2->Fill(6500*ThxR,6500*ThyR);
	    d3->Fill(-6500*(ThxL+ThxR),px[0]+px[1]+px[2]+px[3]);
	    d4->Fill(6500*(ThyL+ThyR),py[0]+py[1]+py[2]+py[3]);
	    
	    h1->Fill(px[0]+px[1]+px[2]+px[3]-6500*(ThxL+ThxR));
	    h2->Fill(6500*(ThyL+ThyR)+py[0]+py[1]+py[2]+py[3]);
	    d5->Fill(-6500*(ThxL+ThxR)+px[0]+px[1]+px[2]+px[3],6500*(ThyL+ThyR)+py[0]+py[1]+py[2]+py[3]);
	    
	    
	    
	 }
     }
      
      
    TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
    c1->Divide(2,2);
    c1->cd(1); d1->Draw("Colz");
    c1->cd(2); d2->Draw("Colz");
    c1->cd(3); d3->Draw("Colz");
    c1->cd(4); d4->Draw("Colz");
    
    
    TCanvas *c2 = new TCanvas("c2","c2", 800,800);
    c2->Divide(2,2);
    c2->cd(1); h1->Draw();
    c2->cd(2); h2->Draw(); 
    c2->cd(3); d5->Draw("Colz");  

}



