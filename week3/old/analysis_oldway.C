#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/old/TOTEM40.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");


     //  Handy functions     ----------------------------------------------------
const Float_t mpi = 0.1396;
const Float_t mrho = 0.770;

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


Int_t count = 0;

     //                      ----------------------------------------------------


void analysis_oldway()
{

    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_dedx[200], trk_p[200];
    static Int_t trk_isK[3], trk_isPi[3], trk_isP[3];
    static Int_t trk_q[200];
    static Int_t ntrk;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q); 
    my_tree->SetBranchAddress("trk_isPi",&trk_isPi);
    my_tree->SetBranchAddress("trk_isP",&trk_isP);
    my_tree->SetBranchAddress("trk_isK",&trk_isK);
    my_tree->SetBranchAddress("trk_dedx",&trk_dedx);
    my_tree->SetBranchAddress("trk_p",&trk_p);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,0.2,2,200,0.2,2);
    
    TH1F *h1 = new TH1F("h1","px",200,-2,2);
    TH1F *h2 = new TH1F("h2","py",200,-2,2);
    TH1F *h3 = new TH1F("h3","pz",200,-4,4);
    
    TH2F *d2 = new TH2F("d2", "dedx", 200, 0,2,200,0,12);
 
    
    int numEnt = my_tree->GetEntries();

    for (int irow = 0; irow< numEnt; irow++)
    {
        my_tree->GetEntry(irow);
        if (ntrk == 4)
	{
	    Float_t px[4];
            Float_t py[4];
	    Float_t pz[4];
	    Float_t pt[4];
	    Float_t q[4];
	    Float_t dedx[4];
	    Float_t eta[4];
	    Float_t phi[4];
	    Float_t p[4];
	    
	    Int_t isPi[4];
	    Int_t isP[4];
	    Int_t isK[4];
	
	    for (int i = 0; i < 4; i++)
	    {   
	        eta[i] = trk_eta[i];
		phi[i] = trk_phi[i];
		pt[i] = trk_pt[i];
	        px[i] = pt[i]*TMath::Cos(phi[i]);
	        py[i] = pt[i]*TMath::Sin(phi[i]);
	        pz[i] = pt[i]*TMath::SinH(eta[i]);
		q[i] = trk_q[i];
		dedx[i] = trk_dedx[i];
		
		isPi[i] = trk_isPi[i];
	        isP[i] = trk_isP[i];
		isK[i] = trk_isK[i];
	    }
	    
	    h1->Fill(px[0]+px[1]+px[2]+px[3]);
	    h2->Fill(py[0]+py[1]+py[2]+py[3]);
	    h3->Fill(pz[0]+pz[1]+pz[2]+pz[3]);
		
	    
	    if (q[0]+q[1]+q[2]+q[3] != 0) 
	    {
	        continue;
	    }
	    
	    Float_t invM1, invM2,invM3, invM4;
	    
	    
	    if (q[0]+q[1] == 0)
	    {   
	        if (q[0] + q[2] == 0)
		{
		    invM1 = calc_InvM(px[0],py[0],pz[0],px[1],py[1],pz[1]);
	            invM2 = calc_InvM(px[2],py[2],pz[2],px[3],py[3],pz[3]);
		
		    invM3 = calc_InvM(px[0],py[0],pz[0],px[2],py[2],pz[2]);
	            invM4 = calc_InvM(px[1],py[1],pz[1],px[3],py[3],pz[3]);
		    
		    d1->Fill(invM1,invM2);
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
	   
	   
//	   if (isP[0] != 2 && isP[1] != 2 && isP[2] != 2 && isP[3] != 2)
//	   {
	      
//	       if ((isK[0] != 2 || dedx[0] < 4) && (isK[1] != 2 || dedx[1] < 4)&& (isK[2] != 2 || dedx[2] < 4) && (isK[3] != 2 || dedx[3] < 4))
//	       {
	       
//	           if (TMath::Abs(invM1-invM2) < 0.2 && dedx[0]> 0.5 && dedx[1] > 0.5 && dedx[2] > 0.5 && dedx[3] > 0.5)
//	           {
	              //d1->Fill(invM1,invM2);
//	           }  
	       
//	           if (TMath::Abs(invM3-invM4) < 0.2 && dedx[0]> 0.5 && dedx[1] > 0.5 && dedx[2] > 0.5 && dedx[3] > 0.5)
//	           {  
	              //d1->Fill(invM3,invM4);
//	           }
	       
//	       }
	       
//	   }
	   
	   if (isK[1] != 2 || dedx[1] < 4)
	   {
	       d2->Fill(TMath::Sqrt(px[1]*px[1]+py[1]*py[1]+pz[1]*pz[1]),dedx[1]);
	   }
		
	 }    


     }
      
      
    TCanvas *c1 = new TCanvas("c1","c1",1800,800);
    c1->Divide(3,1);
    d1->SetTitle("InvM1 / InvM2 ; InvM1 [GeV] ; InvM2 [GeV]");
    
    // gStyle->SetPalette("kCividis");
    c1->cd(1); h1->Draw();
    c1->cd(2); h2->Draw();
    c1->cd(3); h3->Draw();
    
    const char* name1 = "X-projection";
    const char* name2 = "Y-projection";
    
    auto Xhist = d1->ProjectionX(name1);
    auto Yhist = d1->ProjectionY(name2);
    Yhist->Sumw2();
    Xhist->Sumw2();
     
    Xhist->SetTitle("X-projection; Energy [GeV] ; #Counts");
    Yhist->SetTitle("Y-projection; Energy [GeV] ; #Counts");
    
    
    TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
    c2->Divide(2,2);
    
    c2->cd(1); d1->Draw("Colz");
    c2->cd(2); Xhist->Draw("E1");
    c2->cd(3); Yhist->Draw("E1");
    c2->cd(4); d2->Draw("Colz");
    Xhist->SetLineColor(kRed);
    h1->SetLineColor(kBlue);
    Yhist->SetLineColor(kRed);
    
    
    
     

}



