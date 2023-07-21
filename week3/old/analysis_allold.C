#include <tuple>

TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/old/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");

TFile *file2 = new TFile("sorted_allold.root");


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

Int_t best_pair(Float_t invM1,Float_t invM2,Float_t invM3,Float_t invM4, Float_t m)
{
    if (TMath::Abs(invM1-m)+ TMath::Abs(invM2-m) < TMath::Abs(invM3-m) + TMath::Abs(invM4-m))
    {
        return 1;
    }
    
    else
    {
        return 2;
    }
}



     //                      ----------------------------------------------------


void analysis_allold()
{   
    Int_t count0 = 0;
    Int_t count1 = 0;
    Int_t count2 = 0;
 
    TNtuple* my_tuple;
    file2->GetObject("tuple1", my_tuple);
    
    float event;
    float *row_content;

    static Float_t trk_pt[200], trk_eta[200], trk_phi[200], trk_dxy[200], trk_dz[200], trk_dedx[200];
    static Int_t trk_q[200];
    static Int_t ntrk;
    static Float_t zPV;
    
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    my_tree->SetBranchAddress("trk_q",&trk_q);
    my_tree->SetBranchAddress("trk_dxy",&trk_dxy);
    my_tree->SetBranchAddress("zPV",&zPV);
    my_tree->SetBranchAddress("trk_dz",&trk_dz);
    my_tree->SetBranchAddress("trk_dedx",&trk_dedx);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    
    TH2F *d1 = new TH2F("d1","InvM1 / InvM2",200,0.25,1.2,200,0.25,1.2);
    TH2F *d2 = new TH2F("d2", "Pt / InvMass",200,0.4,2.8,200,-0.1,1.8);
    TH2F *d3 = new TH2F("d3", "dxy /InvMass",200,-4,4,200,-0.1,1.8);
    TH2F *d4 = new TH2F("d4", "zPV /InvMass",200,-24,24,200,0.2,1.2);
    TH2F *d5 = new TH2F("d5", "dz /InvMass",200,-4,4,200,0.2,1.2);  
  

    
    int numEnt = my_tree->GetEntries();
    Int_t j = -1;
       
    for (int irow = 0; irow< numEnt; irow++)
    {
        my_tree->GetEntry(irow);
	Int_t events[4];
	
        if (ntrk == 4)
	{   
	    j += 1;
	    my_tuple->GetEntry(j);
            row_content = my_tuple->GetArgs();
            events[0] = row_content[0];
	    
	    j += 1;
	    my_tuple->GetEntry(j);
            row_content = my_tuple->GetArgs();
            events[1] = row_content[0];
	    
	    
	    
	    if (events[0] > 0 || events[1] > 0)
	    {  	       
	       Float_t px[5];
               Float_t py[5];
	       Float_t pz[5];
	       Float_t pt[5];
	       Float_t eta[5];
	       Float_t dxy[5];
	       Float_t dz[5];
	       Float_t dedx[5];
	
	       for (int i = 0; i < 5; i++)
	       {   
	           px[i] = trk_pt[i]*TMath::Cos(trk_phi[i]);
	           py[i] = trk_pt[i]*TMath::Sin(trk_phi[i]);
	           pz[i] = trk_pt[i]*TMath::SinH(trk_eta[i]);
		   pt[i] = trk_pt[i];
		   eta[i] = trk_eta[i];
		   dxy[i] = trk_dxy[i];
		   dz[i] = trk_dz[i];
		   dedx[i] = trk_dedx[i];
		   
	        }
	    
	       Float_t invM1, invM2, invM3, invM4; 
	    
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
		  
		  if (events[0] > 0)
		  {
		      d1->Fill(invM1,invM2);
	              d2->Fill(pt[0]+pt[1]+pt[2]+pt[3], invM1);
		      d2->Fill(pt[0]+pt[1]+pt[2]+pt[3],invM2);
		      d3->Fill(dxy[0]+dxy[1]+dxy[2]+dxy[3],invM1);
                      d3->Fill(dxy[0]+dxy[1]+dxy[2]+dxy[3],invM2);
		      d4->Fill(zPV,invM1);
	              d4->Fill(zPV,invM2);
	              d5->Fill(dz[0]+dz[1]+dz[2]+dz[3],invM1);
	              d5->Fill(dz[0]+dz[1]+dz[2]+dz[3],invM2);
		  }
		  
		  if (events[1] > 0)
		  {
		      d1->Fill(invM3,invM4);
	              d2->Fill(pt[0]+pt[1]+pt[2]+pt[3], invM3);
		      d2->Fill(pt[0]+pt[1]+pt[2]+pt[3],invM4);
		      d3->Fill(dxy[0]+dxy[1]+dxy[2]+dxy[3],invM3);
                      d3->Fill(dxy[0]+dxy[1]+dxy[2]+dxy[3],invM4);
		      d4->Fill(zPV,invM3);
	              d4->Fill(zPV,invM4);
	              d5->Fill(dz[0]+dz[1]+dz[2]+dz[3],invM3);
	              d5->Fill(dz[0]+dz[1]+dz[2]+dz[3],invM4);
		  }
		  
		  
	     }
	 }    


     }
      
      
    TCanvas *c1 = new TCanvas("c1","Sample space",700,700);
    d1->SetTitle("InvM1 / InvM2 ; InvM1 [GeV] ; InvM2 [GeV]");
    
    c1->cd(4);d1->Draw("Colz");
    
    const char* name1 = "X-projection";
    const char* name2 = "Y-projection";
    
    auto Xhist = d1->ProjectionX(name1);
    auto Yhist = d1->ProjectionY(name2);
    Yhist->Sumw2();
    Xhist->Sumw2();
     
    Xhist->SetTitle("X-projection; Energy [GeV] ; #Counts");
    Yhist->SetTitle("Y-projection; Energy [GeV] ; #Counts");
    
    
    TCanvas *c2 = new TCanvas("c2","Projections",1600,600);
    c2->Divide(2,1);
    c2->cd(1); Xhist->Draw();
    c2->cd(2); Yhist->Draw();
    
    TCanvas *c3 = new TCanvas("c3","Correlations with invariant mass",1600,1600);
    c3->Divide(2,2);
    c3->cd(1); d2->Draw("Colz");
    d2->SetTitle("InvM / pT ; pT [GeV] ; InvM [GeV]");
    c3->cd(2); d3->Draw("Colz");
    d3->SetTitle("InvM / dxy ; dxy [] ; InvM [GeV]");
    c3->cd(3); d4->Draw("Colz");
    d4->SetTitle("InvM / zPV ; dis [cm] ; InvM [GeV]");
    c3->cd(4); d5->Draw("Colz");
    d5->SetTitle("InvM / dz ; dz [] ; InvM [GeV]");
   
}



