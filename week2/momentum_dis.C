TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");
    
void momentum_dis()
{

    static Float_t trk_p,trk_pt, trk_eta, trk_phi;
    static Float_t trk_isPi, trk_isK, trk_isP;
    static Float_t ntrk;
    
    my_tree->SetBranchAddress("trk_p",&trk_p);
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    
    my_tree->SetBranchAddress("trk_isPi",&trk_isPi);
    my_tree->SetBranchAddress("trk_isK",&trk_isK);
    my_tree->SetBranchAddress("trk_isP",&trk_isP);
    
    my_tree->SetBranchAddress("ntrk",&ntrk);
    
    TH1F *h1 = new TH1F("h1","Includes P",200,-2,2);
    TH1F *h2 = new TH1F("h2","Includes K",200,-2,2);
    TH1F *h3 = new TH1F("h3","Includes Pi",200,-2,2);
    TH1F *h4 = new TH1F("h4","py",200,-2,2);
    
    TH2F *d1 = new TH2F("d1","py/isP",200,-4,4,4,0,4);
    TH2F *d2 = new TH2F("d1","py/isK",200,-4,4,4,0,4);
    TH2F *d3 = new TH2F("d1","py/isPI",200,-4,4,4,0,4);

    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow< numEnt; irow++)
    {
    my_tree->GetEntry(irow);
        if (ntrk >= 2)
	{
	    Float_t px = trk_pt*TMath::Cos(trk_phi);
	    Float_t py = trk_pt*TMath::Sin(trk_phi);
	    Float_t pz = trk_pt*TMath::SinH(trk_eta);
	
	    if (trk_isP >= 1)
	    {
	       h1->Fill(py);
	    }
	
	    if (trk_isK >= 1)
	    {
	       h2->Fill(py);
	    }
	
	    if (trk_isPi >= 1)
	    {
	        h3->Fill(py);
	     }
	
	    h4->Fill(py);
	
            if(irow < 10)		
	    {
                std::cout << py << std::endl;
	     }

	}
    }
      
//    h1->Sumw2();
//    Double_t scale = 1/h1->Integral();
//    h1->Scale(scale);
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->Divide(2,2);
    c1->cd(1);h1->Draw();
    c1->cd(2);h2->Draw();
    c1->cd(3);h3->Draw();
    c1->cd(4);h4->Draw();
    

}



