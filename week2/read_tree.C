TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");
    
void read_tree()
{

    static Float_t trk_p,trk_pt, trk_eta, trk_phi;
    
    my_tree->SetBranchAddress("trk_p",&trk_p);
    my_tree->SetBranchAddress("trk_pt",&trk_pt);
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    my_tree->SetBranchAddress("trk_phi",&trk_phi);
    
    TH1F *h1 = new TH1F("h1","p",200,0,5);
    TH1F *h2 = new TH1F("h2","pt",200,0,2);
    TH1F *h3 = new TH1F("h3","eta",200,-4,4);
    TH1F *h4 = new TH1F("h4","phi",200,-4,4);

    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow< numEnt; irow++)
    {
        my_tree->GetEntry(irow);
	 
        h1->Fill(trk_p);
	h2->Fill(trk_pt);
	h3->Fill(trk_eta);
	h4->Fill(trk_phi);
	
	
        if(irow < 10)		
	{
            std::cout << trk_p << std::endl;
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



