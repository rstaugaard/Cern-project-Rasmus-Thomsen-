TFile *file1 = new TFile("/eos/cms/store/group/phys_smp/CMS_TOTEM/ntuples/data/TOTEM43.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");
    
void read_tree()
{

    static Float_t trk_q, trk_pt, trk_eta;
    
    my_tree->SetBranchAddress("trk_eta",&trk_eta);
    
    TH1F *h1 = new TH1F("h1","h1",200,-4,4);

    int numEnt = my_tree->GetEntries();

       
    for (int irow = 0; irow< numEnt; irow++)
    {
        my_tree->GetEntry(irow);
	 
        h1->Fill(trk_eta);
	
        if(irow < 10)		
	{
            std::cout << trk_eta << std::endl;
	}
    }
      
//    h1->Sumw2();
//    Double_t scale = 1/h1->Integral();
//    h1->Scale(scale);
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    h1->Draw();
    

}
