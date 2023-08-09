TFile *file1 = new TFile("/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/skim/data2/TOTEM2/log/Studies/perfect_ones_to_be_copied_back/with_protons/Results/Ntuples/summary_00.root","read");
TTree *my_tree = (TTree*)file1->Get("tree");


void dedx_analysis()
{
    Int_t tree_detector_id[3000], tree_chip_id[3000];
    Float_t tree_momentum[3000], tree_pt[3000], tree_deltaE[3000], tree_path_length[3000], tree_eta[3000];
    
    my_tree->SetBranchAddress("tree_detector_id",& tree_detector_id);
    my_tree->SetBranchAddress("tree_chip_id",&tree_chip_id);
    
    my_tree->SetBranchAddress("tree_momentum",&tree_momentum);
    my_tree->SetBranchAddress("tree_pt",&tree_pt); 
    my_tree->SetBranchAddress("tree_deltaE",&tree_deltaE);
    my_tree->SetBranchAddress("tree_path_length",&tree_path_length);
    my_tree->SetBranchAddress("tree_eta",&tree_eta);
    
    for (int irow = 0; irow< numEnt; irow++)
    { 
        my_tree->GetEntry(irow);
	detector_id = tree_detector_id[0
	pt = tree_pt[0];
	p = tree_momentum[0];
	deltaE = tree_deltaE[0];
	dx = tree_path_length[0];
	eta = tree_eta[0];
    
    
    } 


}
