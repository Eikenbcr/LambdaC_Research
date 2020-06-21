R__LOAD_LIBRARY($ROOTSYS/test/libEvent.so)
void copytree()
{
   TString dir = "$ROOTSYS/test/Event.root";
   gSystem->ExpandPathName(dir);
   const auto filename = gSystem->AccessPathName(dir) ? "./Event.root" : "$ROOTSYS/test/Event.root";
   TFile oldfile(LcTopKKNew2.root);
   TTree *oldtree;
   oldfile.GetObject("T", DecayTree;3);
   // Deactivate all branches
   oldtree->SetBranchStatus("*", 0);
   // Activate only four of them
   for (auto activeBranchName : {"event", "fNtrack", "fNseg", "fH"})
      oldtree->SetBranchStatus(activeBranchName, 1);
   // Create a new file + a clone of old tree in new file
   TFile newfile("LcToPpKmKpMagDown.root", "recreate");
   auto newtree = oldtree->CloneTree();
   newtree->Print();
   newfile.Write();
}
