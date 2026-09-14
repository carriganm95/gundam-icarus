std::vector<std::string> GetDetectorSystKnobNames();

void makeParameterRootFile()
{

  std::vector<std::string> detNames = GetDetectorSystKnobNames();

  const unsigned int NKnob = detNames.size();

  // Make the covariance matrix: ((0,1), (1,0))
  TFile *file = new TFile("gundaminput_detectorsyst.root","RECREATE");
  file->cd();

  TMatrixTSym<double> detectorsyst_cov(NKnob);
  for (int i = 0; i < NKnob; i++) {
    detectorsyst_cov(i, i) = 1.0;
  }
  detectorsyst_cov.Write("detectorsyst_cov");

  TObjArray detectorsyst_param_names;
  for(const auto& name: detNames){
    detectorsyst_param_names.Add( new TObjString(name.c_str()) );
    std::cout << "@@ Writing " << name << std::endl;
  }
  file->WriteObjectAny( &detectorsyst_param_names, "TObjArray", "detectorsyst_param_names" );

  TVectorD detectorsyst_param_prior(NKnob);
  TVectorD detectorsyst_param_lb(NKnob);
  TVectorD detectorsyst_param_ub(NKnob);
  for(int i=0; i<NKnob; i++){
    detectorsyst_param_prior[i] = 0.;
    detectorsyst_param_lb[i] = -3.;
    detectorsyst_param_ub[i] = +3.;
  }
  detectorsyst_param_prior.Write("detectorsyst_param_prior");
  detectorsyst_param_lb.Write("detectorsyst_param_lb");
  detectorsyst_param_ub.Write("detectorsyst_param_ub");

  file->Close();
}

std::vector<std::string> GetDetectorSystKnobNames(){

  return 
    {
      "induction_high",
      "induction_low",
      "coherent_high",
      "coherent_low",
      "intrinsic_high",
      "intrinsic_low",
      "lifetime_high",
      "lifetime_low",
      "scintillation",
      "wire_gap",
      "yz_uniformity",
      "recombination",
    };

}