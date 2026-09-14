std::vector<std::string> GetNuMIBeamlineFluxSYstKnobNames();

void GenerateNuMIBeamlineSystCov()
{

  std::vector<std::string> detNames = GetNuMIBeamlineFluxSYstKnobNames();

  const unsigned int NKnob = detNames.size();

  // Make the covariance matrix: ((0,1), (1,0))
  TFile *file = new TFile("gundaminput_numibeamlinesyst.root","RECREATE");
  file->cd();

  TMatrixTSym<double> numibeamlinesyst_cov(NKnob);
  for (int i = 0; i < NKnob; i++) {
    numibeamlinesyst_cov(i, i) = 1.0;
  }
  numibeamlinesyst_cov.Write("numibeamlinesyst_cov");

  TObjArray numibeamlinesyst_param_names;
  for(const auto& name: detNames){
    numibeamlinesyst_param_names.Add( new TObjString(name.c_str()) );
    std::cout << "@@ Writing " << name << std::endl;
  }
  file->WriteObjectAny( &numibeamlinesyst_param_names, "TObjArray", "numibeamlinesyst_param_names" );

  TVectorD numibeamlinesyst_param_prior(NKnob);
  TVectorD numibeamlinesyst_param_lb(NKnob);
  TVectorD numibeamlinesyst_param_ub(NKnob);
  for(int i=0; i<NKnob; i++){
    numibeamlinesyst_param_prior[i] = 0.;
    numibeamlinesyst_param_lb[i] = -3.;
    numibeamlinesyst_param_ub[i] = +3.;
  }
  numibeamlinesyst_param_prior.Write("numibeamlinesyst_param_prior");
  numibeamlinesyst_param_lb.Write("numibeamlinesyst_param_lb");
  numibeamlinesyst_param_ub.Write("numibeamlinesyst_param_ub");

  file->Close();
}

std::vector<std::string> GetNuMIBeamlineFluxSYstKnobNames(){

  return 
    {
    "hysyst_beam_horn_2kA",
    "hysyst_beam_horn1_x_3mm",
    "hysyst_beam_horn1_y_3mm",
    "hysyst_beam_spot_1_3mm",
    "hysyst_beam_spot_1_7mm",
    "hysyst_beam_horn2_x_3mm",
    "hysyst_hpc_0",
    "hysyst_hpc_1",
    "hysyst_hpc_2",
    "hysyst_hpc_3",
    "hysyst_hpc_4",
    "hysyst_hpc_5",
    "hysyst_hpc_6",
    "hysyst_hpc_7",
    "hysyst_hpc_8",
    "hysyst_hpc_9",
    "hysyst_hpc_10",
    "hysyst_hpc_11",
    "hysyst_hpc_12",
    "hysyst_hpc_13",
    "hysyst_hpc_14"
    };

}