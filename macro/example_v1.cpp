//
// Created by mikhail on 8/1/21.
//
void example_v1(){
  auto file = TFile::Open( "~/Correlations/ag158_elliptic_2021_07_27.root" );
  std::vector<std::string> ep_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED" };
  std::vector<std::string> res_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED", "Mf_RESCALED", "Mb_RESCALED" };
  std::vector<std::string> components{"x1x1", "y1y1"};

  auto file_out = TFile::Open( "example.root", "recreate" );
  file_out->cd();

  for( auto qa : ep_vectors ){
    Correlation u_qa( file, "SP/uQ", {"protons_RESCALED", qa}, {"x1x1", "y1y1"});
    auto res_v = Functions::VectorResolutions3S( file, "SP/QQ", qa, res_vectors, components );
    for( auto R1 : res_v ) {
      auto v1 = u_qa * 2 / R1;
      file_out->cd();
      R1.Save("R1."+R1.Title());
      v1.Save("v1."+R1.Title());
    }
  }
  file_out->Close();
  file->Close();
}