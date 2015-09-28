#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <sstream>

#define long_to_long    0.200
#define inner           0.080
#define side            95.000
#define silicon_guard   1.036



double get_overlap(double short_to_short){

  double ladder_direction[16] = {- short_to_short/2 - 4 * side - 3*inner,
                          - short_to_short/2 - 3 * side - 3*inner,
                          - short_to_short/2 - 3 * side - 2*inner,
                          - short_to_short/2 - 2 * side - 2*inner,
                          - short_to_short/2 - 2 * side - 1*inner,
                          - short_to_short/2 - 1 * side - 1*inner,
                          - short_to_short/2 - 1 * side,
                          - short_to_short/2,
                          short_to_short/2,
                          short_to_short/2 + 1 * side,
                          short_to_short/2 + 1 * side + 1*inner,
                          short_to_short/2 + 2 * side + 1*inner,
                          short_to_short/2 + 2 * side + 2*inner,
                          short_to_short/2 + 3 * side + 2*inner,
                          short_to_short/2 + 3 * side + 3*inner,
                          short_to_short/2 + 4 * side + 3*inner};

  double orthogonal_direction[16] =   {- long_to_long/2 - 4 * side - 3*long_to_long,
                              - long_to_long/2 - 3 * side - 3*long_to_long,
                              - long_to_long/2 - 3 * side - 2*long_to_long,
                              - long_to_long/2 - 2 * side - 2*long_to_long,
                              - long_to_long/2 - 2 * side - 1*long_to_long,
                              - long_to_long/2 - 1 * side - 1*long_to_long,
                              - long_to_long/2 - 1 * side,
                              - long_to_long/2,
                              long_to_long/2,
                              long_to_long/2 + 1 * side,
                              long_to_long/2 + 1 * side + 1*long_to_long,
                              long_to_long/2 + 2 * side + 1*long_to_long,
                              long_to_long/2 + 2 * side + 2*long_to_long,
                              long_to_long/2 + 3 * side + 2*long_to_long,
                              long_to_long/2 + 3 * side + 3*long_to_long,
                              long_to_long/2 + 4 * side + 3*long_to_long};
  /*
  double llong           = 0.200  + 2 * silicon_guard;
  double sshort          = 0.080  + 2 * silicon_guard;
  double sside            = 95.000 - 2 * silicon_guard;

  double tungsten_coordinates[8] =   {- llong/2 - 4 * sside - 3*sshort,
                                      - llong/2 - 2 * sside - 2*llong,
                                      - llong/2 - 2 * sside - 1*sshort,
                                      - llong/2,
                                        llong/2,
                                        llong/2 + 2 * sside + 1*sshort,
                                        llong/2 + 2 * sside + 2*llong,
                                        llong/2 + 4 * sside + 3*sshort,
  };
  */
  //for(int i=0; i<8; i++){
  //  printf("tung: %f \n",tungsten_coordinates[i]);
  //}



  int overlap = 0;
  for(int i=0; i<40000; i++){
    double position = 0.01 * i + 0.001 - 400.;
    for(int ladder=0; ladder<4; ladder++){
      if(ladder_direction[ladder*2]  +silicon_guard > position) continue;
      if(ladder_direction[ladder*2+1]-silicon_guard < position) continue;
      if(orthogonal_direction[ladder*2]  +silicon_guard > position) continue;
      if(orthogonal_direction[ladder*2+1]-silicon_guard < position) continue;

      //if(tungsten_coordinates[((int)(ladder/2))*2]   > position) continue;
      //if(tungsten_coordinates[((int)(ladder/2))*2+1] < position) continue;

      overlap+=1;
      break;
    }
  }

  printf("overlap: %f mm\n",overlap*0.01);
  return overlap*0.01;

}

void plot(){
  TH1F* histo = new TH1F("overla_optimization","overla_optimization", 101, -10.1, 999.9 );

  for (int i=-1; i<101; i++){
    double short_to_short = i*0.01;
    result = get_overlap(short_to_short);
    double relative_coverage = result*result / ( (side - 2* silicon_guard)*(side - 2* silicon_guard)*16);
    histo->Fill(short_to_short * 1000., relative_coverage);
  }

  std::ostringstream histotitle;
  histotitle<<"long_to_long = "<<long_to_long*1000.<<" [#mum]";
  histo->SetTitle((histotitle.str()).c_str());
  histo->SetStats(false);
  histo->GetXaxis()->SetTitle("short_to_short [#mum]");
  histo.Draw();
}
