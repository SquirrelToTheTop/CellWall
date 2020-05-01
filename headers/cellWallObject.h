#ifndef CELLWALLMONOLAYER_CLASS
#define CELLWALLMONOLAYER_CLASS

class CellWallMonolayer{

 public:
  CellWallMonolayer(int, int);  
  ~CellWallMonolayer();

 private:
  int nstrands = 0; // number of strands
  int npgstrand = 0; // number of peptidoglycans per strand

  int total_npg = 0; // total number of peptidoglycans

};

#endif