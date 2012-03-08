#ifndef _BUMP_H_
#define _BUMP_H_

class SectorMap;
class BeamLine;

class Bump{
 private:  
  double Q_x;

  double decr[4]; 
  double defl[4]; 
  SectorMap *kickers[4]; 
  double b_21;
  double b_I1;
  double b_I2;
  double b_31;
  double b_32;
  double b_41;
  double b_42;
  double b_43;
  double d_I1;
  double d_I2;   

  double  a0B1, a1B1, a2B1, bB1;  // S11MB1
  double  a0B2, a1B2, a2B2, bB2;  // S12MB2
  double  a0B3, a1B3, a2B3, bB3;  // S01MB3
  double  a0B4, a1B4, a2B4, bB4; // S03MB4

 public:  
  Bump() {} 
  Bump(double Q); 
  ~Bump() {}  
  void BumpSp(BeamLine* bl,unsigned &max_inj, int id, double amp0, double ampp0, double delAmp); 
  void BumpModi(BeamLine* bl,double amp0);
  void decrement(); 
  void decrement2(double amp0, double ampp0, int id); 
  void decrementModi(double amp0);  
  
};

#endif  // _BUMP_H_
