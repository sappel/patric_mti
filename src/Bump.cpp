// This class defines a local orbit bump for the injection of the beam.
// 4 kickers and a marker called INJ for the injection must be contained
// in the madx twiss file.
// created by Stefan Paret (SP), May 2010    
// modify by Sabrina Appel (SA), Okt 2011

using namespace std;

#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <mpi.h>
#include <PhyConstants.h>

typedef vector<double> vektor;

#include "Bump.h"
#include "SectorMap.h"

Bump::Bump(double Q){	
	Q_x=Q;
}

void Bump::BumpSp(BeamLine* bl, unsigned &max_inj, int id, double amp0, double ampp0, double delAmp){  
  
  list<SectorMap>::iterator el = bl->get_first_element();
  const list<SectorMap>::iterator ending = bl->get_end_element();
  

  // find injection point, injection kickers and their lattice functions
  double beta_I, psi_I, alpha_I, psi[4], beta[4];
  while(el != ending){
    if(el->get_name() == "\"INJ\""){
      beta_I = el->get_betx();
      psi_I = el->get_mux();
      alpha_I = el->get_alpx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S01MB3\""){
      kickers[0] = &*el;
      psi[0] = el->get_mux();
      beta[0] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S03MB4\""){
      kickers[1] = &*el;
      psi[1] = el->get_mux();
      beta[1] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S11MB1\""){
      kickers[2] = &*el;
      psi[2] = (el->get_mux()-Q_x*2.*PI);
      beta[2] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S12MB2\""){
      kickers[3] = &*el;
      psi[3] = (el->get_mux()-Q_x*2*PI);
      beta[3] = el->get_betx();
      break;
    }
    ++el;
  }
  if(el == ending){
    cout<<"Error: Not all kickers for local orbit bump could be found or mark for injection point is missing.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  // calculate deflection angles for local orbit bump adjusted to incoming beam
  b_21 = sqrt(beta[3]*beta[2])*sin(psi[3]-psi[2]);
  b_I1 = sqrt(beta_I*beta[2])*sin(psi_I-psi[2]);
  b_I2 = sqrt(beta_I*beta[3])*sin(psi_I-psi[3]);
  b_31 = sqrt(beta[0]*beta[2])*sin(psi[0]-psi[2]);
  b_32 = sqrt(beta[0]*beta[3])*sin(psi[0]-psi[3]);
  b_41 = sqrt(beta[1]*beta[2])*sin(psi[1]-psi[2]);
  b_42 = sqrt(beta[1]*beta[3])*sin(psi[1]-psi[3]);
  b_43 = sqrt(beta[1]*beta[0])*sin(psi[1]-psi[0]);
  d_I1 = sqrt(beta[2]/beta_I)*cos(psi_I-psi[2]) - alpha_I/beta_I*b_I1;
  d_I2 = sqrt(beta[3]/beta_I)*cos(psi_I-psi[3]) - alpha_I/beta_I*b_I2;

  defl[2] =(d_I2*amp0 - b_I2*ampp0)/b_21;  // S11MB1
  defl[3] =(-d_I1*amp0 + b_I1*ampp0)/b_21;  // S12MB2
  defl[0] =(-b_41*defl[2] - b_42*defl[3])/b_43;  // S01MB3
  defl[1] =(b_31*defl[2] + b_32*defl[3])/b_43;  // S03MB4

  if(id == 0){
    cout<<"Injection parameters: " <<", amp0/m=" << amp0 << ", delAmp/m=" << delAmp <<
      ", ampp0= " << ampp0<< ", Q_x=" << Q_x << ", max_inj=" << max_inj << endl;
    cout << "beta_I/m=" << beta_I << ", alpha_I=" << alpha_I << endl;  
     cout << "Deflection angles/(mm mrad): ";
    for(short i=0; i<4; ++i)
      cout << defl[i]*1000. << ", ";
    cout << endl;
    for(short i=0; i<4; ++i){  // check against maximal deflection
      if(defl[i] > 0.0084)
	cout << "WARNING: Required deflection angle in kicker " << i
	     << " exceeds limit of 8.4 mrad.\n";
    }  
	for(short i=0; i<4; ++i){  // check against zero deflection SA
      if(defl[i] == 0.00)
	cout << "WARNING: Required deflection angle in kicker " << i
	     << " 0  mrad.\n";         
    }
    if(max_inj<1){
      cout<<"Error: Beam too large for injection. Execution aborted.\n";
      cout.flush();
      MPI_Abort(MPI_COMM_WORLD, 0);
    }
  }

  double ab = delAmp/amp0;  // relative decrement
  for(short i=0; i<4; ++i){
    decr[i] = ab*defl[i];  // set reduction per revolution
    kickers[i]->get_K(1) = defl[i];  // set deflection angle
  }
  

  // compensate first decrementation at beginning of loop in Main
  kickers[0]->get_K(1) += decr[0];
  kickers[1]->get_K(1) += decr[1];
}  

void Bump::BumpModi(BeamLine* bl,double amp0){    
    
    // constant coefficient of the polynomial for each bump magnet, soure D. Ondreka SA
    a0B1=0.0472005,  a1B1=0.00281663, a2B1=0., 		   bB1=0.0625;  // S11MB1
    a0B2=0.08084513, a1B2=-0.000875,  a2B2=0., 		   bB2=0.0625;  // S12MB2
    a0B3=-1.198325,  a1B3=0.7003125,  a2B3=-0.094026,  bB3=0.075;   // S01MB3
    a0B4=-0.815207,  a1B4=0.354425,   a2B4=-0.0358875, bB4=0.05875; // S03MB4    

    list<SectorMap>::iterator el = bl->get_first_element();
    const list<SectorMap>::iterator ending = bl->get_end_element();

	// find injection point, injection kickers and their lattice functions
  double beta_I, psi_I, alpha_I, psi[4], beta[4];
  while(el != ending){
    if(el->get_name() == "\"INJ\""){
      beta_I = el->get_betx();
      psi_I = el->get_mux();
      alpha_I = el->get_alpx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S01MB3\""){
      kickers[0] = &*el;
      psi[0] = el->get_mux();
      beta[0] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S03MB4\""){
      kickers[1] = &*el;
      psi[1] = el->get_mux();
      beta[1] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S11MB1\""){
      kickers[2] = &*el;
      psi[2] = (el->get_mux()-Q_x*2.*PI);
      beta[2] = el->get_betx();
      break;
    }
    ++el;
  }
  while(el != ending){
    if(el->get_name() == "\"S12MB2\""){
      kickers[3] = &*el;
      psi[3] = (el->get_mux()-Q_x*2*PI);
      beta[3] = el->get_betx();
      break;
    }
    ++el;
  }
  if(el == ending){
    cout<<"Error: Not all kickers for local orbit bump could be found or mark for injection point is missing.\n";
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
	
	if(Q_x<4.289){
	defl[2] =  (a0B1+(a1B1+a2B1*Q_x)*Q_x)*amp0; // S11MB1
    defl[3] =  (a0B2+(a1B2+a2B2*Q_x)*Q_x)*amp0; // S12MB2
    defl[0] =  (a0B3+(a1B3+a2B3*Q_x)*Q_x)*amp0; // S01MB3
    defl[1] =  (a0B4+(a1B4+a2B4*Q_x)*Q_x)*amp0; // S03MB4  
    }
    if(Q_x>=4.289){
	defl[2] =  bB1*amp0; // S11MB1
	defl[3] =  bB2*amp0; // S12MB2
	defl[0] =  bB3*amp0; // S01MB3
	defl[1] =  bB4*amp0; // S03MB4
    } 
	
	for(short i=0; i<4; ++i){
    kickers[i]->get_K(1) = defl[i];  // set deflection angle
    }


  // compensate first decrementation at beginning of loop in Main
  kickers[0]->get_K(1) += decr[0];
  kickers[1]->get_K(1) += decr[1];
}


void Bump::decrement(){
  // works only as long as the addresses of the elements pointed to by kickers[] remain valid 
  if(kickers[0]->get_K(1) != 0.)
    for(short i=0; i<4; ++i)
      if(kickers[i]->get_K(1) > decr[0])
	kickers[i]->get_K(1) -= decr[i];
      else
	kickers[i]->get_K(1) = 0.;
}   

void Bump::decrement2(double amp0, double ampp0, int id){ 
  //SA
  defl[2] = (d_I2*amp0 - b_I2*ampp0)/b_21;  // S11MB1
  defl[3] = (-d_I1*amp0 + b_I1*ampp0)/b_21;  // S12MB2
  defl[0] = (-b_41*defl[2] - b_42*defl[3])/b_43;  // S01MB3
  defl[1] = (b_31*defl[2] + b_32*defl[3])/b_43;  // S03MB4   

  
  for(short i=0; i<4; ++i){
	//decr[i] = ab*defl[i];  // set reduction per revolution
    kickers[i]->get_K(1) = defl[i];  // set deflection angle
  }
  
	
}   

void Bump::decrementModi(double amp0){ 
   //SA

   if(Q_x<4.289){
   defl[2] =  (a0B1+(a1B1+a2B1*Q_x)*Q_x)*amp0; // S11MB1
   defl[3] =  (a0B2+(a1B2+a2B2*Q_x)*Q_x)*amp0; // S12MB2
   defl[0] =  (a0B3+(a1B3+a2B3*Q_x)*Q_x)*amp0; // S01MB3
   defl[1] =  (a0B4+(a1B4+a2B4*Q_x)*Q_x)*amp0; // S03MB4  
   }
  if(Q_x>=4.289){
  	defl[2] =  bB1*amp0; // S11MB1
  	defl[3] =  bB2*amp0; // S12MB2
  	defl[0] =  bB3*amp0; // S01MB3
   	defl[1] =  bB4*amp0; // S03MB4
   } 
  
  for(short i=0; i<4; ++i){
	//decr[i] = ab*defl[i];  // set reduction per revolution
    kickers[i]->get_K(1) = defl[i];  // set deflection angle
  }
  
	
}
