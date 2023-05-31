#include "PhysicsList.hh"


  
MyPhysicsList::MyPhysicsList()
{
  new FTFP_BERT;
  ReplacePhysics(new G4EmStandardPhysics_option4());
  //RegisterPhysics (new G4EmStandardPhysics());
  RegisterPhysics (new G4OpticalPhysics());
}

MyPhysicsList::~MyPhysicsList()
{
}

