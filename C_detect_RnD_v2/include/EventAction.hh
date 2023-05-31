#ifndef CANDYEventAction_h
#define CANDYEventAction_h 1

#include "G4UserEventAction.hh"


#include "G4AnalysisManager.hh"


#include "RunAction.hh"
#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4UnionSolid.hh"
/// Event action class
///




class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction* runAction);
  ~EventAction() override;
  
  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;
  
  void AddCount_Whole(G4int Count1) { fCount1 += Count1; }
  void AddCount_Positron(G4int Count2) { fCount2 += Count2; }
  void AddCount_notPositron(G4int Count3) { fCount3 += Count3; }
 
  void AddCount_SiPM(G4int N[12], G4int P){
    N[P] += 1;
  }
  void AddCount_PMT(G4int N[2], G4int P){
    N[P] += 1;
  }
  G4int f_SiPM_Count[12] = {0,};
  G4int f_PMT_Count[2] = {0,};

private:
  RunAction* fRunAction = nullptr;
  G4int   fCount1 = 0;
  G4int   fCount2 = 0;
  G4int   fCount3 = 0;
  G4int   fEventID=0;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
