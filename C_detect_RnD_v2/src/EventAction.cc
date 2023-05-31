#include "EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fCount1 = 0;
  fCount2 = 0;
  fCount3 = 0;
  for(int i=0;i<12;i++){
    f_SiPM_Count[i] = {0,};
  }
  for(int i=0;i<2;i++){
    f_PMT_Count[i] = {0,};
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  fEventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(0, fEventID);

  analysisManager->FillNtupleDColumn(8, fCount2);
  analysisManager->FillNtupleDColumn(9, fCount3);
  
  for(G4int i=0; i<2; i++){
    analysisManager->FillNtupleDColumn(i+10, f_PMT_Count[i]);    
  }
  for(G4int i=0; i<12; i++){
    analysisManager->FillNtupleDColumn(i+12, f_SiPM_Count[i]);
  }

  analysisManager->AddNtupleRow();  
  fRunAction->printEventproc();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

