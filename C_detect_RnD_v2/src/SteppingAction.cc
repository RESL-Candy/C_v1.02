#include "SteppingAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// collect energy deposited in this step
//G4double edepStep = step->GetTotalEnergyDeposit();
//aStep->GetTrack()->GetParentID() == 1 /양전자 조건
//aStep->GetTrack()->GetParentID() != 1 /양전자 아닌 조건


void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  if (!fSiPMScoringVolume && !fPMTScoringVolume) {
    const MyDetectorConstruction* detConstruction
      = static_cast<const MyDetectorConstruction*>
      (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    fSiPMScoringVolume = detConstruction->GetSiPMVolume();
    fPMTScoringVolume = detConstruction->GetPMTVolume();
  }
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  // aStep->GetTrack()->GetCurrentStepNumber() == 0 
  if (aStep->GetTrack()->GetCurrentStepNumber() == 1 && aStep->GetTrack()->GetTrackID() == 1)
    //1번 스텝에 첫출발 입자를 잡기 
    {
      // fill ntuple
      analysisManager->FillNtupleDColumn(1, aStep->GetPreStepPoint()->GetKineticEnergy());

      analysisManager->FillNtupleDColumn(2,  aStep->GetPreStepPoint()->GetPosition().x()/cm);
      analysisManager->FillNtupleDColumn(3,  aStep->GetPreStepPoint()->GetPosition().y()/cm);
      analysisManager->FillNtupleDColumn(4,  aStep->GetPreStepPoint()->GetPosition().z()/cm);

      analysisManager->FillNtupleDColumn(5, aStep->GetPreStepPoint()->GetMomentumDirection().x());
      analysisManager->FillNtupleDColumn(6, aStep->GetPreStepPoint()->GetMomentumDirection().y());
      analysisManager->FillNtupleDColumn(7, aStep->GetPreStepPoint()->GetMomentumDirection().z());

    }

  if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "opticalphoton"){
    if(aStep->GetTrack()->GetCurrentStepNumber() == 1){
      fEventAction->AddCount_Whole(1);
      if(aStep->GetTrack()->GetParentID() == 1){
        fEventAction->AddCount_Positron(1);
      }
      else{
        fEventAction->AddCount_notPositron(1);
      }
    }
  }
  // get volume of the current step
  G4LogicalVolume* volume
    = aStep->GetPreStepPoint()->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fSiPMScoringVolume && volume != fPMTScoringVolume ){
    return;
    } 




  G4int copyNo = 0;

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////SiPM scoring/////////////////////////////////////////
  //
  if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "opticalphoton"){
    if (volume == fSiPMScoringVolume){
      //if(aStep->GetTrack()->GetParentID() == 1){
        copyNo = aStep->GetPreStepPoint()->GetTouchableHandle()
        ->GetVolume()->GetCopyNo();  
        // fill ntuple
        fEventAction->AddCount_SiPM(fEventAction->f_SiPM_Count, copyNo);
        //analysisManager->FillNtupleDColumn(copyNo+12,  1);
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);        
      //}
    }
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////PMT scoring//////////////////////////////////////////
  //
    else if (volume == fPMTScoringVolume){
      copyNo = aStep->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetCopyNo();  
      fEventAction->AddCount_PMT(fEventAction->f_PMT_Count, copyNo);
      // fill ntuple
      //analysisManager->FillNtupleDColumn(copyNo+10,  1);
      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }
  
  
  
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
