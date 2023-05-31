#include "RunAction.hh"






RunAction::RunAction()
{
  
  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("Candy");
  analysisManager->SetVerboseLevel(1);
  //analysisManager->SetNtupleMerging(true);
  // Note: merging ntuples is available only with Root output


  // Creating ntuple
  //
  //Creating Histo
  /*
  analysisManager->CreateH1("SiPMScoring","SiPMscoring", 12, 0.5, 12.5);
  analysisManager->CreateH1("PMTScoring","PMTscoring", 2, 0.5, 2.5);
  analysisManager->CreateH1("Whole Scoring","Whole scoring", 100001, -0.5, 100000.5);
  analysisManager->CreateH1("Positron Scoring","Positron scoring", 48001, -0.5, 48000.5);
  analysisManager->CreateH1("non P Scoring","non P scoring", 48001, -0.5, 48000.5);
  analysisManager->CreateH1("Positron Energy Scoring","Positron Energy scoring", 200, 0., 4.);
  analysisManager->CreateH1("Positron z_momentum Scoring","Positron z_momentum scoring", 200, -10., 10.);
  analysisManager->CreateH1("Positron z_position Scoring","Positron z_position scoring", 230, -11.5, 11.5);
  */

  
  analysisManager->CreateNtuple("CANDY_tree", "Event_by_Event");
  analysisManager->CreateNtupleDColumn("Event_ID");
  analysisManager->CreateNtupleDColumn("Positron_Energy");
  //Primary Postion
  analysisManager->CreateNtupleDColumn("Positron_Position_x");
  analysisManager->CreateNtupleDColumn("Positron_Position_y");
  analysisManager->CreateNtupleDColumn("Positron_Position_z");
  //Primary Momentum direction
  analysisManager->CreateNtupleDColumn("Positron_Momentum_Direction_x");
  analysisManager->CreateNtupleDColumn("Positron_Momentum_Direction_y");
  analysisManager->CreateNtupleDColumn("Positron_Momentum_Direction_z");
  //Photon scoring 
  analysisManager->CreateNtupleDColumn("Photon_number_from_positron");
  analysisManager->CreateNtupleDColumn("Photon_number_from_nonpositron");
  //PMT Scoring
  analysisManager->CreateNtupleDColumn("Left_PMT");
  analysisManager->CreateNtupleDColumn("Right_PMT");
  //SiPM Scoring
  analysisManager->CreateNtupleDColumn("1_SiPM");//0
  analysisManager->CreateNtupleDColumn("2_SiPM");//1
  analysisManager->CreateNtupleDColumn("3_SiPM");//2
  analysisManager->CreateNtupleDColumn("4_SiPM");//3
  analysisManager->CreateNtupleDColumn("5_SiPM");//4
  analysisManager->CreateNtupleDColumn("6_SiPM");//5
  analysisManager->CreateNtupleDColumn("7_SiPM");//6
  analysisManager->CreateNtupleDColumn("8_SiPM");//7
  analysisManager->CreateNtupleDColumn("9_SiPM");//8
  analysisManager->CreateNtupleDColumn("10_SiPM");//9
  analysisManager->CreateNtupleDColumn("11_SiPM");//10
  analysisManager->CreateNtupleDColumn("12_SiPM");//11

  analysisManager->FinishNtuple();




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* )
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "Candy.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
  // G4String fileName = "B4.xml";
  analysisManager->OpenFile(fileName);
  G4cout << "Open File " << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

 
  const MyDetectorConstruction* detConstruction
   = static_cast<const MyDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  
  // Run conditions

  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  
  const PrimaryGeneratorAction* generatorAction
   = static_cast<const PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  
  
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

   G4cout
       << G4endl
       << " The run consists of " << nofEvents << " "<< runCondition
       << G4endl
       << "------------------------------------------------------------"
       << G4endl;

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
  // G4cout<<"루트파일 만들어짐"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::printEventproc()
{

    fEvent += 1;
  if(fEvent % 200 == 0)
    {
      G4cout<<"처리 Events : "<<fEvent<<"개"<<G4endl;
    }
  
}
