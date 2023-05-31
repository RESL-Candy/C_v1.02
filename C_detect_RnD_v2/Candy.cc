#include <iostream>
#include "G4RunManager.hh"
/////
#include "G4RandomTools.hh"
#include "Randomize.hh"
/////


////////////////////////////////
////////// 다룰 부분 ////////////
////////////////////////////////
#include "Construction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
////////////////////////////////
////////////////////////////////
////////////////////////////////



#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4Vismanager.hh"

#include "G4SteppingVerbose.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char** argv)// argc 에는 정보의 갯수, argv에는 입력값이 들어감
{
  long seed = 12345; // Replace with your desired seed value
  CLHEP::HepRandom::setTheSeed(seed);
  // 그냥 실행하면 argc가 1로 저장됨 
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }
  
  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();
    //  G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  

  
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////우리가 코드를 짜고 다루는 부분//////////////////////////////
  /////////////////////////////// 초기화 //////////////////////////////////////
  //
  // Detector construct
  // 어떤 물질로 어떤모양을 만들지
  MyDetectorConstruction* detectorconstruction = new MyDetectorConstruction();
  runManager->SetUserInitialization(detectorconstruction);

  // Physics List
  // 어떤 물리법칙이 적용될지
  G4VModularPhysicsList* physicsList = new MyPhysicsList();
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  // 다른 다양한 설정들 (초기입자, step, event, run) 등등
  ActionInitialization* actioninit = new ActionInitialization();
  runManager->SetUserInitialization(actioninit);
  
  
  // runManager->Initialize();
  
 
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////







  
  // Initialize visualization
  G4VisManager *visManager = new G4VisExecutive("Quiet");// 쓸때없는 warning message 안뜨게 
  visManager ->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();



  //UI command Line
  /*
  UImanager->ApplyCommand("/vis/open OGL");
  UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
  UImanager->ApplyCommand("/vis/drawVolume");
  UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
  UImanager->ApplyCommand("/vis/verbose warnings");
  UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
  UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge");
  UImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true");
  UImanager->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 3");
  UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
  //UImanager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set e+ yellow");
  //delete ui;
  */
  //UImanager->ApplyCommand(process/optical/cerenkov/setStackPhotons false); //cherencov 안하기 

  //매크로 만들고 나서 해야함  
  // Process macro or start UI session
  //
  
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");

    
    G4cout<<"/////////////////////////////////////////////////////////////"<<G4endl;
    G4cout<<"//                                                         //"<<G4endl;
    G4cout<<"//                                                         //"<<G4endl;
    G4cout<<"//                                                         //"<<G4endl;
    G4cout<<"//    ///////      /////////      ///////     //           //"<<G4endl;
    G4cout<<"//    //    //     //           //     ///    //           //"<<G4endl;
    G4cout<<"//    //     //    //           //       //   //           //"<<G4endl;
    G4cout<<"//    //     //    //           ///           //           //"<<G4endl;
    G4cout<<"//    //   ///     /////////     //           //           //"<<G4endl;
    G4cout<<"//    //////       //             ////////    //           //"<<G4endl;
    G4cout<<"//    //   //      //                    //   //           //"<<G4endl;
    G4cout<<"//    //    //     //                    //   //           //"<<G4endl;
    G4cout<<"//    //     //    //           //      ///   //           //"<<G4endl;
    G4cout<<"//    //      //   /////////     /////////    ///////////  //"<<G4endl;
    G4cout<<"//                                                         //"<<G4endl;
    G4cout<<"//                                                         //"<<G4endl;
    G4cout<<"//                                                         //"<<G4endl;
    G4cout<<"/////////////////////////////////////////////////////////////"<<G4endl;
    G4cout<<"                                                             "<<G4endl;
    G4cout<<"                            HELLO                            "<<G4endl;


    
    ui->SessionStart();

    ///이 밑에 뭐 입력하면 오류남 조심///
 
    
    
    delete ui;
  }
  

  
  delete visManager;
  delete runManager;//이거 안하면 예상치못한 종료 오류

  G4cout<<"/////////////////////////////////////////////////////////////"<<G4endl;
  G4cout<<"//                                                         //"<<G4endl;
  G4cout<<"//                                                         //"<<G4endl;
  G4cout<<"//                                                         //"<<G4endl;
  G4cout<<"//      //////     ///      ///   // //////   //     //    //"<<G4endl;
  G4cout<<"//    ///         ////     ////  // //    //  //    //     //"<<G4endl;
  G4cout<<"//    //         // //    // // // //     //  //  //       //"<<G4endl;
  G4cout<<"//    //        //////   //  //// //      //  ////         //"<<G4endl;
  G4cout<<"//    ///      //   //  //   /// //     //    //           //"<<G4endl;
  G4cout<<"//     ////// //    // //    // ////////     //            //"<<G4endl;
  G4cout<<"//                                                         //"<<G4endl;
  G4cout<<"//                                                         //"<<G4endl;
  G4cout<<"//                                                         //"<<G4endl;
  G4cout<<"/////////////////////////////////////////////////////////////"<<G4endl;
  G4cout<<"                                                             "<<G4endl;
  G4cout<<"                         Good bye :)                         "<<G4endl;

  
  return 0;
}
  
