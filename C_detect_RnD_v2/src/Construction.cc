#include "Construction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorConstruction::MyDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorConstruction::~MyDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* MyDetectorConstruction::Construct()
{
  //////////////////   materials   //////////////////////////
  //
  //
  /////complete material////
  //
  G4NistManager *nist = G4NistManager::Instance();
  G4Material* Mat_Air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* Mat_Teflon = nist->FindOrBuildMaterial("G4_TEFLON");
  //
 

  
 


  
  /////Compound material (need to compound)///
  //
  // Define the elements that make up LAB
  G4Element* El_C = nist->FindOrBuildElement("C");
  G4Element* El_H = nist->FindOrBuildElement("H");

  ////Define the chemical formula for LAB and densities of its constituents///
  //
  fMat_LAB = new G4Material("LAB", 0.853*g/cm3, 2);
  fMat_LAB->AddElement(El_C, 18);
  fMat_LAB->AddElement(El_H, 30);
  //
  ///////////////Define LAB's scintillation properties/////////// 
  //
  std::vector<G4double> lab_Energy = { 2.67 * eV, 2.92 * eV, 3.26 * eV };
  std::vector<G4double> lab_SCINT = { 0.25, 0.5 , 1.0 };//emission rate
  std::vector<G4double> lab_RIND  = { 1.52, 1.53, 1.55 };
  std::vector<G4double> lab_ABSL  = { 1000. * cm, 1000. * cm, 1000. * cm };
  fmat_Prop_table_LAB = new G4MaterialPropertiesTable();
  fmat_Prop_table_LAB->AddProperty("SCINTILLATIONCOMPONENT1", lab_Energy, lab_SCINT);//emission rate 
  // fmat_Prop_table_LAB->AddProperty("SCINTILLATIONCOMPONENT2", lab_Energy, lab_SCINT);
  fmat_Prop_table_LAB->AddProperty("RINDEX", lab_Energy, lab_RIND);
  fmat_Prop_table_LAB->AddProperty("ABSLENGTH", lab_Energy, lab_ABSL);
  fmat_Prop_table_LAB->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  fmat_Prop_table_LAB->AddConstProperty("RESOLUTIONSCALE", 1.0);
  fmat_Prop_table_LAB->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
  //fmat_Prop_table_LAB->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  //fmat_Prop_table_LAB->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  // fmat_Prop_table_LAB->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
  fMat_LAB->SetMaterialPropertiesTable(fmat_Prop_table_LAB);
  // Set the Birks Constant for the LXe scintillator
  fMat_LAB->GetIonisation()->SetBirksConstant(0.151 * mm / MeV);

  ///////////////Define Teflon's scintillation properties///////////
  //
  //
  /*
  fmat_Prop_table_Teflon = new G4MaterialPropertiesTable();

  std::vector<G4double> teflon_Energy = { 2.67 * eV, 2.92 * eV, 3.26 * eV };
  std::vector<G4double> teflon_RIND  = { 2.00, 2.00, 2.00 };
  std::vector<G4double> teflon_ABSL  = { 1000. * cm, 1000. * cm, 1000. * cm };
  fmat_Prop_table_Teflon->AddProperty("RINDEX", teflon_Energy, teflon_RIND);
  fmat_Prop_table_Teflon->AddProperty("ABSLENGTH", teflon_Energy, teflon_ABSL);  
  Mat_Teflon->SetMaterialPropertiesTable(fmat_Prop_table_Teflon);
  */
  
  
  
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //////////////////  Detector Define   //////////////////////////
  //
  //

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  
  ///////////// World /////////////////
  //
  //
  G4double world_Size = 3. * m;
  
  G4Box* solide_World = new G4Box("Worldsol",                            //its name
			      .5 * world_Size,                        //its size
			      .5 * world_Size,
			      .5 * world_Size);
  
  G4LogicalVolume* logical_World = new G4LogicalVolume(solide_World,  //its solid
						       Mat_Air,       //its material
						       "Worldlog");      //its name
  
  G4VPhysicalVolume* Phys_World = new G4PVPlacement(0,                //no rotation
						    G4ThreeVector(),  //at (0,0,0)
						    logical_World,    //its logical volume
						    "Worldphy",       //its name
  						    0,                //its mother  volume
						    false,            //no boolean operation
						    0,                //copy number
						    checkOverlaps);   //overlaps checking

  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  ////// Teflon cylinder /////
  G4double innerRadius = 0.;
  G4double outerRadius = 0.055 * m;
  G4double height = 0.45 * m;
  
  G4Tubs* teflonSol = new G4Tubs("Teflonsol",
				 innerRadius,
				 outerRadius,
				 .5 * height,
				 0.,         //tub 에서 시작 각도 
				 2. * M_PI);
  
  G4LogicalVolume* teflonLog = new G4LogicalVolume(teflonSol,
						   Mat_Teflon,
						   "Teflonlog");
  
  G4ThreeVector teflonPos = G4ThreeVector(0., 0., 0.); // position of Teflon cylinder
  G4VPhysicalVolume* teflonPhy = new G4PVPlacement(nullptr,
						   teflonPos,
						   teflonLog,
						   "Teflonphy",
						   logical_World, 
						   false,
						   0,
						   checkOverlaps);

  
  

  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  ////// LAB LS /////

   
  G4double LS_innerRadius = 0.;
  G4double LS_outerRadius = 0.05 * m;
  G4double LS_height = 0.1 * m;
  
  G4Tubs* LS_Sol = new G4Tubs("LAB_Sol",
			      LS_innerRadius,
			      LS_outerRadius,
			      .5 *LS_height,
			      0.,
			      2. * M_PI);

  
  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  ////// LAB LS GUIDE/////  
  G4double Guide_Axial = 0.;
  G4double Guide_Down_R = 25.*mm;
  G4double Guide_Up_R = 50.*mm;
  G4double Guide_Height = 50*mm;

  G4Cons * Guide_Sol = new G4Cons("LAB_Guide",
				  Guide_Axial,
				  Guide_Up_R,
				  Guide_Axial,
				  Guide_Down_R,
				  Guide_Height,
				  0.,
				  2. * M_PI);

  
  G4UnionSolid* Union_LS_Sol = new G4UnionSolid("Union_LS_Sol",
                  LS_Sol,
                  Guide_Sol,
                  0,
                  G4ThreeVector(0,0,(LS_height+Guide_Height*2.)*0.5));
                  
  G4LogicalVolume* Union_LS_Log = new G4LogicalVolume(Union_LS_Sol,
						     fMat_LAB,
						     "Union_LS_Log");
  
  G4double thetaX = 180*deg;
  G4double thetaY = 0*deg;
  G4double thetaZ = 0*deg;
  G4RotationMatrix * Rotation_Matrix_Union = new G4RotationMatrix();
  Rotation_Matrix_Union->rotateX(thetaX);

  G4VPhysicalVolume* Union_Plus_LS_Phys = new G4PVPlacement(0,
		    G4ThreeVector(0., 0., 0.5*(LS_height)),
		    Union_LS_Log,
		    "Union_Plus_LS_Phys",
		    teflonLog,  
		    false,
		    0,
		    checkOverlaps);

  G4VPhysicalVolume* Union_Minus_LS_Phys = new G4PVPlacement(Rotation_Matrix_Union,
		    G4ThreeVector(0., 0., -0.5*(LS_height)),
		    Union_LS_Log,
		    "Union_Minus_LS_Phys",
		    teflonLog,  
		    false,
		    1,
		    checkOverlaps);     
  
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// PMT ////////////////////////////
  //

  G4double PMT_height = 0.001 * m;

  G4double PMT_pos_z=0.;
  
  G4Tubs* PMT_Sol = new G4Tubs("PMT_Sol",
			       0.,
			       Guide_Down_R,
			       .5 *PMT_height,
			       0.,
			       2. * M_PI);

  fPMT_SD_LOG = new G4LogicalVolume(PMT_Sol,
				    fMat_LAB,
				    "PMT_Log");

  // position of PMT
  

  for(G4int i=0; i<2; i++)
    {
      PMT_pos_z=(Guide_Height*4.+LS_height*2+PMT_height)*0.5*(2*i-1);     
      //G4VPhysicalVolume* PMT_Phy =
      new G4PVPlacement(nullptr,
			G4ThreeVector(0., 0., PMT_pos_z),
			fPMT_SD_LOG,
			"PMT_phy",
			teflonLog,  
			false,
			i,
			checkOverlaps);
      
    }

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////// SD ////////////////////////////
  //

  G4double SD_width = 0.006 * m;
  G4double SD_height = 0.001 * mm;
  G4int SiPM_N=6;
  G4double sipm_pos=0;
  
  G4Box* SiPM_sol = new G4Box("SiPM_sol",.5*SD_width,.5*SD_height,.5*SD_width);
  
  
  fSipm_SD_LOG = new G4LogicalVolume(SiPM_sol,fMat_LAB,"SiPM_LOG");

  
  for(G4int i=0; i<SiPM_N; i++)
    {
      sipm_pos = (LS_height)/2.*((2.*i+1.)/SiPM_N-1.);
      //G4VPhysicalVolume* SiPM_Phy =
      new G4PVPlacement(0,
			G4ThreeVector(0., 0., sipm_pos),
			fSipm_SD_LOG,
			"SiPM_PHY",
			Union_LS_Log,
			false,
			i,
			checkOverlaps);
    }
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////광섬유///////////////////////////////////////////////

  G4double Optical_Fiber_Length = 48*mm;
  G4double Optical_Fiber_R = 0.5*mm;
  G4int Optical_Fiber_N=6;
  G4double Optical_Fiber_Pos=0;


  G4Tubs* Optical_Fiber_Sol = new G4Tubs("Optical_Fiber_Sol",
  0.,
  Optical_Fiber_R,
  .5*Optical_Fiber_Length,
  0,
  360.*deg);
  
  fOptical_Fiber_LOG = new G4LogicalVolume(Optical_Fiber_Sol,fMat_LAB,"Optical_Fiber_Log");
  G4RotationMatrix* Rotation_Matrix_Optical_Fiber = new G4RotationMatrix();
  Rotation_Matrix_Optical_Fiber -> rotateX(90.0*deg);

  for(G4int i=0; i<Optical_Fiber_N; i++)
    {
      Optical_Fiber_Pos = (LS_height)/2.*((2.*i+1.)/Optical_Fiber_N-1.);
      new G4PVPlacement(Rotation_Matrix_Optical_Fiber,
			G4ThreeVector(0.,(.5*LS_outerRadius-.5*SD_height),Optical_Fiber_Pos),
			fOptical_Fiber_LOG,
			"Optical_Fiber_Phys",
			Union_LS_Log,
			false,
			i,
			checkOverlaps);

      new G4PVPlacement(Rotation_Matrix_Optical_Fiber,
			G4ThreeVector(0.,-(.5*LS_outerRadius-.5*SD_height),Optical_Fiber_Pos),
			fOptical_Fiber_LOG,
			"Optical_Fiber_Phys",
			Union_LS_Log,
			false,
			i,
			checkOverlaps);      
     
    }



  

  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////Surface define////////////////////////////////////
  //
  //
  //
  //
  //
  //////Define Teflon surface/////////////
  //

  /*   //Geant4 Application manual의 Listing 73에 dielectric_dielectric의 경우가 적혀있음
  //또한 teflon의 material property 또한 지정해줘야함.
  //그래서 LUT모델을 씀.
  std::vector<G4double> ephoton = { 2.67 * eV, 2.92 * eV, 3.26 * eV }; //입사에너지
  std::vector<G4double> reflectivity = { 1.0, 1.0, 1.0 }; //반사비율
  std::vector<G4double> efficiency = { 0.0, 0.0, 0.0 }; //투과율
  G4MaterialPropertiesTable* Teflon_Surface = new G4MaterialPropertiesTable();
  Teflon_Surface->AddProperty("REFLECTIVITY", ephoton, reflectivity);
  Teflon_Surface->AddProperty("EFFICIENCY", ephoton, efficiency);
  */




  
  G4OpticalSurface* Tef_roughSurf = new G4OpticalSurface("tef_roughSurf");
  Tef_roughSurf->SetType(dielectric_LUTDAVIS); //이거 하면 LUT를 사용하여 데이터를 갖고와서 씀. 
  Tef_roughSurf->SetModel(DAVIS);// LUT DAVIS 모델은 DAVIS 고정
  Tef_roughSurf->SetFinish(RoughTeflon_LUT);//이거 다른 세팅이면 반사가 안됨. --> 해결!  

  new G4LogicalBorderSurface("tef_roughSurf",
			     Union_Plus_LS_Phys, teflonPhy,
			     Tef_roughSurf);
  
  new G4LogicalBorderSurface("tef_roughSurf",
			     Union_Minus_LS_Phys, teflonPhy,
			     Tef_roughSurf); 

 
  //Tef_roughSurf->SetSigmaAlpha(0.1);//0~1 : grounded ~ polished
  //Tef_roughSurf->SetMaterialPropertiesTable(Teflon_Surface);



  
  //G4LogicalSkinSurface* Tef_skinSurface = new G4LogicalSkinSurface("TargetSurface",
  //								   teflonLog,Tef_roughSurf);

  //위의 skin과 border의 차이는 단순 표면 그리고 경계의 차이이다.










  
  logical_World->SetVisAttributes(G4VisAttributes::GetInvisible());//world 안보이게 함 
  
  return Phys_World;
  
}
/*
void MyDetectorConstruction::ConstructSD()
{
  MySensitiveDetector* SD = new MySensitiveDetector("SiPM_SD");
  fSipm_SD_LOG->SetSensitiveDetector(SD);
}
*/
