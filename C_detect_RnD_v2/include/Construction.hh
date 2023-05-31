#ifndef CandyConstruction_h
#define CandyConstruction_h

#include "G4VUserDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"


#include "G4MaterialTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4VisAttributes.hh" //World 안보이게 하려고 

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"

//#include "G4VSolid.hh"
//#include "G4PhysicalConstants.hh"
//#include "G4SDManager.hh"
//#include "G4ProcessManager.hh"  // For managing processes in a material
//#include "G4Scintillation.hh"  // For the scintillation process
//#include "G4RotationMatrix.hh"
//#include "G4ThreeVector.hh"



class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  MyDetectorConstruction();
  ~MyDetectorConstruction() override;
  
  G4VPhysicalVolume* Construct() override;

  G4LogicalVolume* GetSiPMVolume() const { return fSipm_SD_LOG; }
  G4LogicalVolume* GetPMTVolume() const { return fPMT_SD_LOG; }
  G4LogicalVolume* GetOpticla_Fiber_Volume() const { return fOptical_Fiber_LOG; }

protected:
  G4LogicalVolume* fSipm_SD_LOG = nullptr;
  G4LogicalVolume* fPMT_SD_LOG = nullptr;
  G4LogicalVolume* fOptical_Fiber_LOG = nullptr;

private:
  G4MaterialPropertiesTable* fmat_Prop_table_LAB;
  //G4MaterialPropertiesTable* fmat_Prop_table_Teflon;
  G4Material* fMat_LAB;

 

  
  // virtual void ConstructSD();
};
#endif
