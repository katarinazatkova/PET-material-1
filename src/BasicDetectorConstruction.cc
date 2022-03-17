//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "BasicDetectorConstruction.hh"
#include "BasicPETSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSDoseDeposit.hh"

#include "G4SDManager.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//

//G4ThreadLocal

//

BasicDetectorConstruction::BasicDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//

BasicDetectorConstruction::~BasicDetectorConstruction()
{
}

//

G4VPhysicalVolume* BasicDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//

void BasicDetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes = false;
  
  // LSO
  /*
  G4Element*  O = man->FindOrBuildElement("O" , isotopes);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5);
  */
  
  // LYSO
  G4Element*  O = man->FindOrBuildElement("O" , isotopes);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);  
  G4Element* Ce = man->FindOrBuildElement("Ce", isotopes);		 
  G4Element*  Y = man->FindOrBuildElement("Y" , isotopes);

  G4Material* LYSO = new G4Material("LYSO", 7.1*g/cm3, 5);
  LYSO->AddElement(Lu, 71.43*perCent);
  LYSO->AddElement(Y, 4.03*perCent);
  LYSO->AddElement(Si, 6.37*perCent);
  LYSO->AddElement(O, 18.14*perCent);
  LYSO->AddElement(Ce, 0.02*perCent);
}

//

G4VPhysicalVolume* BasicDetectorConstruction::DefineVolumes()
{
  // we'll need some air
  G4NistManager* nist = G4NistManager::Instance();
   // pick material by commenting out the rest
  G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_MUSCLE_STRIATED_ICRU");
  /*
   G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
  G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_SKIN_ICRP");
  G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_BLOOD_ICRP");
  G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_ADIPOSE_TISSUE_ICRP");
 */
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* tube_mat   = nist->FindOrBuildMaterial("LYSO");

  // PET dimensions
  G4double PET_in_rad = (78.6/2)*cm, PET_out_rad = PET_in_rad + (CrystLength)*cm, PET_length = 1.94*m; // vary these

  // world size
  G4double world_dim = 4*m;

  //
  // World
  //
  G4Box* solidWorld =
  new G4Box("World",
    0.5*world_dim, 0.5*world_dim, 0.5*world_dim);

// fill it with air

  G4LogicalVolume* logicWorld =
  new G4LogicalVolume(solidWorld,
                      default_mat,
                      "World");

// put it somewhere

  G4VPhysicalVolume* physWorld =
  new G4PVPlacement(0,  // no rotation
                      G4ThreeVector(),  // set origin
                      logicWorld,  // logical volume
                      "World",    // name
                      0,    // mother volume
                      false,  // boolean
                      0,    // copy number
                      fCheckOverlaps);    // checks for volume overlaps

// make cylinder

  G4Tubs* solidTube =
    new G4Tubs("Tube", PET_in_rad, PET_out_rad, 0.5*PET_length, 0., twopi);

// fill it with LYSO

  G4LogicalVolume* logicTube =
  new G4LogicalVolume(solidTube,
                      tube_mat,
                      "Tube");

// put it in the world

new G4PVPlacement(0,                       // no rotation
                  G4ThreeVector(),         // set origin
                  logicTube,           // logical volume
                  "Tube",              // name
                  logicWorld,              // mother  volume
                  false,                   // boolean operation
                  0,                       // number
                  fCheckOverlaps);         // checking overlaps
  
  //
  // patient
  //
  G4double patient_radius = (76/2)*cm;
  G4double patient_dZ = 1.94*m;  
    
  G4Tubs* solidPatient =
    new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi);
      
//   G4LogicalVolume* logicPatient =                         
//     new G4LogicalVolume(solidPatient,        //its solid
//                         phantom_mat,         //its material
//                         "Patient");        //its name
  G4LogicalVolume* logicPatient =                         
    new G4LogicalVolume(solidPatient,        //its solid
                        phantom_mat,         //its material
                        "Patient");        //its name
               
  //
  // place patient in world
  //
  /*                  
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicPatient,            //its logical volume
                    "Patient",               //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 
  
  */
  return physWorld;
}

//

void BasicDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // Sensitive detector
  auto detectorSD
    = new BasicPETSD("detectorSD", "DetectorHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(detectorSD);
  SetSensitiveDetector("Tube",detectorSD);

  // Make phantom a sensitive detector 

  
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  // the next line was missing which was causing a segmentation fault
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  //
  G4VPrimitiveScorer* primitiv2 = new G4PSEnergyDeposit("edep");
  patient->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("Patient",patient);
  
  
  /*
  auto phantomSD
    = new BasicPETSD("phantomSD", "PatientHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(phantomSD);
  SetSensitiveDetector("Patient",phantomSD);
  */
 
}

//
