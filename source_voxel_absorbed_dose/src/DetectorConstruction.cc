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
// $Id$
//
/// \file VoxelSValues.cc
/// \brief Implementation of the VoxelSValues class

// includes for application written by end user
#include "DetectorConstruction.hh"
// use of stepping action to set the accounting volume
#include "SteppingAction.hh"


//g4 includes (written by g4 team)
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelSValuesDetectorConstruction::VoxelSValuesDetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelSValuesDetectorConstruction::~VoxelSValuesDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* VoxelSValuesDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Phantom parameters
  //  30 cm water cube
  G4double phantom_sizeXY = 30*cm, phantom_sizeZ = 30*cm;
  G4Material* phantom_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*phantom_sizeXY;
  G4double world_sizeZ  = 1.2*phantom_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Phantom
  //  
  G4Box* solidPhantom =    
    new G4Box("Phantom",                    //its name
        0.5*phantom_sizeXY, 0.5*phantom_sizeXY, 0.5*phantom_sizeZ); //its size
      
  G4LogicalVolume* logicPhantom =                         
    new G4LogicalVolume(solidPhantom,            //its solid
                        phantom_mat,             //its material
                        "Phantom");         //its name

  G4VisAttributes* phantomVisAtt = new G4VisAttributes(G4Colour::Green());
  phantomVisAtt->SetForceWireframe(true);
  logicPhantom->SetVisAttributes(phantomVisAtt);

               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicPhantom,                //its logical volume
                    "Phantom",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     
  // source voxel
  //   We explicitly create a source voxel for scoring purposes in line with example B1.
  //
  
  G4Material* sourceVoxelMat = nist->FindOrBuildMaterial("G4_WATER");
  G4ThreeVector pos1 = G4ThreeVector(0*mm, 0*cm, 0*cm);
        
  // Box shape for isotropic voxel
  G4double dx=3*mm;
  G4Box* solidVoxel = 
    new G4Box( "solidVoxel", 0.5*dx, 0.5*dx, 0.5*dx);

  G4LogicalVolume* logicVoxel =                         
    new G4LogicalVolume(solidVoxel,         //its solid
                        sourceVoxelMat,          //its material
                        "solidVoxel");           //its name
  G4VisAttributes* voxelVisAtt = new G4VisAttributes(G4Colour::White());
  voxelVisAtt->SetForceSolid(true);
  logicVoxel->SetVisAttributes(voxelVisAtt);
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicVoxel,             //its logical volume
                    "solidVoxel",                //its name
                    logicPhantom,            //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
  // Set scoring volume to stepping action 
  // (where we will account energy deposit)
  //
  VoxelSValuesSteppingAction* steppingAction = VoxelSValuesSteppingAction::Instance(); 
  ////steppingAction->SetVolume(logicShape1);
  steppingAction->SetVolume(logicVoxel);

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
