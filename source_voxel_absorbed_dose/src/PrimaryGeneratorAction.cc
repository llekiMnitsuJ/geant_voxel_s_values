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
//  Modified by Justin Mikell for isotropic spatial and directional sampling
//  inside the source voxel. 
//
// $Id$
//
/// \file VoxelSValuesPrimaryGeneratorAction.cc
/// \brief Implementation of the VoxelSValuesPrimaryGeneratorAction class

//includes written for application by end user
#include "PrimaryGeneratorAction.hh"


//includes written by g4 team
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelSValuesPrimaryGeneratorAction* VoxelSValuesPrimaryGeneratorAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const VoxelSValuesPrimaryGeneratorAction* VoxelSValuesPrimaryGeneratorAction::Instance()
{
// Static acces function via G4RunManager 

  return fgInstance;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelSValuesPrimaryGeneratorAction::VoxelSValuesPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(0.01*MeV);
  
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelSValuesPrimaryGeneratorAction::~VoxelSValuesPrimaryGeneratorAction()
{
  delete fParticleGun;
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelSValuesPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the beginning of each event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Phantom volume
  // from G4LogicalVolumeStore.
  

  G4double dXYZ = 0;
  G4LogicalVolume* voxelLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("solidVoxel");
  G4Box* voxelBox = NULL;
  if ( voxelLV ) voxelBox = dynamic_cast<G4Box*>(voxelLV->GetSolid());
  if ( voxelBox ) {
    dXYZ = voxelBox->GetXHalfLength()*2.;
  }  
  else  {
    G4cerr << "source voxel shape not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  //add energy sampling here...

  //sample position within the source voxel
  G4double x0 = dXYZ * (G4UniformRand()-0.5);
  G4double y0 = dXYZ * (G4UniformRand()-0.5);
  G4double z0 = dXYZ * (G4UniformRand()-0.5);
  
  //sample isotropic direction 
  G4double cosTheta = -1.0 + 2.0*G4UniformRand();
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4double sinTheta = sqrt(1. - cosTheta*cosTheta);
  
  G4double px0 = sinTheta*cos(phi);
  G4double py0 = sinTheta*sin(phi);
  G4double pz0 = cosTheta;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px0,py0,pz0));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

