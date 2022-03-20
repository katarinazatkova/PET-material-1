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

#include "BasicPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include <CLHEP/Random/RandGaussQ.h>

// PGA
// based on Tangle2's back2back photons

BasicPrimaryGeneratorAction::BasicPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
}

//

BasicPrimaryGeneratorAction::~BasicPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//

void BasicPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  
  // Choose to comment out either case 1 or case 2 (see where they start and end)
  
  // THIS IS WHERE CASE 1 STARTS
  
  
  // Case 1: check that hits' distribution matches the expected
  // Gauss distribution of parameters mu = 0 and sigma = 0.25 degrees. 
  // By choice, the hits will go towards the opposite x-axis. One of them will
  // go straight towards (-1,0,0) direction. The other one will suffer the non-
  // collinearity effect from the annihilation event. It also contains different 
  // collinearity cases for testing and resolution analysis purposes.
  
  // NOTE: 
  // When running this part, comment out the patient of the detector from BasicDetectorConstruction.cc file.
  // Only 'position', 'direction' and 'collinearity' integer parameters need to be changed to obtain different cases.
  
  // Choose from the following cases:
  // position 0:
    // fixed source position at (0,0,0)
    // direction 0: 
      // Implementation of different beam cases for +x direction photon (-x direction photon is fixed):
      // collinearity 0 - back-to-back
      // collinearity 1 - fixed direction at 5 deg wrt photon 1
      // collinearity 2 - fixed angle around beam axis
      // collinearity 3 - with Gaussian angular distribution 
    // direction 1:
      // Implementation for oblique events
      // collinearity 0 - collinear, back-to-back
      // collinearity 1 - noncollinear, with Gaussian angular distribution 
    // direction 2: 
      // Implementation for isotropic distribution
      // collinearity 0 - collinear, back-to-back
      // collinearity 1 - noncollinear, Gaussian angular distribution
  // position 1:
    // variable source position with isotropic beam directions
    // collinearity 0 - collinear, back-to-back
    // collinearity 1 - noncollinear, Gaussian angular distribution
  

  // only these parameters need to be changed
  int position = 0;
  int direction = 1;
  int collinearity = 1; 

  
  if (position == 0){
    // Get a vertex
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    
    fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->SetParticleEnergy(511*keV);
  

    if (direction == 0){
      // We work with gamma rays of 511 keV with the origin at (0,0,0) vector
      
      G4ThreeVector dir0 = G4ThreeVector(-1, 0, 0);
      // Photon 1
      fParticleGun->SetParticleMomentumDirection(dir0);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      // Photon 2
      if (collinearity == 0){
        fParticleGun->SetParticleMomentumDirection(-dir0);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 1){
        G4double phi = 5 * twopi/360;
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(std::cos(phi), std::sin(phi), 0));
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 2){         
        G4double phi = 5 * twopi/360;    
        G4double theta = twopi * G4UniformRand();

        // Direction for the second photon: 
        // It goes in positive x-direction with a slight solid angle for y and z non-zero values.
        // Spherical polar coordinates are used
        G4ThreeVector dircol2 = G4ThreeVector(std::cos(phi), 
                      std::sin(phi) * std::cos(theta), 
                      std::sin(phi) * std::sin(theta));  

        fParticleGun->SetParticleMomentumDirection(dircol2);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
      
      else{
        // Forming a cone of solid angle proportional to the gauss distribution angle
        // with mu = 0 and sigma = 0.25 degrees.
        G4double phi_gauss = twopi * G4RandGauss::shoot(0,0.25) / 360; 
        G4double theta = twopi * G4UniformRand(); 
        
        // Direction for the second photon: 
        // It goes in positive x-direction with a slight solid angle proportional
        // to the gauss distribution angle with mu = 0 and sigma = 0.25 degrees.
        // Spherical polar coordinates are used again
        G4ThreeVector dircol3 = G4ThreeVector(std::cos(phi_gauss), 
                      std::sin(phi_gauss) * std::cos(theta), 
                      std::sin(phi_gauss) * std::sin(theta));  

        fParticleGun->SetParticleMomentumDirection(dircol3);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }


    else if (direction == 1){
      // for simplicity y is kept at zero
      
      // for slightly oblique events
      //G4ThreeVector dir = G4ThreeVector(-1, 0, -1);

      // for highly oblique events that reach almost the end of the detector
      // detector length to diameter ratio is ~2.5, 
      // to ensure that the particles hit the detector even in the noncollinear
      // beam case, choose z to x ratio of e.g. 2.3
      G4ThreeVector dir1 = G4ThreeVector(-1, 0, -2.3);

      // Photon 1
      fParticleGun->SetParticleMomentumDirection(dir1);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      // Photon 2
      if (collinearity == 0){
        fParticleGun->SetParticleMomentumDirection(-dir1);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 1){
        // finding the angle about which to rotate
        G4double angle = dir1.angle(G4ThreeVector(-1, 0, 0));

        // setting angles
        G4double phi_g = twopi * G4RandGauss::shoot(0,0.25) / 360;
        G4double theta_g = twopi * G4UniformRand(); 

        // noncollinear vector placed around (1,0,0)
        G4ThreeVector antidir = G4ThreeVector(std::cos(phi_g), 
                      std::sin(phi_g) * std::cos(theta_g), 
                      std::sin(phi_g) * std::sin(theta_g));
        
        // this vector is then rotated around y axis
        antidir.rotate(-angle,G4ThreeVector(0,1,0));

        fParticleGun->SetParticleMomentumDirection(antidir);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }


    else if (direction == 2){
      G4double theta = twopi * G4UniformRand();
      G4double phi = twopi * G4UniformRand();
      
      G4ThreeVector dir2 = G4ThreeVector(std::cos(theta) * std::sin(phi), std::cos(theta) * std::cos(phi), 
                    std::sin(theta));
      
      // Photon 1
      fParticleGun->SetParticleMomentumDirection(dir2);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      // Photon 2
      if (collinearity == 0){
        G4ThreeVector dir2nogauss = G4ThreeVector(-std::cos(theta) * std::sin(phi), -std::cos(theta) * std::cos(phi),
                    -std::sin(theta));
        
        fParticleGun->SetParticleMomentumDirection(dir2nogauss);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 1){
        G4double gauss_dev = G4RandGauss::shoot(0, 0.25) * twopi / 360;
        
        G4ThreeVector dir2gauss = G4ThreeVector(-std::cos(theta + gauss_dev) * std::sin(phi), -std::cos(theta + gauss_dev) * std::cos(phi),
                    -std::sin(theta + gauss_dev));
        
        fParticleGun->SetParticleMomentumDirection(dir2gauss);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }
  }


  else if (position == 1){
    // uniform radiation spread throughout the cylinder
    G4double PET_radius = (78.6/2)*cm;
    G4double z_max = 1.94*m;
    G4double alpha_max = twopi;
  
    // Introducing uniform randomness
    G4double r = PET_radius * G4UniformRand();
    G4double z = z_max * G4UniformRand() - 0.5 * z_max;
    G4double alpha = alpha_max * G4UniformRand();
  
    // Apply the spherical polar coordinates
    G4ThreeVector radiationOrigin = G4ThreeVector(r * std::cos(alpha), r * std::sin(alpha), z); 

    //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    //CLHEP::HepRandom::setTheSeed(0);

    // angles  
    G4double theta = twopi * G4UniformRand();
    G4double phi = twopi * G4UniformRand();
    G4double gauss_dev = G4RandGauss::shoot(0,0.25) * twopi / 360;
    
    // Now use the spherical polar coordinates two produce two beams in nearly opposite directions
    G4ThreeVector dirpos1 = G4ThreeVector(std::cos(theta) * std::sin(phi), std::cos(theta) * std::cos(phi), 
                  std::sin(theta));
    
    fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
      
    // Photon 1
    fParticleGun->SetParticlePosition(radiationOrigin);
    fParticleGun->SetParticleEnergy(511*keV);
    fParticleGun->SetParticleMomentumDirection(dirpos1);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Photon 2
    if (collinearity == 0){
      G4ThreeVector dirpos1nogauss = G4ThreeVector(-std::cos(theta) * std::sin(phi), -std::cos(theta) * std::cos(phi),
                  -std::sin(theta));
      
      fParticleGun->SetParticleMomentumDirection(dirpos1nogauss);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }

    else if (collinearity == 1){
      G4ThreeVector dirpos1gauss = G4ThreeVector(-std::cos(theta + gauss_dev) * std::sin(phi), -std::cos(theta + gauss_dev) * std::cos(phi),
                  -std::sin(theta + gauss_dev));
      
      fParticleGun->SetParticleMomentumDirection(dirpos1gauss);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }
  }

  // THIS IS WHERE CASE 1 ENDS
  
  
  // THIS IS WHERE CASE 2 STARTS
  /*
  
  // Case 2: now design a uniform radiation spread throughout the human phantom. The patient can be commented
  // out or not. The patient and detector are modelled as cyllinders, so cyllindrical polar coordinates will 
  // be used. Use the maximum allowed values.
  
  G4double PET_radius = 15*cm;
  G4double z_max = 1.95*m;
  G4double alpha_max = twopi;
  
  // Now introduce uniform randomness for getting the gamma-radiation inside the human phantom
  
  G4double r = PET_radius * G4UniformRand();
  G4double z = z_max * G4UniformRand() - 0.5 * z_max;
  G4double alpha = alpha_max * G4UniformRand();
  
  // Apply the spherical polar coordinates
  
  G4ThreeVector radiationOrigin = G4ThreeVector(r * std::cos(alpha), r * std::sin(alpha), z);
  
  // Now use again the non-collinearity effect. This time, however, the origin for each event
  // becomes variable. Introduce the random angles for the spherical polar coordinates and the
  // additional gaussian variable
  
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  // CLHEP::HepRandom::setTheSeed(0);
  
  G4double theta = twopi * G4UniformRand();
  G4double phi = twopi * G4UniformRand();
  G4double gauss_dev = G4RandGauss::shoot(0,0.25) * twopi / 360;
  // G4double gauss_dev = G4RandGauss::shoot(0, 0.25);
  // G4double gauss_dev = 0;
  
  G4cout << "Value of theta angle: " << theta << G4endl;
  G4cout << "Value of phi angle: " << phi << G4endl;
  G4cout << "Value of the gaussian deviation is: " << gauss_dev << G4endl;
  
  // Now use the spherical polar coordinates two produce two beams in nearly opposite directions
  
  G4ThreeVector photonDir = G4ThreeVector(std::cos(theta) * std::sin(phi), std::cos(theta) * std::cos(phi), 
  					   std::sin(theta));
  G4ThreeVector photonAntiDir = G4ThreeVector(-std::cos(theta + gauss_dev) * std::sin(phi), -std::cos(theta + gauss_dev) * std::cos(phi),
  						-std::sin(theta + gauss_dev));
  						
  // Set the two outcoming photons
  
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  
  // Photon 1
  
  fParticleGun->SetParticlePosition(radiationOrigin);
  fParticleGun->SetParticleEnergy(511*keV);
  fParticleGun->SetParticleMomentumDirection(photonDir);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  // Photon 2
  
  fParticleGun->SetParticlePosition(radiationOrigin);
  fParticleGun->SetParticleEnergy(511*keV);
  fParticleGun->SetParticleMomentumDirection(photonAntiDir);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  */
  // THIS IS WHERE CASE 2 ENDS
  
}
