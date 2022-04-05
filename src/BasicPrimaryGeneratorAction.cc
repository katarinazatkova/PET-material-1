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

// Intitializing the seeds 
G4int seed1 = 0;
G4int seed2 = 10000000;
G4int seed3 = 20000000;
G4int seed4 = 30000000;
G4int seed5 = 40000000;


void BasicPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  
  // Choose to comment out either case 1 or case 2 (see where they start and end)
  
  // THIS IS WHERE CASE 1 STARTS
  
  
  // Case 1: check that hits' distribution matches the expected
  // Gaussian distribution of parameters mu = 0 and sigma = (0.25/2.35) degrees.
  // This is because FWHM = 2.35*sigma = 0.25 degrees.  
  // By choice, the hits will go towards the opposite x-axis. One of them will
  // go straight towards (-1,0,0) direction. The other one will suffer the non-
  // collinearity effect from the annihilation event. It also contains different 
  // collinearity/noncollinearity cases for testing and resolution analysis purposes.
  
  // NOTE: 
  // When running this part, comment out the patient of the detector from BasicDetectorConstruction.cc file.
  // Only 'position', 'direction' and 'collinearity' integer parameters need to be changed to obtain different cases.
  
  // Choose from the following cases:
  // position 0:
    // fixed source position at (0,0,0)
    // direction 0: 
      // Implementation of different beam cases for +x direction photon (-x direction photon is fixed):
      // collinearity 0 - collinear, back-to-back
      // collinearity 1 - non-collinear, fixed direction at (0.25/2.35)*deg wrt photon 1
      // collinearity 2 - non-collinear, with Gaussian angular distribution 
    // direction 1:
      // Implementation for oblique events
      // collinearity 0 - collinear, back-to-back
      // collinearity 1 - non-collinear, with Gaussian angular distribution 
    // direction 2: 
      // Implementation for isotropic distribution
      // collinearity 0 - collinear, back-to-back
      // collinearity 1 - non-collinear, Gaussian angular distribution
  // position 1:
    // variable source position with isotropic beam directions
    // collinearity 0 - collinear, back-to-back
    // collinearity 1 - non-collinear, Gaussian angular distribution
  

  // Only these parameters need to be changed
  G4int position = 0;
  G4int direction = 0;
  G4int collinearity = 0;

  // We work with gamma rays of 511 keV 
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  fParticleGun->SetParticleEnergy(511*keV); 

  
  if (position == 0){
    // Vertex position
    G4double x0  = 0*cm, y0  = 0*cm, z0  = 0*cm;
    
    // Gamma rays have their origin at (0,0,0)
    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));

    if (direction == 0){
      G4ThreeVector dir0 = G4ThreeVector(-1, 0, 0);

      // Photon 1
      fParticleGun->SetParticleMomentumDirection(dir0);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      // Photon 2
      if (collinearity == 0){
        fParticleGun->SetParticleMomentumDirection(- dir0);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 1){
        G4double phi = (0.25/2.35) * twopi/360;
        G4ThreeVector dircol1 = G4ThreeVector(std::cos(phi), std::sin(phi), 0);

        fParticleGun->SetParticleMomentumDirection(dircol1);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
      
      else{
        // The second photon goes in positive x direction with a slight solid angle proportional 
        // to the Gaussian distribution angle with mu = 0 and sigma = (0.25/2.35) degrees.
        // This is because FWHM = 2.35*sigma = 0.25 degrees.
        G4double phi_gauss = twopi * G4RandGauss::shoot(0, (0.25/2.35)) / 360; 
        G4double theta = twopi * G4UniformRand(); 

        // Spherical polar coordinates are used 
        G4ThreeVector dircol2 = G4ThreeVector(std::cos(phi_gauss), 
                      std::sin(phi_gauss) * std::cos(theta), 
                      std::sin(phi_gauss) * std::sin(theta));  

        fParticleGun->SetParticleMomentumDirection(dircol2);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }

    else if (direction == 1){      
      // For simplicity y component is kept at zero
      // Choose if the simulation should be run for slightly or highly oblique events:

      // For slightly oblique events
      G4ThreeVector dir1 = G4ThreeVector(-1, 0, -1);

      // For highly oblique events that reach almost the end of the detector
      // detector length to diameter ratio is ≈ 2.5. 
      // To ensure that the particles hit the detector even in the noncollinear
      // beam case, choose z to x ratio of e.g. 2.4.
      //G4ThreeVector dir1 = G4ThreeVector(-1, 0, -2.4);

      // Photon 1
      fParticleGun->SetParticleMomentumDirection(dir1);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      // Photon 2
      if (collinearity == 0){
        fParticleGun->SetParticleMomentumDirection(- dir1);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 1){
        // Finding the angle about which to rotate
        G4double angle = dir1.angle(G4ThreeVector(-1, 0, 0));

        // Setting angles
        G4double phi_g = twopi * G4RandGauss::shoot(0,(0.25/2.35)) / 360;
        G4double theta_g = twopi * G4UniformRand(); 

        // Noncollinear vector placed around (1,0,0) direction
        G4ThreeVector antidir1 = G4ThreeVector(std::cos(phi_g), 
                      std::sin(phi_g) * std::cos(theta_g), 
                      std::sin(phi_g) * std::sin(theta_g));
        
        // This vector is then rotated around y axis
        antidir1.rotate(- angle, G4ThreeVector(0, 1, 0));

        fParticleGun->SetParticleMomentumDirection(antidir1);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }

    else if (direction == 2){
      // Photons 1 will point in random directions
      CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
      CLHEP::HepRandom::setTheSeed(seed1);
      G4double theta = twopi * G4UniformRand();
      seed1 += 1;
      
      CLHEP::HepRandom::setTheSeed(seed2);
      G4double phi = twopi * G4UniformRand();
      seed2 += 1;

      G4ThreeVector dir2 = G4ThreeVector(std::sin(theta) * std::cos(phi),
                    std::sin(theta) * std::sin(phi),
                    std::cos(theta));
      // Photon 1
      fParticleGun->SetParticleMomentumDirection(- dir2);
      fParticleGun->GeneratePrimaryVertex(anEvent);

      // Photon 2
      if (collinearity == 0){
        fParticleGun->SetParticleMomentumDirection(dir2);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }

      else if (collinearity == 1){
        // Setting angles
        G4double phi_g2 = twopi * G4RandGauss::shoot(0,(0.25/2.35)) / 360;
        G4double theta_g2 = twopi * G4UniformRand();        
        
        G4ThreeVector antidir2 = G4ThreeVector(std::cos(phi_g2), 
                      std::sin(phi_g2) * std::cos(theta_g2), 
                      std::sin(phi_g2) * std::sin(theta_g2));
        G4ThreeVector pos = G4ThreeVector(1, 0, 0);
        
        G4ThreeVector v = pos.cross(dir2);
        G4double a = pos.angle(dir2);
        
        antidir2.rotate(a, v);
        
        fParticleGun->SetParticleMomentumDirection(antidir2);
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }
  }


  else if (position == 1){
    // Uniform radiation spread throughout the cylinder
    G4double PET_radius = (78.6/2)*cm;
    G4double z_max = 1.94*m;
    G4double alpha_max = twopi;
    
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());


    // Introducing uniform randomness
    CLHEP::HepRandom::setTheSeed(seed1);
    G4double r = PET_radius * G4UniformRand();
    //seed1 += 1;

    //CLHEP::HepRandom::setTheSeed(seed2);
    G4double z = z_max * G4UniformRand() - 0.5 * z_max;
    //seed2 += 1;
    
    //CLHEP::HepRandom::setTheSeed(seed3);
    G4double alpha = alpha_max * G4UniformRand();
    //seed3 += 1;
  
    // Apply the spherical polar coordinates
    G4ThreeVector radiationOrigin = G4ThreeVector(r * std::cos(alpha), r * std::sin(alpha), z); 
    

    CLHEP::HepRandom::setTheSeed(seed4);;
    G4double theta_pos1 = twopi * G4UniformRand();
    seed4 += 1;
    
    CLHEP::HepRandom::setTheSeed(seed5);
    G4double phi_pos1 = twopi * G4UniformRand();
    seed5 += 1;

    G4ThreeVector dirpos1 = G4ThreeVector(std::sin(theta_pos1) * std::cos(phi_pos1),
                         std::sin(theta_pos1) * std::sin(phi_pos1),
                         std::cos(theta_pos1));


    // Setting the radiation origin
    fParticleGun->SetParticlePosition(radiationOrigin);
    
    // Photon 1
    fParticleGun->SetParticleMomentumDirection(- dirpos1);
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Photon 2
    if (collinearity == 0){
      fParticleGun->SetParticleMomentumDirection(dirpos1);
      fParticleGun->GeneratePrimaryVertex(anEvent);
      }

    else if (collinearity == 1){ 
      // Setting angles
      G4double phi_gpos1= twopi * G4RandGauss::shoot(0,(0.25/2.35)) / 360;
      G4double theta_gpos1 = twopi * G4UniformRand(); 

      G4ThreeVector antidirpos1 = G4ThreeVector(std::cos(phi_gpos1), 
                      std::sin(phi_gpos1) * std::cos(theta_gpos1), 
                      std::sin(phi_gpos1) * std::sin(theta_gpos1));
      G4ThreeVector pos1 = G4ThreeVector(1, 0, 0);

      G4ThreeVector v_pos1 = pos1.cross(dirpos1);
      G4double a_pos1 = pos1.angle(dirpos1);
        
      antidirpos1.rotate(a_pos1, v_pos1);
        
      fParticleGun->SetParticleMomentumDirection(antidirpos1);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }
  }

  // THIS IS WHERE CASE 1 ENDS
  
  
  // THIS IS WHERE CASE 2 STARTS
  /*
  
  // Case 2: now design a uniform radiation spread throughout the human phantom. The patient can be commented
  // out or not. The patient and detector are modelled as cyllinders, so cyllindrical polar coordinates will 
  // be used. Use the maximum allowed values. 
  // This part of the code was not used for the spatial resolution analysis. 
  // Note that the detector dimensions used here are not the ones of the uEXPLORER.
  
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
