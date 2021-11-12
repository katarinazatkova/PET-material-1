# Alexandru Mihai Hau
Jun 2021 - September 2021

# Research Project on Total Body PET Scanner

This project is about comparing the sensitivity and signal-to-noise ratio of a Positron Emission Tomography (PET) Scanner at different lengths (those used nowadays in industry, of around 10-15cm length, and the Total Body PET Scanners at 2m length). There are two steps involved in the project, both of them of great importance. For the first simulation, a small beam on positive and negative x-axis has been considered, with angular deviation of 180 degrees and the additional gaussian distribution of mean 0 and standard deviation 0.25 degrees. By locating the hits on detector, the purpose of this step is to check that computer simulations match the theory. Moreover, this step is important for future work in spatial resolution. The second step of the project involves creating random radiation distribution in order to simulate the sensitivity and SNR parameters.

## Geometry definition

A cyllinder has been simulated as a PET detector. The cyllinder has been chosen at two lengths: 15 cm and 2 meters. A patient has also been placed inside. The patient can be commented out in the BasicDetector.cc file, by removing the lines of introducing the G4PhysicalVolume patient object into the physical world. The material from which the patient is made can also be changed. For the first step, the detector should be declared very thick (30 cm thickness is a good choice), and the patient should be commented out. For the second part, thickness should be smaller, in order to consider high-energy beams which escape the detector.

## Primary Generator Action

This file is used for creating the outgoing beam which represents the radiation inside the human body from the positron-electron annihilation event. Two gamma rays of 511 keV are created. Moreover, this is the file where the two steps of the problem are implemented.

## Event Action

This file represents the methods implemented at the end of each event. Any event with deposited energy greater than 0.9 MeV is considered a good event (an event which helps in reconstructing the location of the positron - electron annihilation event). Each event consists of a collection of hits (G4HitsCollection), whose parameters can be accessed through BasicPetHit.cc and BasicPetHit.hh files. Warning: Do not change the AddEdep method of the hits collection. At the end of each event, the desired parameters are added into the histogram for event collection of the run

## Run Action

The histrogams and plots are created for the parameters to be analysed. Moreover, the sensitivity and SNR are printed out.
