# Artificial-Auditory-Simulation
This project is a software simulation of a cochlear implant, developed using   C/C++ and MATLAB to process audio and generate simulated stimulation signals without physical hardware.

This is a project report for an "Artificial Auditory Stimulation Device with C/C++ Implementation" submitted by A Sai Vardhan, M Srinivas, and V Srivathsa to Vasavi College of Engineering in May 2025. The project was completed to fulfill the academic requirements for a Bachelor of Engineering in Electronics & Communication Engineering.

# Project Overview
The project is a software-based prototype of a cochlear implant system. It was developed entirely in C and C++, and it simulates the key auditory processing mechanisms of a cochlear implant without the need for physical hardware. This makes the technology more accessible for academic research and eliminates the complexity and cost of using actual implants.

The simulation uses digital signal processing (DSP) techniques to process audio input:
Audio Acquisition: Captures audio from pre-recorded sources or in real time.
Bandpass Filtering: Decomposes the audio signal into different frequency bands, mimicking the cochlea's function.
Envelope Extraction: Extracts the amplitude envelope of each frequency band.
Amplitude Compression: Reduces the dynamic range of each envelope to manage loudness perception.
Stimulation Pattern Generation: Creates time-series pulse sequences based on the processed signals, which mimic the electrical stimulation patterns used in real cochlear implants.


# Implementation and Results
The project was implemented using two different platforms:
MATLAB: Used for algorithm prototyping, signal analysis, and visualization. It was used to process both recorded and real-time audio.
C/C++: Used for a more efficient and lightweight implementation, which was executed within a VirtualBox environment to ensure hardware independence.

The results show that the signal processing pipeline successfully simulates the core functions of a cochlear implant. The MATLAB implementation provided effective visualization and flexibility, while the C/C++ version proved to be efficient for offline analysis and data logging. The project demonstrates the feasibility of software-only auditory prosthetic simulation, laying the groundwork for future research in the field.

# General Details
pro-1.c it is recorded audio file processing in c/c++

real1.c it is a real time audio file processing in c/c++

project2.m it is recorded audio file processing in MATLAB

real4.m it is a real time audio file processing in MATLAB

for these recorded audio file processing we used sample.wav file
