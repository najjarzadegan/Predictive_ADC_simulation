Predictive ADC Simulator (MATLAB GUI)
This repository contains a specialized MATLAB-based simulation tool designed for modeling and analyzing Predictive Analog-to-Digital Converters. It provides a high-level interface to evaluate how hardware constraints and prediction algorithms impact ADC performance.

🚀 Features

Customizable Hardware Parameters: Adjust sub-ADC and sub-DAC resolutions (bits).
Define sampling rates and reference voltages.
Prediction Algorithms: Toggle between different prediction methods to compare tracking accuracy.
Time-Domain Analysis: Visualize real-time signal tracking and the behavior of internal nodes.
Frequency-Domain Analysis: Perform FFTs on the output to analyze SNR, ENOB, and harmonic distortion.
Intermediate Node Probing: Access waveforms from internal stages of the predictive loop for debugging and optimization.

🛠 Installation & Usage

Clone the repository: git clone https://github.com/najjarzadegan/Predictive_ADC_simulation.git
Open MATLAB and navigate to the project folder.
Run the application: Predictive_Two_step.mlapp
Configure & Simulate: Enter your desired resolutions and methods, then click Run Simulation.

📋 Requirements

MATLAB (R2020a or later recommended)
Signal Processing Toolbox

## Citation

If you use this simulator or the associated research in your work, please cite it as follows:

1. **Software:** Mohammad Najjarzadegan. (2026). *Predictive ADC Simulator MATLAB GUI* (Version 1.0) MATLAB. GitHub.(https://github.com/najjarzadegan/Predictive_ADC_simulation)
2. **Journal Paper:** Mohammad Najjarzadegan. (Submitted 2026). A Predictive Analog-to-Digital Converter (ADC) Architecture: Analysis, Design, and Application to Battery Cell Voltage Monitoring. *IEEE Transactions on Instrumentation and Measurement*.
3. **PhD Thesis:** Mohammad Najjarzadegan. (2026). *A high-performance multi-channel impedance spectroscopy system for fuel cell and battery diagnostic* PhD dissertation, University of British Columbia.
