# OpenUQFOAM

## Overview
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nullam nec purus nec nunc.

## Project Structure
```
OpenUQFOAM
├── data
│   ├── traces       # Sample traces from inverse UQ
│   ├── plots        # Plots generated from analysis (should be improved)
│   └── <other-data-files>
├── notebooks
│   ├── 01_sampling.ipynb     # Sampling step, where the simulator is called to generate samples and outputs stored
│   ├── 02_propagation.ipynb  # Propagation step, where the saved samples are used to propagate the uncertainty
├── src
│   ├── utils.py           # Utility functions
│   └── <source-files>
├── simulator              # Simulator (submodule), where the openFOAM simulation runs
│   └── <simulator-files>
├── requirements.txt       # Project dependencies
```

## Setup Instructions
1. Clone the repository:
   ```
   git clone <repository-url>
   ```
2. Navigate to the project directory:
   ```
   cd OpenUQFOAM
   ```
3. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

