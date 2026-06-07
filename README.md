# Shine DalGuardYES
**Computational design of RBS architectures for orthogonal anti–Shine-Dalgarno ribosome systems**  
*IGEM — Université Evry Paris-Saclay*  
###### By [@PaulVerot03](https://github.com/PaulVerot03),  [@VoxelMC](https://github.com/VoxelMC), [@cuajiniquil](https://github.com/cuajiniquil), [@ritar18](https://github.com/ritar18), [@georgyzaouk](https://github.com/georgyzaouk), [@serenapandzou](https://github.com/serenapandzou), [@talissakassably](https://github.com/talissakassably) & [@crakshay1](https://github.com/crakshay1)

-----------------------------------------------------
**Shine DalGuardYES** designs orthogonal ribosome binding sites from a user-defined anti-Shine-Dalgarno sequence. Unlike standard RBS calculators, which evaluate natural ribosome systems, our model performs inverse design: it searches possible RBS motifs and spacer lengths to find sequences that maximize recognition by the engineered ribosome while minimizing recognition by the endogenous ones.

#### [Presenting DalGuardYES: Computational design of RBS for orthogonal anti–Shine-Dalgarno ribosomes. (Slides only)](https://canva.link/uuod8ao0w4xcydx)
-------

## Project Overview

### Preview
<img width="1920" height="1047" alt="image" src="https://github.com/user-attachments/assets/f2ed2f93-02c0-4419-b051-5477a9083be9" />
<img width="1920" height="1047" alt="image" src="https://github.com/user-attachments/assets/6d476235-1508-496e-8498-1910194801e2" />

### Workflow
<img width="499" height="743" alt="image" src="https://github.com/user-attachments/assets/afeb3068-dfb6-4162-966b-7fa35461ab0d" />

## Repository Structure

```text
  DalGuardYES-UI/                         # Source code
  ├── candidates.py                       # Initial candidates 
  ├── orbs_duplex.py                      # Core detection 
  ├── riboguard_ga_engine_clean.py        # Genetic algorithm application
  ├── riboguard_streamlit.py              # Dashboard
  └── requirements.txt                    # Requirements
  _Archives/                              # Old code
```
