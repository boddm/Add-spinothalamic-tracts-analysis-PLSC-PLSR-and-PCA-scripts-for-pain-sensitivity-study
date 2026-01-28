# Spinothalamic Tract Microstructure Analysis Pipeline

Analysis scripts for the paper **"Spinothalamic Tract Microstructure as a Common Neural Substrate of Pain Sensitivity Across Modalities."**

## System Requirements

**Software**
- MATLAB R2018b or newer
- Statistics and Machine Learning Toolbox

**Computer**
- Windows 10/11, macOS 10.14+, or Linux
- 8 GB RAM minimum (16 GB recommended)
- 2 GB free disk space

**Tested on**: MATLAB R2018b-R2022a, Windows 10, macOS Monterey, Ubuntu 20.04

---

## Installation Guide

1. Download or clone this repository
2. Open MATLAB and navigate to the project folder
3. Run:
   ```matlab
   addpath(genpath('1_PLSC_Script'))
   addpath(genpath('2_Prediction_Script'))
   savepath
   ```
4. Verify: `which myPLS_main` should show the file path

**Install time**: ~2-5 minutes on a typical desktop

---

## Demo

### Quick Start with Sample Data

**Demo 1: PLSC Analysis**
```matlab
cd 1_PLSC_Script
% Edit myPLS_inputs.m to load sample data from 0_Data/
myPLS_main
```

**Demo 2: PLSC+PLSR Prediction**
```matlab
cd 2_Prediction_Script/1_PLSC_PLSR_original
mkdir example/1/Processing example/1/Results_PLSR
copyfile('../../0_Data/data1_pain.mat', 'example/1/Processing/data.mat')
S1_PLSC_PLSR_real
```

**Demo 3: PLSC+PCA Analysis**
```matlab
cd 2_Prediction_Script/2_PLSC_PCR_original
S2_PLSC_PCR_real
```

**Expected run time** (Intel i7, 16GB RAM):
- Small dataset (50 subjects): 5-10 min
- Medium dataset (100 subjects): 15-30 min
- Large dataset (200+ subjects): 45-90 min

---

## Instructions for Use

### Step 1: Prepare Your Data

Create a `.mat` file with:
- `brain_data`: subjects × imaging features matrix
- `beh_data`: subjects × behavioral variables matrix
- `diagnosis` (optional): grouping vector

```matlab
save('data.mat', 'brain_data', 'beh_data', 'diagnosis')
```

### Step 2: Configure Analysis

**For PLSC**: Edit `1_PLSC_Script/myPLS_inputs.m`
- Set data path: `load('your/data.mat')`

**For PLSC+PLSR/PCR**: Edit the script in `2_Prediction_Script/`
- Set paths to your data and output folders

### Step 3: Run Analysis

```matlab
% PLSC
cd 1_PLSC_Script
myPLS_main

% PLSC+PLSR
cd 2_Prediction_Script/1_PLSC_PLSR_original
S1_PLSC_PLSR_real

% PLSC+PCR
cd 2_Prediction_Script/2_PLSC_PCR_original
S2_PLSC_PCR_real
```

% Prediction testing
cd 2_Prediction_Script/3_PLSC_PCR_original_noCV
S3_1_PLSC_PCR_final_model
S3_2_test_final_model_newmethod_for_mult

### Step 4: Results

**Output files**: Saved as `.mat` files in your output directory

**Visualizations**: Charts and plots display automatically

**Key metrics**:
- PLSC: Components, singular values, bootstrap ratios, p-values
- PLSR/PCR: R² values, cross-validation scores, feature weights

### Troubleshooting

| Problem | Solution |
|---------|----------|
| Out of memory | Reduce nPerms/nBootstraps, close other programs |
| File not found | Check file paths, use full paths |
| Dimension mismatch | Ensure brain_data and beh_data have same number of subjects |
| Running slowly | Reduce permutations, use smaller test dataset |
| No significant results | Check data quality, increase sample size |

---

## Citation

Please cite our paper when using these tools:
```
[Paper citation to be added upon publication]
```

## License & Contact

[License information to be added]

For questions: [contact information]

**Last Updated**: January 2026
