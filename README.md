# EEG Aesthetics

This repository contains code for the presentation and analysis of EEG data collected during aesthetic perception experiments. The project explores the neural correlates of aesthetic experience using rapid visual presentation (RSVP) paradigms and representational similarity analysis (RSA).

## Table of Contents

- [Overview](#overview)  
- [Features](#features)  
- [Folder Structure](#folder-structure)  
- [Data Storage](#data-storage)  
- [Installation](#installation)  
- [Usage](#usage)  
- [Contributing](#contributing)  
- [License](#license)  

## Overview

This project investigates the relationship between brain activity and aesthetic preferences by:

- Presenting artwork stimuli under rapid presentation conditions  
- Preprocessing and analyzing EEG data  
- Performing representational similarity analysis (RSA) and correlation analyses  
- Studying associations between style, category, and liking ratings  

## Features

- RSVP task code for artwork stimuli presentation  
- Preprocessing pipeline including ICA  
- RSA with EEG data and model representations (3 Hz, 20 Hz)  
- Correlation analyses (simple, partial, style-category-liking)  
- Statistical modeling and FDR correction  

## Folder Structure
<details>
<summary><strong>📁 Folder Structure</strong></summary>
eeg_aesthetics/
├── RSVP_task/ # Code for stimuli presentation
│ ├── task_artworks_VG_faster_without_task.m # 50 ms presentation task
│ ├── task_artworks_VG_v2.m # 150 ms presentation task
│ ├── training_task_artworks_VG_faster_without_task.m # Training for 50 ms task
│ └── training_task_artworks_VG_v2.m # Training for 150 ms task
│
├── Analysis/ # Scripts for data analysis
│ ├── compute_average_RSA_eeg_vs_models_20HZ.m
│ ├── compute_average_RSA_eeg_vs_models_3HZ.m
│ ├── compute_average_correlations.m
│ ├── compute_style_category_liking_correlations.m
│ ├── fdr_bh.m
│ └── run_participant_analysis.m
│
├── data/ # Placeholder for local EEG datasets
├── results/ # Output of analyses (if locally stored)
├── notebooks/ # Jupyter notebooks (if any)
├── requirements.txt # Python/MATLAB dependencies
└── README.md # Project documentation
</details>


## Data Storage

- **Raw and processed EEG data:**  
  Stored securely on:
  - The lab computer  
  - University of Giessen Drive under the `Sanjeev` folder  

- **Participant log file:**  
  Contains details for all participants (e.g., participant number, condition, notes).  

- **Results:**  
  All results are stored in Google Drive, organized by participant:  
  [Results Google Drive Folder](https://drive.google.com/drive/u/0/folders/1cZUsK1hHitXB75Y79ct4GYoSpFLYd9Lt)  
  Each participant’s results are under folders named `eeg_aesthetics_<participant_number>`.  

## Installation

Clone the repository:

```bash
git clone https://github.com/vaishali2806/eeg_aesthetics.git
cd eeg_aesthetics

<details>
<summary><strong>🚀 Usage</strong></summary>

1. **Stimulus presentation:**  
   Run the appropriate MATLAB script in the `RSVP_task` folder:
   - `task_artworks_VG_faster_without_task.m` (50 ms)
   - `task_artworks_VG_v2.m` (150 ms)  
     Training versions are also available.

2. **Data analysis:**  
   Use `run_participant_analysis.m` for preprocessing, ICA, and core analyses.  
   Additional scripts for RSA and correlation analyses are provided in the `Analysis` folder.

3. **Results:**  
   View or download participant-specific results from the [Google Drive folder](https://drive.google.com/drive/u/0/folders/1cZUsK1hHitXB75Y79ct4GYoSpFLYd9Lt).

</details>

