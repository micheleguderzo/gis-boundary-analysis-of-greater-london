# GIS boundary analysis of Greater London – Spatial Econometrics Project

![Language](https://img.shields.io/badge/code-R-blue?logo=r&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green?logo=open-source-initiative)
![Status](https://img.shields.io/badge/status-finished-success?style=flat&logo=github)

This repository contains the final project for the **Spatial Econometrics in R** course, carried out during the academic year 2021/2022 at the Faculty of Economic Sciences, University of Warsaw.

## 📌 Project Overview

This project presents a **geospatial analysis** of Greater London based on **GIS boundary data**, aiming to understand whether and how the current structure of boroughs is influenced by their geographical location. By using shapefiles of London boroughs and wards, the analysis examines the spatial relationships and structural characteristics that define London's urban geography.

## 🧠 Key Topics

- Spatial data import and analysis  
- Choropleth mapping and spatial visualization  
- OLS regression and log-linear regression  

## 🗂️ Repository Structure

├── README.md     # This file  
├── LICENSE       # Project license (MIT)  
├── .gitignore    # Git ignore file for excluded files/folders  
├── code/         # R scripts for spatial analysis, distance calculations, and regression modeling  
├── data/         # Shapefiles and spatial datasets from the London Datastore  
├── plots/        # Choropleth maps and scatter plot
└── report/       # Final report (PDF written in LaTeX)

## 📊 Dataset

The dataset was obtained from the official [London Datastore website](https://data.london.gov.uk/), and includes key GIS boundary files for Greater London. Among the available shape files, only two were used for this project:

- **London Wards** (with City of London merged into a single area)
- **London Boroughs**

## 🛠 Tools & Libraries

- R  
- RStudio  
- LaTeX  
- R packages used:
  - `spdep`
  - `rgdal`
  - `maptools`
  - `sp`
  - `RColorBrewer`
  - `classInt`
  - `GISTools`
  - `maps`
  - `geosphere`
  - `moments`
  - `dplyr`
  - `tseries`
  - `epiDisplay`
  - `tidyverse`
  - `car`
  - `lmtest`
  - `MASS`
  - `olsrr`
  - `haven`
  - `ggplot2`
  - `ggpubr`

## 📄 Final Report

You can read the full report with methodology, analysis, and conclusions in [`report/final_report.pdf`](report/final_report.pdf).

## 👤 Author

**Michele Guderzo**  

## 📝 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

*This project was developed for educational purposes only.*
