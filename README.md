# Breast Cancer Gene Expression Analysis (HTAN Data Jamboree, Dec 4-7, 2023)

## Introduction
This project aims to identify differentially expressed genes in breast cancer using single-cell RNA sequencing data (scRNA). Understanding these gene expressions is pivotal in breast cancer research, diagnosis, and treatment planning.

## Goals and Objectives
The goal of this project is to identify differentially expressed genes in breast cancer using single-cell RNA sequencing data or bulk RNA sequencing data.  
- To process and analyze RNA sequencing data for breast cancer.
- Identify key differentially expressed genes associated with breast cancer.  
![overall_figure](assets/overall_figure.png)

## Technologies and Frameworks Used
- Python 3.x
- Datat proecessing (Scanpy)
- Differential gene expression analysis tools (DESeq2) ?
- Network inference analysis tools (GENIE3)
- Data visualization libraries (Matplotlib, Seaborn)

## Installation and Setup
1. Clone the repository:
   `git clone https://github.com/NCI-HTAN-Jamborees/Differential-Gene-Expression/tree/main](https://github.com/NCI-HTAN-Jamborees/Differential-Gene-Expression.git`
2. Install required packages:
   `pip install requirements.txt`

## Usage
`0_download_data.ipynb`: Jupiter notebook that helps pull data from Synapse and upload onto the Cancer Genomics Cloud.

## Data
This project uses scRNA sequencing data (level 3&4) from [HTAN](https://humantumoratlas.org/). The two datasets of interest are HTAN HTAPP and HTAP WUSTL.  
Different contrast groups were utilized for the different datasets, namely age and diagnosis for HTAPP and race and disease progression/recurrence for WUSTL.

## Results and Analysis
The analysis outputs a list of differentially expressed genes, which can be found in 

## License
This project is licensed under the [MIT](https://github.com/NCI-HTAN-Jamborees/Differential-Gene-Expression/blob/main/LICENSE).

## References

Huynh-Thu et al, PLoS One 2010: "Inferring regulatory networks from expression data using tree-based methods"

## Team Members

- **Team Lead**: Rebecca Sansale
- **Data Analysts**:
  - Giovanni Dall’Olio
  - Shivani Guturu
  - Minji Kim (also Illustrator)
  - Jiaying Lai
- **Tech Lead & Analyst**: Nirvana Nursimulu
- **Writer & Advisor**: Tapasree Roy Sarkar
