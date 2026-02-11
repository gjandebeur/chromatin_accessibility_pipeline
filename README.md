

## Chromatin Accessibility and Microexon Splicing Regulation in Alzheimer's Disease

A 6-week rotation project investigating chromatin regulation of microexon splicing in Alzheimer's disease using ENCODE RUSH AD data.


### Project Overview 

This project develops a SACS (Splicing-Associated Chromatin Score) scoring system to investigate microexon splicing regulation in Alzheimer's disease using ENCODE RUSH AD data. The system integrates multiple chromatin marks (DNase-seq, H3K27ac, H3K27me3, CTCF) to predict microexon accessibility. Furthermore, this project investigated the differential transcription factor binding motifs using the HOMER pipeline, identifying key TF associated with AD (Alzheimer's Disease) or NCI (no cognitive impairment) pathology.

This projects supporting findings
-AD-associated splicing changes likely occur through trans-acting factors (Sox family transcription factors) rather than cis-regulatory chromatin modifications, as many AD dysregulated microexon genes reside within the dynamically regulated SACS score (0/2) rather than highly accessible (SACS 4).

<img width="844" height="418" alt="image" src="https://github.com/user-attachments/assets/85c7ccc0-c907-4c06-8a29-6295dbeeb31a" />

### Key Findings 
(all performed using public ENCODE RUSH AD database)

SACS = 4 microexons show ~2-4x enrichment in chromatin accessibility vs. matched controls
Differentially spliced microexons (|dPSI| > 0.10) in AD include USP13, CERKL, STRADA, ANKRD30B, ABCC1, CACNA1H
Sox family TFs (Sox7, Sox17, Sox4) show 2-4x enrichment at AD-dysregulated microexons
HIF-1β motifs enriched 2.2x at accessible splice sites


## Data Requirements 

Input Data - *Not Included*
Microexon annotations - Microexonator
Genome annotation - GTF file
ChIP-seq/DNase-seq Data - available on [ENCODE](encodeproject.org/)

### Where to obtain public data?
ENCODE: https://www.encodeproject.org/ (brain DNase-seq, ChIP-seq)
RUSH AD Center: https://www.radc.rush.edu/ (AD brain epigenomics)



### Tools used in this analysis 
HOMER for motif analysis: http://homer.ucsd.edu/homer/
VASTTOOLS (optional) for differential splicing: https://github.com/vastgroup/vast-tools


## Quick Start

Step 1: Edit config.txt to correct input paths (microexon data, chip/dnase-seq data, and output results path)

Step 2: run the pipeline!

```
Rscript run_sacs_analysis.R
```


### Key Steps of Analysis
1. Data Preprocessing
   -filter microexons by percent spliced in (PSI) > 0.8, generate flanking genomic windows ~250nt each side of microexons with matched GC content.

2. Chromatin Mark Processing
   -process ChIP-seq peaks with **strand-aware** logic (forward & reverse)

3. DNase Processing
   -extract BigWig signal files & identify top 20% accessible microexons
   -find DNase hypersensitivity peak overlaps

4. SACS Integration
   -merge chromatin components and calculate SACS scores
   -validate with control regions and annotate with gene information

5. Differential Splicing
   -filter for dPSI > 0.10 between AD and control samples
   -test chromatin accessibility differences at differentially spliced MEs
   -Statistical testing with FDR correction

6.Motif Enrichment
  -prepare BED files for HOMER analysis to identify AD vs NCI increased chromatin accessibility regions


#### Contact
Author : Gabe Jandebeur
*Rotation Lab*: Wren Lab
Institution: University of Oklahoma Health Campus / Oklahoma Medical Research Foundation

Acknowledgements
-Wren Lab Members for their assistance and support
-ENCODE Consortium for the countless reference datasets
