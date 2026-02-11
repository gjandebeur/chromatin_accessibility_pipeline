

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
