# PathoVar
The Variant Annotation Tool for Pathoscope pipeline

## Why This? Aren't there a dozen SNP Annotation Platforms out there?
Most SNP Tools are made for people and model organisms, not pathogenic bacteria, viruses, or such. Many tools for discovering drug involvement are oriented on starting with the drug class in mind. 

## Integrating Annotations
 - DrugBank
 - Comprehensive Antibiotic Resistance Database

## Functions
 - Call variants using `samtools`
 - Calculate per nucleotide coverage for each genome called against
 - Map variants into genes in genomes called against
 - Annotate each gene with protein homology associated with drug interaction
 - Annotate the effect of each variant using snpEff

