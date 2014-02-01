# PathoVar
The Variant Annotation Tool for Pathoscope pipeline

## Why This? Aren't there a dozen SNP Annotation Platforms out there?
Most SNP Tools are made for people and model organisms, not pathogenic bacteria, viruses, or such. Many tools for discovering drug involvement are oriented on starting with the drug class in mind. 

## TODO
Things that ought to be done when we have a spare minute

- Allow multiple programmatic schemes for filtering variants 
- Make notes of redundancy between PATRIC, VIPR, and GenBank
- Take the maximum number of threads as a parameter, and use multiprocessing/subprocess to separate SNP calling and annotation parsing to be on *max_threads* processes at a given time.