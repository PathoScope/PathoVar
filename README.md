# PathoVar
The Variant Annotation Tool for Pathoscope pipeline

## Why This? Aren't there a dozen SNP Annotation Platforms out there?
Most SNP Tools are made for people and model organisms, not pathogenic bacteria, viruses, or such. 

## TODO
Things that ought to be done when we have a spare minute

. `locate_varianty.py` currently reparses the GenBank XML files every time. We could save time on re-annotation if we save the resulting objects to disk. Maybe translate the final form used to annotate to dictionaries and use `json` or `cpickle` to save it to disk.
. Change variant frequency to `Minor Allele Frequency` or `MAF`.
. Allow multiple programmatic schemes for filtering variants. 
. Make notes of redundancy between PATRIC, VIPR, and GenBank. Run a batch study to check this out.
. Integrate Immune Epitope Database Web Service.
. Shrink the `locate_variant` GENE tag.
. Collate the genes with variants into a summary file.
. Figure out if a visualization is low-hanging fruit.