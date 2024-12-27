Note: Used code from course materials.

The code first downloads the dataset and untars the folder.

It then reads the three files from which data would be required. 
data_mrna_seq_v2_rsem.txt
data_clinical_patient.txt
data_cna.txt

We then match the ids to create the metadata using the CNA level of ERBB2. 

After this, the data is normalized using DESeq2.

We then obtain the top 10 differentially expressed genes. 

After that, we create the prinicpal component analysis (PCA) plot.

Then a pathway enrichment analysis is performed. 

We then create the heatmap.
