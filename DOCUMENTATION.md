# Documentation PBLMM

PBLMM is a cross-platform application for statistical analysis of proteomics datasets. It applies linear mixed effects models to assess differential protein expression, taking peptide level information into account.

## Usage

### Install

Run the provided installer for MacOS or Windows. You need to have python installed on your computer. The application contains a button to install necessary dependencies, so you don´t need to install them manually. 

The Windows Installer will install the application under /USER/AppData/Local/PBLMM. There you can find the Result files.

Dependencies:
    numpy
    pandas 
    statsmodels
    scipy
    matplotlib

### Input data

PBLMM is directly compatible with peptide or PSM files from ProteomeDiscoverer with its default naming of the columns. PBLMM is primarily label based, meaning the algorithm identifies the relevant columns from the column header names. If different names are used, the user can simply provide the custom names in the input form of the GUI making it compatible with any possible input data file, as long as they contain the relevant data columns.

The input file needs to be a tab-delimited text file in the wide format:
![Image](Untitled.jpg)

Note: The column names can differ but need to be specified. The input file can also contain other columns that are ignored by the algorithm.

The quantification columns should share a string that defines a quantification column, as seen as "Abundance" in the above example.

Additionally you need a tab delimited file containing the treatment specifications according to the order in your input file: E.g.

Control Control Control Treatment1  Treatment1  Treatment1 ...

Note: The design file should only contain the treatments, no column headers or anything else. The treatments should be separated by tabs.

### Calculations

PBLMM will automatically calculate the differential expression of all possible combinations of treatments.

In the implemented statistical model, the expression of each protein is separately modelled by a linear mixed effects model:
y(i)= ß0+βXi+ui+εi
Where yi denotes expression of peptide i, β_0 is the individual protein’s global intercept, βXi is the linear combination of indicator variables encoding categorical experimental conditions, ui is the additive random intercept of peptide i with ui ~ N(0,σ2u), and εi are residual errors with εi ~ N(0,σ2ε). Note that this collapses to ordinary linear regression when there are no multiple peptide measurements per protein.

### Results

The Results file will contain summed protein quantifications for each accession (to provide replicate quantifications for plotting). This quantification is NOT used for statistical analysis and is just used for visualization purposes. 

Additionally the Result file will contain columns for calculated P values, fold changes and q values (Benjamini Hochberg FDR). The naming of the columns is as follows:

p_value_Control_TreatmentX, fold_change_Control_TreatmentX etc.

The naming corresponds to the log2 ratio of the latter treatment versus the first treatment. For the above example, that would mean log2(TreatmentX/Control). A positive value interprets as more protein in TreatmentX compared to the Control.

The Result file is located in the installation folder of the application.