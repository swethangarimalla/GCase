# GCase
This script was created to analyze the GCase Enyzme assay output. 

Input: Machine (which machine?? which software version) output for the BCA & GC enzyme assays and the Tissue Inventory 

Output: 
(1) A prism file per tissue containing the mean units per milligram, sample number, and dose group
(2) An excel file with the sample number, tissue, mean units per milligram, and dose group
(3) An excel file with the raw data, calculated columns (including Grubbs test data), and associated tissue inventory data (dose group)


Pipeline:

(1) Load and combine the BCA and GC files as excel workbooks into one dataframe since there can be multiple BCA files and GC files - one for each plate run. *Query user for folder name.* 
  (1.1) Open each excel file for both BCA and GC files and read Dilution sheets in 
        each workbook
  (1.2) For each assay (BCA and GC), bind the data from all of the worksheets into one         dataframe. The result at this step will be two data frames, each with data   
        combined across across workbooks - one for BCA and one for GC
  (1.3) Add "B_" to BCA column names in the BCA dataframe.
  (1.4) Add "G_" to GC column names in the GC dataframe.
  (1.5) Remove all empty rows from both data frames
  (1.6) Make sure B_Abs and G_RFU are both recognized as numeric for further
        calculations
  (1.7) Perform the Grubbs test for the replicates for B_Abs data
  (1.8) Perform the Grubbs test for the replicates for the G_RFU data
  (1.9) Merge the two dataframes with the remaining values
(2) Load and process tissue inventory file *Tissue Inventory should be in the user input folder*
  (2.1) Read in Tissue Inventory file. The name of this file should contain "Tissue Manifest"
  (2.2) Collect "Sample.ID", "Group", "Gender", and "TissueType" columns from this file.
  (2.3) Merge the BCA-GC dataframe (1.9) with the tissue inventory (2.2) by Sample.ID.
(3) Output barplots and files
  (3.1) Plot data as bar plots with mean as the level of the bar and sem as error bars.
  (3.2) Create a prism file per tissue. Each file contains the mean units/mg, Sample ID, and dose group
  (3.3) Create an R file with all calculated columns including Grubbs data
  (3.4) Create an R file containing only the Sample.ID, mean units/mg, tissue, and dose group. 
