# BIOS694_Final_Project
This is the final project of course BIOS694 from McGill, which aims to critically summarize, review and replicate the results in the paper.
In BIOS694_Final_Project repository, there is a dataset called file2754771351f4.arff which is the Yeast dataset used to do the replication.
I have also attached the original paper here for easier access.
Also, there is a PDF called results_comparison that has three tables which aims compare the multi-label classification results (Table3) in the paper and my own replication results (respectively use 5/15 positive examples from each label as training data).
Finally, the replication R file is uploaded.

For reproducibility, please download the .arff dataset and replace your own path to the dataset and run the whole code, then it should give you the table 2 results shown in the results_comparison PDF. Changing code line 143 and 144, by replacing number 5 by 15 and rerunning the code, you will get replication results as using 15 positive examples from each label, as shown in table 3 in results_comparison PDF.
