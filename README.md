# Asite_LP_method
The following text describes the process to run the Linear Programming method to identify the A-site positions within ribosome protected fragments. 
The analysis is carried out in two steps. 
Step 1: An alignment file is processed to create read count files according to fragment length. These files will be used in Step 2 to run the Linear Programming algorithm. Genes are selected based on filtering criteria specified by the user.
Step 2: Using the read count files created in Step 1, Linear Programming algorithm is run and the output is a A-site offset table specific to fragment size and frame. 
...
