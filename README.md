This is the case study application of batch effects analysis protocol, 
mostly implemented using `proBatch` package.

Here we investigate batch effects and their correction in three datasets:
1. InterLab study [1];
2. PanCancer study [2];
3. AgingMice study [3].

Each dataset analysis is stored in its own folder, and is designed to be launched from that folder. 

In each folder, you will find:
1. Data (empty, files to be downloded from PRIDE);
2. Analysis;
3. Plots (empty, to be filled during execution of the analysis).

Functions, common for all three datasets, are stored in `lib` folder.

For more details on the analysis and structure of respective folder, we recommend to revert to the respective README within each folder.

The structure of this code and the folder are inspired by fSVA manuscript code [4] and 
`ProjectTemplate` R package [5].

1. Reference to InterLab study;
2. Reference to PanCancer study;
3. Reference to Aging Mice study (if available);
4. Reference to PVCA/ Jeff Leek SVA code and paper;
5. ProjectTemplate URL.