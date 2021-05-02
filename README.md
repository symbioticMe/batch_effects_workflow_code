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

1. Collins BC, Hunter CL, Liu Y, et al. Multi-laboratory assessment of reproducibility, qualitative and quantitative performance of SWATH-mass spectrometry. Nat Commun. 2017;8(1):291. Published 2017 Aug 21. doi:10.1038/s41467-017-00249-5
2. Sajic T, Liu Y, Arvaniti E, et al. Similarities and Differences of Blood N-Glycoproteins in Five Solid Carcinomas at Localized Clinical Stage Analyzed by SWATH-MS. Cell Rep. 2018;23(9):2819‐2831.e5. doi:10.1016/j.celrep.2018.04.114
3. Reference to Aging Mice study (if available);
4. Reference to PVCA/ Jeff Leek SVA code and paper;
5. ProjectTemplate URL.