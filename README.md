# SplNet
Splicing Network Reconstruction (Papasaikas, Tejedor Mol. Cell 2015)



-LABCHIPS.tar.gz: Contains the raw labhip data (.txt, .xlsx files) and the processed summary statistics (Z-scores, P-values).

-prepare_z_extra_data_2016.pl is a perl parser that takes as input a raw .txt LABCHIP file and calculates the summary statistics (this is a convoluted script, and any contact with it should be avoided if possible :-) ).

-SplNet.R contains main functions and core code for the network reconstruction starting from LabChip data
