# TopAffy
TopAffy is an approach to infer binding preferences for transcription factors without related binding data. TopAffy constructs a  graph to represent DBD protein sequences  and learns binding preferences of neighbouring amino acid pairs to every possible k-mer DNA sequence from PBM data. Given DBD protein sequences and corresponding PBM binding data for some of these sequences, TopAffy explains the binding data using a linear model that adds the binding preferences of neighbouring amino acids weighted by the "purity" of those neighbouring amino acids in TF binding domains.

TopAffy consists of two Python3 scripts: [TopAffyTrain.py] (TopAffyTrain.py) and [TopAffyPredict.py] (TopAffyPredict.py).

To run first TopAffyTrain.py, one needs to modify the following lines in the script:
* Line 162: variable SeqFile is a csv containing all the sequences to train the model (i.e., for which there are PBM binding data)
* Line 163: variable HMMFile is a csv containing all the sequences used to build the graph
* Line 164: variable profileFile  is a csv of all the 8-mer profiles used to train the model. TFs need to be in the same order as SeqFile
* Line 170: outFile is the name attached to all files created by TopAffy

TopAffyTrain.py will create the emission table named “emission” + outFile + “.csv”, and the sequence graph named outfile + ".csv". 

To predict 8-mer binding profile for some TFs, one has to run TopAffyPredict.py. To run TopAffyPredict.py one needs to modify the following lines in the script:
* Line 69: variable SeqFile is a csv containing all the sequences for which to predict  binding preferences 
* Line 70: filename which should be the same as outFile above
* Line 71: variable hmmSeq is a csv containing all the sequences used to build the graph

This will output the predicted 8-mer profiles for the sequences in SeqFile in that order.


If you use TopAffy, please cite:

Zier-Vogel R., Peña-Castillo L. (2020) TopAffy: Inferring binding preferences for transcription factors without related binding data.
