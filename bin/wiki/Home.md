# Project Cluster  

## 2018.06.0x  
* Now that thesis is out of the way, let's focus on this!  
* Using the top10pc pruning, have managed to get more clusters (which is amazing after all this time)!  

## 2018.03.27  
* I have been overly occupied with lab work until late evenings so have not had much time lately for this project. But Easter is coming up and I have no plans, so coding time!  
* Have read the papers suggested a while ago but have been struggling to implement the formulae onto the code. Let's see...  

## 2018.02.21  
* MCL takes in a similarity matrix so defo should do 1-val as otherwise it would have been distance.  
* The result is the same: just 1 cluster. I have no idea why.  

## 2018.02.20  
* Now that for HPA I have the Euclidean distance and the normalised version, I converted this into an MCL input.  
* But surprise, surprise, after running MCL on it I got only .... *1* cluster.  
* So I tried again, but this time I chose a random set of 2000 genes, redid the whole Euclidean distance calculation and normalisation and ran it with MCL, and nada.  
* Wondered if I should have done 1-val for each normalised so did that and ran MCL, och ingenting igen. Jag undrar varför?! :/  

## 2018.02.19  
* I think I got the covariance and therefore the inverse. Next step to get the partial correlation.  
* What troubled me is that calculating the inverse of such matrix gave me very different results with numpy, Matlab, and R.  
* The multiprocessing script was edited to save the calculations into a file straight away.  

### 2018.02.18  
* Enabled multiprocessing so processing the normalised matrix should be faster.  
* Partial correlation matrix, from Euclidean DM or from 19600 genes x 37 tissues?  

### 2018.02.16
* Script for matrix normalisation created, but perhaps needs optimisation
especially the way the nested list, which may be greedy with memory.  
* Script that converts a normalised matrix into a format accepted by MCL has
also been created.  
* Still looking into partial correlation matrix.  
* Time aim still remains.  

### 2018.02.14  
Create scripts that:  
* Normalises matrix,
* Create partial correlation matrix via inverse covariance matrix?,  
* Convert matrix to MCL input format,  
Time aim: mid v8.
