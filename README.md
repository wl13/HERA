# HERA-nonhpc
# Introduction
This is a modified version of HERA pipeline in order to run locally instead of a job scheduling system. 

Please refer to the original site https://github.com/liangclab/HERA for any detailes about the HERA pipeline.

The modified version has not been extensively tested so use it with cautious.


# Installation

All same as original post expect download and use HERA-nonhpc instead of original scripts

The running of HERA requires a few other software programs. 
1. Downloading and installing bwa-0.7.10
   
   git clone https://github.com/lh3/bwa.git  
   cd bwa; make.
2. Downloading and installing DALIGNER

   https://github.com/thegenemyers/DALIGNER
3. Downloading and installing DAZZ_DB

   https://github.com/thegenemyers/DAZZ_DB
   
 # The example of running HERA
 
 Assume the working directory is "~/"
 
 1. Download and unzip HERA, and then generate a folder named "HERA-nonhpc-master". 
    
    unzip HERA-nonhpc-master.zip
    
 2. Set scripts to executable.
    
    cd HERA-nonhpc-master/ && chmod 750 *
    
 3. Generate the folder name "Test" and files for testing.
 
    cd ../ && mkdir Test
    
    cd Test/ && cp ../HERAv1.0-master/pipeline.sh ./
    
    mv ../HERA-nonhpc-master/*.fasta ./
    
 4. Modify the configuration part in pipeline.sh.

    set genome_seq=~/Test/Test_Genome.fasta
    
    set Corrected_Pacbio=~/Test/Test_CorrectedPacbio.fasta
    
    set Working_Script=~/HERA-nonhpc-master/
    
    set DAZZ_DB=~/DAZZ_DB-master/
    
    set DALIGNER=~/DALIGNER-master/
    
    set MinPathNum=3
    
      
 5. Run the pipeline-nonhpc.sh
 
    sh pipeline-nonhpc.sh
    


   
# Changes

The modified pipeline only use the PERL scripts as the binary versions works abnormally at my end. A few perl scripts are edited to ensure the pipeline could be run properly.




# Citing HERA

Du, H., Liang, C. (2018). Assembly of chromosome-scale contigs by efficiently resolving repetitive sequences with long reads. bioRxiv    doi: https://doi.org/10.1101/345983

