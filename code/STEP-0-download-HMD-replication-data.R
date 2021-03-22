############################################################
# Download and unzip the directory containing the HMD
# data used for the analysis. This replication dataset
# is stored on the HMD site, and is available to any
# registered HMD user. 
# 
# After (free) registration, insert your HMD username and
# password into the code below.
# 
# This program
#  1. downloads the zipped data from HMD
#  2. creates a new directory called ../HMD-input-data/  
#     (a sibling of the code/ directory that contains this script)
#  3. unzips the downloaded file to create three subdirectories
#     in the new ../HMD-input-data/ directory, called 
#         ../HMD-input-data/Deaths_1x10
#         ../HMD-input-data/Exposures_1x10
#         ../HMD-input-data/fltper_1x10
#     each of which contains a large number of HMD text files
#     that will be used by other programs    
############################################################

## RUN THIS PROGRAM WITH THE code/ directory (where this script lives)
## AS THE WORKING DIRECTORY

library(R.utils)

my_HMD_username = '...insert your HMD username here...'
my_HMD_password = '...insert your HMD password here...'

HMDfile='https://www.mortality.org/hmd/REPLICATION/schmertmann_replication_testfile.zip'

downloadFile(url=HMDfile, 
             filename='tmp.zip',
             username= my_HMD_username, 
             password= my_HMD_password)

dir.create('../HMD-input-data')

unzip(zipfile='tmp.zip', exdir='../HMD-input-data')

file.remove('./tmp.zip')