# BATSallTime_LTT
Want BV calculation from CTD casts through all time at BATS; use 2db resolution data\
Krista Longnecker, 4 April 2025\
Krista Longnecker, 10 April 2025

### Some details:
This was done in MATLAB in order to best use Ruth Curry's code organizing the BATS CTD data and for making calculations of BV frequency (using the SEAWATER library from CSIRO). 

I saved a CSV with the results, but that file is too large for GitHub so it is sitting in the BIOS-SCOPE Google Drive. If you are reading this here, chances are pretty good you know how to find the file.

After talking with Ruth, plot this in log10 space and filter the BV data (using get_bvfilt.m)

The primary m-file is written in MATLAB: ([getBATSdata_LTT_KLv2.m](https://github.com/redbluewater/BATSallTime_LTT/blob/main/getBATSdata_LTT_KLv2.m)), and you can follow the link if you want to see what I did. That m-file calls a few other m-files, all of which are here in GitHub.

Here is a JPG with the results - two subplots (top is using the 2db CTD data and bottom is the plot from the data in the Excel file):
<img src = "https://github.com/redbluewater/BATSallTime_LTT/blob/main/BATS_bvfrq_allYears.jpg" width="100%" height="100%">
