# CMIP6_ICM-CSIC
Where three pioneers explain how they discovered data from the CMIP6 world and, in their infinite kindness, gifted this knowledge to the ICM.

1/ Data selection was made there: https://aims2.llnl.gov/search
- Activity ID: ScenarioMIP (model predictions for different ICPP scenarios of climate change)
- Resolution: 100 km
- Variant Label: r1i1p1f1 (this is based on an email from Pablo Ortega)
- Frequency: yr (for most variable) and mon (for temperature/thetao)
- Variable ID: https://redmine.dkrz.de/attachments/download/5773/AR6%20WG1%20priority%20variables.xlsx
- Result Type: Originals Only (avoid duplicates)
=> display 20 files per page (for more files the wget script won't be generated), add to cart, download wget script, copy the http links to new text file, do this for all files for each variable
=> I then subset the http links based on the presence of year 2100 in the filename (this is consistent across files)

2) On ada, I then run the download with command:
```for i in `awk '{print $2}' thetao_2100.txt`; do wget $i; done```
4) the processing of the data is then performed with the script "cmip6_process.R".

