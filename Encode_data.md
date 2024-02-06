# The file show us the process to download bed file for peaks from Encode website###

## search Encode database with filters on 02/05/2024
![image](https://github.com/ww1021pp/database_data/assets/60449311/dd25df85-f7ed-4bfe-a4ed-d20973113279)
![image](https://github.com/ww1021pp/database_data/assets/60449311/e12bdb7d-7df0-4f60-8f82-808f9f8f2bc8)

finnaly, I got 62 datasets for mouse TF files.

### Then click Downloads and fille download a file named files.txt in directory (/home/chunjie/Data/DataBase_project0126)
open terminal and use the code to generate the download command for each files in the files.txt
```
while read -r line
 do echo "wget $line" >>download_files.sh
 done <files.txt
gunzip *gz ## to uncompressed .gz files.
```
