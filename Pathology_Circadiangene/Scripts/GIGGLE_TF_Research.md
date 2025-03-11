Download data from Cistrom and create GIGGLE index with top 1K peaks in each bedfiles based on peak score.

go to https://db3.cistrome.org/browser/#/download download mm10 TF data

```
cd ~/Data/DataBase_project0126/Cistrome_ChIPTF/Liver_TF/202501_Cistrom/mm10_tranfac

for i in $(ls mm10_tranfac/*bed); do sample=`basename $i`; echo -e "sort -rn -k5 $i| head -1000 | sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > directory_Top1k_file/${sample}.gz" >>Cistrom_Top1kPeaks.sh; done
```

then use parallel to get the .gz file
```
parallel -j10 :::: Cistrom_Top1kPeaks.sh
```

#### Creat GIGGLE search index 


```
ulimit -Sn 16384
###Top 1K####
~/Software/giggle/bin/giggle index -i "directory_Top1k_file/*gz" -o Top1k_Giggle_index -f -s
```

1. use Step3 get the bed file of gain/lost Rhythmic genes to do giggle search

and the results are in correspdong directory
```
for i in $(ls *bed) 
do
 bgzip -c $i > ${i}.gz
done

for i in $(ls *bed.gz)
do
 file=$(basename $i .bed.gz)
 ~/Software/giggle/bin/giggle search -i ~/Data/DataBase_project0126/Cistrome_ChIPTF/Liver_TF/202501_Cistrom/Top1k_Giggle_index -q $i -s >${file}_GIGGLE_res_Top1K.xls
done

for i  in *.xls; do  libreoffice --headless --convert-to csv "$i" ; done
```

## Gro-seq motifs analysis####

```
#Lost Circadian in Case
cd /workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/LostRhyInCase/Lost_Enhancer.BedFiles

for i in $(ls *bed); do filename=`basename $i .bed`; bedtools intersect -a $i -b /workspace/rsrch1/tmp/DataBase_datatable0618/Gro-seq_Liver/Fangbin/Groseq_WTRhythmicPeaks_mm10.bed >./LostRhy_Overlap_Groseq/Groseq_${filename}.bed; done

cd LostRhy_Overlap_Groseq

 for i in $(ls *bed)
 do
  filename=`basename $i .bed`
  echo -e "findMotifsGenome.pl $i mm10 LostRhy_HomerMotifs_res/${filename} -mask" >>Lost_homerMotif.sh
  done



## for gain Circadian in Case same process with Lost
cd /workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/Gain_Enhancer.BedFiles

for i in $(ls *bed)
 do
 filename=`basename $i .bed`
 bedtools intersect -a $i -b /workspace/rsrch1/tmp/DataBase_datatable0618/Gro-seq_Liver/Fangbin/Groseq_WTRhythmicPeaks_mm10.bed >GainRhy_Overlap_GroSeq/Groseq_${filename}.bed
 done
 
 cd GainRhy_Overlap_GroSeq
 
for i in $(ls *bed)
do 
  filename=`basename $i .bed`
  echo -e "findMotifsGenome.pl $i mm10 GainRhy_HomerMotifs_res/${filename} -mask" >>Gain_homerMotif.sh
  done
```
 
