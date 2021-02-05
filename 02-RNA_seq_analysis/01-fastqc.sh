

pdir="/data/NHLBI_BCB/Mehdi_P/RLee/RLee-RNA-seq-rlee-v2"

for i in `cat fastqs.txt` ; do
   sbatch --mem=2g  $pdir/02-fastqc/sbatch-fastqc.sh  ${i}  $pdir/01-fastqs  $pdir/02-fastqc
done
 



########### unzip them
for i in `cat $pdir/fastq_sids` ; do
   cd $pdir/02-fastqc/$i/
   unzip ${i}_R1_001_fastqc.zip
   unzip ${i}_R2_001_fastqc.zip
   cd -
done





# get the adaptors
fastqc_dir="$pdir/02-fastqc"
rm -f tmp1

for i in `cat $pdir/fastq_sids` ; do
  cd $fastqc_dir/$i/${i}_R1_001_fastqc/
  csplit fastqc_data.txt '/>>/' {*}
  seqn=`grep "Overrepresented sequences" xx* | head -1 | cut -f 1 -d':'`
  cat $seqn | sed '1,2d' | awk '{print $0}' | cut -f 1,4 | cut -f 1 -d',' >> ../../tmp1

  cd $fastqc_dir/$i/${i}_R2_001_fastqc/
  csplit fastqc_data.txt '/>>/' {*}
  seqn=`grep "Overrepresented sequences" xx* | head -1 | cut -f 1 -d':'`
  cat $seqn | sed '1,2d' | awk '{print $0}' | cut -f 1,4 | cut -f 1 -d',' >> ../../tmp1
done

cd $fastqc_dir
count=1
for i in `cat tmp11 | cut -f 1 | sort | uniq ` ; do
  echo "> overrep $count" >> tmp2
  echo $i >> tmp2
   (( count += 1 ))
done

rm -f tmp1



# append the default adaptors
cat  /data/NHLBI_BCB/bin/BBMap_36.11/bbmap/resources/adapters.fa tmp2 > adapters.fa

rm -f tmp2


