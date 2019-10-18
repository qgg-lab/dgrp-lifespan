#! /bin/bash
# ============================================================

# =================
# = bwa alignment =
# =================

# get command line options
args=`getopt -o "b:,m:,x:,r:,c:,q:,n:,d" -l "bwa:,sam:,index:,ref:,conv:,fastq:,ncpu:,dir:" -- "$@"`
echo "command arguments given: $args"
eval set -- "$args"

# parse arguments
while true;
do
  case $1 in
    -b|--bwa)
       bwa=$2
       shift 2;;

    -m|--sam)
       sam=$2
       shift 2;;

    -x|--index)
       index=$2
       shift 2;;

    -r|--ref)
       ref=$2
       shift 2;;

    -c|--conv)
       conv=$2
       shift 2;;

    -q|--fastq)
       fastq=$2
       shift 2;;

    -n|--ncpu)
       ncpu=$2
       shift 2;;
    
    -d|--dir)
       dir=$2
       shift 2;;

    --)
       shift
       break;;
  esac
done

# check if files and programs exist
# bwa
bwacheck=`command -v $bwa`
samcheck=`command -v $sam`

if [[ $bwacheck ]]
then
  echo "using $bwacheck"
else
  echo "cannot find bwa"
  exit 1
fi

if [[ $samcheck ]]
then
  echo "using $samcheck"
else
  echo "cannot find samtools"
  exit 1
fi

# bwa index
if [[ -e $index.amb && -e $index.ann && -e $index.bwt && -e $index.pac && -e $index.sa ]]
then
  echo "using $index.*"
else
  echo "cannot find one of $index.amb .ann .bwt .pac .sa"
  exit 1
fi

# reference sequence
if [[ -e $ref.fasta && -e $ref.fasta.fai && -e $ref.dict ]]
then
  echo "using $ref.*"
else
  echo "cannot find one of $ref.fasta .fasta.fai .dict"
  exit 1
fi

# check if all fastq files exist
# the fastq file needs to follow this format
# rgid sm ln pu pair(?) illumina(I)/sanger(S) file1 file2
if [[ -e $fastq ]]
then
  echo "testing file existence for fastq files"
  while read line
  do
    pair=`echo $line | awk '{print $5}'`
    file1=`echo $line | awk '{print $7}'`
#    file1path=`eval realpath $file1`
    if [[ ! -e $file1 ]]
    then
      echo "cannot find $file1"
      exit 1
    fi
    if [[ $pair = "pair" ]]
    then
      file2=`echo $line | awk '{print $8}'`
#      file2path=`eval realpath $file2`
      if [[ ! -e $file2 ]]
      then
        echo "cannot find $file2"
        exit 1
      fi
    fi
  done < $fastq
else
  echo "cannot find the fastq input file list"
	exit 1
fi

# ============================================================
# main program

if [[ ! -e $dir ]]
then
  echo "$dir does not exist, making one"
  mkdir $dir
fi

while read line
do

  rgid=`echo $line | awk '{print $1}'`
  sm=`echo $line | awk '{print $2}'`
  lb=`echo $line | awk '{print $3}'`
  pu=`echo $line | awk '{print $4}'`
  pair=`echo $line | awk '{print $5}'`
  file1=`echo $line | awk '{print $7}'`
  file2=""
  qual=`echo $line | awk '{print $6}'`

  if [[ $qual == "I" ]]
  then
    gunzip -c $file1 | perl $conv illumina2std > $dir/file1.fq
  else
    gunzip -c $file1 > $dir/file1.fq
  fi

  if [[ $pair == "pair" ]]
  then
    file2=`echo $line | awk '{print $8}'`
    if [[ $qual == "I" ]]
    then
      gunzip -c $file2 | perl $conv illumina2std > $dir/file2.fq
    else
      gunzip -c $file2 > $dir/file2.fq
    fi
    $bwa mem -M -t $ncpu -R "@RG\tID:$rgid\tPL:illumina\tLB:$lb\tSM:$sm\tPU:$pu" $index $dir/file1.fq $dir/file2.fq 2> $dir/$rgid.mem.log | $sam view -bS -t $ref.fai - 2> $dir/$rgid.samview.log | $sam sort -m 20000000000 - $dir/$rgid.map 2> $dir/$rgid.samsort.log
    rm $dir/file1.fq $dir/file2.fq    
  else
    $bwa mem -t $ncpu -R "@RG\tID:$rgid\tPL:illumina\tLB:$lb\tSM:$sm\tPU:$pu" $index $dir/file1.fq 2> $dir/$rgid.mem.log | $sam view -bS -t $ref.fai - 2> $dir/$rgid.samview.log | $sam sort -m 20000000000 - $dir/$rgid.map 2> $dir/$rgid.samsort.log
    rm $dir/file1.fq
  fi

done < $fastq
