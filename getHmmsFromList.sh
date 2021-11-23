while getopts S:L:M:i:d: flag
do
    case "${flag}" in
        S) mfccsTrainList=${OPTARG};;
        L) labInputDir=${OPTARG};;
        M) hmmOutputDir=${OPTARG};;
		i) prototypeHMM=${OPTARG};;
		d) wordList=${OPTARG};;
    esac
done

WORDS=$(cat $wordList)

for WORD in $WORDS
do
	HInit -S $mfccsTrainList -l $WORD -L $labInputDir -M $hmmOutputDir -o $WORD -T 1 $prototypeHMM
done
	