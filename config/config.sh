echo -n "Enter Deliverables Directory PATH which must contain a folder named symlink : "
read DIR

FASTQ_DIR=${DIR}/symlink
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Directory '$FASTQ_DIR' not found."
    exit 1
fi

if [ -f ${DIR}/config.yaml ]; then
    echo "Error: config.yaml already exists in ${DIR}. Plese remove and execute again"
    exit 1
else
    echo "config.yaml not found. Proceeding with the script..."
    # Continue with your script logic here
fi

cd ${DIR}
wget https://raw.githubusercontent.com/sarangian/grenepipe/master/config/config.yaml

if [ -f ${DIR}/samples.tsv ]; then
    echo "Error: samples.tsv already exists in ${DIR}. Plese remove and execute again"
    exit 1
else
    echo "samples.tsv not found. Proceeding with the script..."
    # Continue with your script logic here
fi


echo "sample\tunit\tplatform\tfq1\tfq2" > ${DIR}/samples.tsv
for i in `ls ${FASTQ_DIR} | grep -E 'R1'`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R1_001.fastq.gz: "
read R1_extn

for i in `ls ${FASTQ_DIR} | grep -E 'R2'`;do echo $i;done
echo -n "Please look into the read names and provide the extension after sample name Example: _R2_001.fastq.gz: "
read R2_extn

for i in `ls ${FASTQ_DIR} | grep ${R1_extn}`;do 
sample=`basename $i ${R1_extn}`
FQ1=${FASTQ_DIR}/${sample}${R1_extn}
FQ2=${FASTQ_DIR}/${sample}${R2_extn}
echo "${sample}\t1\tILLUMINA\t${FQ1}\t${FQ2}" >> "${DIR}/samples.tsv"
done
