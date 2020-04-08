SOURCE=$1
CONVERT_DATA_SCRIPT=$2
NAME_TEMPLATE=${3:-bac8_0.05_test}


for CUR_BAC in Bacillus Enterococcus Escherichia Lactobacillus Listeria Pseudomonas Salmonella Staphylococcus
do
  CUR_DIR=$SOURCE/$NAME_TEMPLATE${CUR_BAC}
  echo "Converting ${CUR_BAC}"
  python3 $CONVERT_DATA_SCRIPT $CUR_DIR/train.hdf5 $CUR_DIR --chunks 24000 --min-run 2
done
echo "Done"
