SOURCE=$1
NAME_TEMPLATE=${$2:bac8_0.05_test}

for CUR_BAC in Bacillus Enterococcus Escherichia Lactobacillus Listeria Pseudomonas Salmonella Staphylococcus
do
  CUR_DIR=$SOURCE/$NAME_TEMPLATE${CUR_BAC}
  ./scripts/convert-data $CUR_DIR/train.hdf5 $CUR_DIR --chunks 24000 --min-run 2
done
