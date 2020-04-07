#!/bin/bash
BAC8_TO_READ=$1
SOURCE=$2
DEST=$3
PERCENTAGE=$4

mkdir DEST

for CUR_BAC in Bacillus Enterococcus Escherichia Lactobacillus Listeria Pseudomonas Salmonella Staphylococcus
do
  echo "Making $CUR_BAC"
  CUR_DIR=$DEST/bac8_${PERCENTAGE}_test$CUR_BAC
  mkdir -p $CUR_DIR
  python3 ./bin/create_dataset.py bac8 split $BAC8_TO_READ $SOURCE $CUR_BAC -p $PERCENTAGE --output $CUR_DIR/train.hdf5 --test_output $CUR_DIR/test.hdf5
done

echo "Done!"
#CUR_BAC=Bacillus
#echo "Making $CUR_BAC"
#python3 bac8 split $BAC8_TO_READ $SOURCE $CUR_BAC -p $PERCENTAGE --output DEST/bac8_$PERCENTAGE_test
#
#echo "Making Enterococcus"
#python3 bac8 split $BAC8_TO_READ $SOURCE Enterococcus -p $PERCENTAGE
#
#echo "Making Escherichia"
#python3 bac8 split $BAC8_TO_READ $SOURCE Escherichia -p $PERCENTAGE
#
#echo "Making Lactobacillus"
#python3 bac8 split $BAC8_TO_READ $SOURCE Lactobacillus -p $PERCENTAGE
#
#echo "Making Listeria"
#python3 bac8 split $BAC8_TO_READ $SOURCE Listeria -p $PERCENTAGE
#
#echo "Making Pseudomonas"
#python3 bac8 split $BAC8_TO_READ $SOURCE Pseudomonas -p $PERCENTAGE
#
#echo "Making  Salmonella"
#python3 bac8 split $BAC8_TO_READ $SOURCE Salmonella -p $PERCENTAGE
#
#echo "Making Staphylococcus"
#python3 bac8 split $BAC8_TO_READ $SOURCE Staphylococcus -p $PERCENTAGE

