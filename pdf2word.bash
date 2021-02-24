#!/bin/sh
mkdir tmp
cp test.part3.pdf tmp
cd tmp
pdftoppm *  -r 150 ocrbook
for i in *.ppm; do convert "$i" "`basename "$i" .ppm`.tif"; done
for i in *.tif; do tesseract "$i" "`basename "$i" .tif`" ; done
for i in *.txt; do cat $i >> ${name}.txt; echo "[pagebreak]" >> pdf-ocr-output.txt; done
mv pdf-ocr-output.txt ..
rm *ppm *tif
cd ..
rmdir tmp
