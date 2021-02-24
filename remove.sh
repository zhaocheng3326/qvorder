#!/bin/bash
ITEM=locfit
R CMD REMOVE -l /home/zhaoch/R/x86_64-pc-linux-gnu-library/3.0 $ITEM
sudo R CMD REMOVE -l /usr/local/lib/R/site-library $ITEM
sudo R CMD REMOVE -l /usr/lib/R/site-library  $ITEM
sudo R CMD REMOVE -l /usr/lib/R/library  $ITEM

