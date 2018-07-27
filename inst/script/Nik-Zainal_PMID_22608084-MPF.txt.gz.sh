#!/bin/bash


# the MPF data file containing somatic SNVs for 21 tumors can be downloaded from:

wget -nd https://github.com/friend1ws/pmsignature/raw/devel/inst/extdata/Nik_Zainal_2012.mutationPositionFormat.txt.gz

# we opted for a shorter file name:

mv Nik_Zainal_2012.mutationPositionFormat.txt.gz Nik-Zainal_PMID_22608084-MPF.txt.gz
