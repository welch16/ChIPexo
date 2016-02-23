#!/bin/sh

## performs chip exo analysis with peakzilla

codedir=$HOME/Desktop/Docs/Code/lib/peakzilla
indir=/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo/bedfiles
outdir=/p/keles/ChIPexo/volume6/K12/other_methods/peakzilla

python $codedir/peakzilla.py $indir/edsn1311_Sig70.sort.bed > $outdir/edsn1311_Sig70_peaks.tsv
python $codedir/peakzilla.py $indir/edsn1314_Sig70.sort.bed > $outdir/edsn1314_Sig70_peaks.tsv
python $codedir/peakzilla.py $indir/edsn1317_Sig70.sort.bed > $outdir/edsn1317_Sig70_peaks.tsv
python $codedir/peakzilla.py $indir/edsn1320_Sig70.sort.bed > $outdir/edsn1320_Sig70_peaks.tsv
python $codedir/peakzilla.py $indir/edsn931_Sig70.sort.bed > $outdir/edsn931_Sig70_peaks.tsv
python $codedir/peakzilla.py $indir/edsn933_Sig70.sort.bed > $outdir/edsn933_Sig70_peaks.tsv
