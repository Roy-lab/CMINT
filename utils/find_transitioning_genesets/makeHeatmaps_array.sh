if [ $# -lt 3 ]
then
	echo "Too few arguments $#"
	echo "Usage: ./makeHeatmaps.sh outputdir numberofclusters heatmap_script_path"
	exit
fi
OUTDIR=$1
FIGDIR=$OUTDIR/figs
if [ -d $FIGDIR ]
then
	echo "Please remove the figs directory in $OUTDIR"
	exit
fi
mkdir $FIGDIR
MAXCID=`expr $2 - 1`
HALFCID=`echo $MAXCID| awk '{printf("%d\n",$1/2)}' `
HEATMAPPER=$3
#First make the birdseye view of the figures
export MAXCID
export HALFCID
echo "cat $OUTDIR/ordered_clusterset_means.txt| $HEATMAPPER -vC=\"-2:(250,250,250);-1:(255,255,255);0:(90,180,172);$HALFCID:(245,245,245);$MAXCID:(216,179,101) -2:(0,0,255);0:(255,255,255);2:(255,0,0)\" -vStrokeC=\"-\" -vStrokeSC=\"black\" -vFontSize=8 > $FIGDIR/ordered_clusterset_means.svg"
cat $OUTDIR/ordered_clusterset_means.txt| $HEATMAPPER -vC="-2:(250,250,250);-1:(255,255,255);0:(90,180,172);$HALFCID:(245,245,245);$MAXCID:(216,179,101) -2:(0,0,255);0:(255,255,255);2:(255,0,0)" -vStrokeC="-" -vStrokeSC="black" -vFontSize=8 -vL="moduleIDs marklevel"> $FIGDIR/ordered_clusterset_means.svg
#Now make the figures for each of the clusters of reasonable size
for FNAME in `ls $OUTDIR/clusterset*.txt`
do
	#echo $FNAME
	FIGNAME=${FNAME/.txt/.svg}
	FIGNAME=${FIGNAME/$OUTDIR}
	FIGNAME=$FIGDIR/$FIGNAME
 	echo "cat $FNAME | $HEATMAPPER -vC=\"-2:(250,250,250);-1:(255,255,255);0:(90,180,172);$HALFCID:(245,245,245);$MAXCID:(216,179,101) -2:(0,0,255);0:(255,255,255);2:(255,0,0)\" -vStrokeC=\"-\" -vStrokeSC=\"black\" -vFontSize=6 > $FIGNAME"
 	cat $FNAME | $HEATMAPPER -vC="-2:(250,250,250);-1:(255,255,255);0:(90,180,172);$HALFCID:(245,245,245);$MAXCID:(216,179,101) -2:(0,0,255);0:(255,255,255);2:(255,0,0)" -vStrokeC="-" -vStrokeSC="black" -vFontSize=6 -vL="moduleID marklevel" > $FIGNAME
done

 echo "cat $OUTDIR/all_assign.txt | $HEATMAPPER -vC=\"-2:(250,250,250);-1:(255,255,255);0:(90,180,172);$HALFCID:(245,245,245);$MAXCID:(216,179,101) -2:(0,0,255);0:(255,255,255);2:(255,0,0)\" -vStrokeC=\"-\" -vStrokeSC=\"black\" -vFontSize=6 > $FIGDIR/all_assign.svg"
cat $OUTDIR/all_assign.txt | $HEATMAPPER -vC="-2:(250,250,250);-1:(255,255,255);0:(90,180,172);$HALFCID:(245,245,245);$MAXCID:(216,179,101) -2:(0,0,255);0:(255,255,255);2:(255,0,0)" -vStrokeC="-" -vStrokeSC="black" -vFontSize=6 -vL="moduleID marklevel" > $FIGDIR/all_assign.svg
