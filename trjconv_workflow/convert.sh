set -e
GMX=gmx_2019.3

TPR=$PWD/run000/md1.tpr
INDEX=$PWD/run000/index.ndx

MDXTC=extend.part0002
FOLDER=run001

FIT_GROUP=DPPC_POPC
EXTRACT_GROUP=DPPC_POPC


SKIP=1

# Define names of output files:
OUT_SKIP=$MDXTC-skip-$SKIP.xtc
OUT_WHOLE=whole.xtc
OUT_NO_JUMP=no-jump.xtc
OUT_SAME_BOX=in-box.xtc
OUT_FIT=fit.xtc
OUT_EXTRACT=$EXTRACT_GROUP.xtc

cd $FOLDER

# Remove files if they already exists:
for file in $OUT_SKIP $OUT_WHOLE $OUT_NO_JUMP $OUT_SAME_BOX $OUT_FIT $OUT_CLL $OUT_EXTRACT
do
    if [ -f $file ] ; then
        echo Removing file: $file
        rm $file
    fi
done

CURRENT=$MDXTC.xtc


if [ "$SKIP" -gt "1" ]; then
    # Down-sample the source XTC:
    echo System | $GMX trjconv -f $MDXTC.xtc -o $OUT_SKIP -s $TPR -skip $SKIP
    CURRENT=$OUT_SKIP
fi


# Make whole:
echo System | $GMX trjconv -f $CURRENT -o $OUT_WHOLE -s $TPR -pbc whole
CURRENT=$OUT_WHOLE

# Remove jumps:
echo System | $GMX trjconv -f $CURRENT -o $OUT_NO_JUMP -s $TPR -pbc nojump
CURRENT=$OUT_NO_JUMP

# Fit to initial
echo $FIT_GROUP System | $GMX trjconv -f $CURRENT -o $OUT_FIT -s $TPR -fit rot+trans -n $INDEX
CURRENT=$OUT_FIT

# Move everything into the same box:
echo System | $GMX trjconv -f $CURRENT -o $OUT_SAME_BOX -s $TPR -ur compact -pbc mol
CURRENT=$OUT_SAME_BOX

# Extract a single group:
echo $EXTRACT_GROUP | $GMX trjconv -f $CURRENT -o $OUT_EXTRACT -s $TPR -n $INDEX
CURRENT=$OUT_EXTRACT

cd ..
