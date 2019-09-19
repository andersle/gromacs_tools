set -e
GMX=gmx_2019.3

# Look for all .edr files:
find . -name '*.edr' | while read EDR_FILE; do
    echo "Processing energy file: $EDR_FILE"
    out_dir=$(dirname -- "${EDR_FILE}")
    base_name=$(basename -- "${EDR_FILE}")
    file_name="${base_name%.*}"
    OUT=$out_dir/box-$file_name.xvg
    echo $OUT
    $GMX energy -f $EDR_FILE -o $OUT << EOF
Box-X
Box-Y
Box-Z
EOF
done
