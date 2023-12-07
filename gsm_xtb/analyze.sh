rm -rf analysis
mkdir analysis
if [ -d $1 ]; then
    mawk -f scripts/extract_path.awk $1/stringfile.xyz0001
    mawk -f scripts/analyze_path.awk $1/stringfile.xyz0001
    cp scripts/align.py analysis
    cd analysis
    echo -ne Aligning...
    python align.py
    cd ..
    echo done.
else
    echo $1 not found
fi
