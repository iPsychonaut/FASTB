cd "${1:-.}"  # directory holding g5.fasta g50.fasta g200.fasta
FB=${FB:-fastb}  # path to the fastb command
best() {  # best of 3: wall seconds and peak RSS KB for a command writing to $OUT
  local w=999 r=0
  for i in 1 2 3; do
    /usr/bin/time -f "%e %M" -o t.txt bash -c "$1" 2>/dev/null
    read e m < t.txt
    if (( $(echo "$e < $w" | bc) )); then w=$e; r=$m; fi
  done
  echo "$w $r"
}
echo "| file | tool | threads | bytes | encode s | encode RSS MB | decode s | decode RSS MB |"
echo "|---|---|---|---|---|---|---|---|"
for f in g5 g50 g200; do
  for p in 1 4 8 16; do
    read te re < <(best "pigz -6 -k -f -p $p $f.fasta")
    read td rd < <(best "pigz -dc -p $p $f.fasta.gz > out.fa")
    echo "| $f | pigz -6 | $p | $(stat -c %s $f.fasta.gz) | $te | $((re/1024)) | $td | $((rd/1024)) |"
  done
  read te re < <(best "$FB encode $f.fasta -o $f.fastb")
  for p in 1 4 8 16; do
    read td rd < <(best "$FB cat -p $p $f.fastb > out.fb")
    echo "| $f | fastb | $p | $(stat -c %s $f.fastb) | $te | $((re/1024)) | $td | $((rd/1024)) |"
  done

  $FB cat -p 8 $f.fastb > out.fb; cmp out.fa out.fb && echo "| $f | (fastb cat output identical to pigz -dc) | | | | | | |"
done
rm -f t.txt out.fa out.fb
echo; nproc; lscpu | grep "Model name"; python3 --version; ~/fastb-venv/bin/python -c "import numpy;print('numpy',numpy.__version__)"
