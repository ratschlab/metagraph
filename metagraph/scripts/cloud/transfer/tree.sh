prefix=$1
for i in $(seq -f "%03g" 0 999); do
  for f in $prefix$i*gz; do
    ## Check if the glob gets expanded to existing files.
    ## If not, f here will be exactly the pattern above
    ## and the exists test will evaluate to false.
    [ -e "$f" ] && found=1 || found=0
    break
  done
  if (( found == 0)); then
    continue
  fi
  mkdir -p "$prefix$i"
  mv -- "$prefix$i"*.gz "$prefix$i/"
done