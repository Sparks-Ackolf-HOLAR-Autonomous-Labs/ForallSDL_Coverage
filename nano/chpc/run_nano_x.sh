#!/bin/bash
# Nano campaign, one composition: reuse the party-4a 24^3 VCA h5 if it
# exists, else build it (x = 0.1, 0.3), then phononinfo at the five
# pre-registered temperatures. Outputs: ~/mat_theory/nano/modes/x<X>/
set -euo pipefail
X="$1"; GRID=24
module load almabte 2>/dev/null
B=~/mat_theory/party4_bte
OUT=~/mat_theory/nano/modes/x$X
mkdir -p "$OUT"; cd "$OUT"
H5=$(ls $B/work_x$X/h5/*_${GRID}_${GRID}_${GRID}.h5 2>/dev/null | head -1 || true)
if [ -z "$H5" ]; then
  SI=$(python3 -c "print(round(1-$X, 6))")
  cat > vca.xml <<XML
<alloy>
  <materials_repository root_directory="$B/repo"/>
  <compound name="Si" mixfraction="$SI"/>
  <compound name="Ge" mixfraction="$X"/>
  <gridDensity A="$GRID" B="$GRID" C="$GRID"/>
  <target directory="h5"/>
  <overwrite/>
</alloy>
XML
  VCAbuilder vca.xml > vca.log 2>&1
  H5=$(ls h5/*_${GRID}_${GRID}_${GRID}.h5)
fi
echo "$H5" > h5_used
for T in 300 500 700 900 1100; do phononinfo "$H5" $T > pi_$T.log 2>&1; done
mv "$(dirname "$H5")"/*.phononinfo . 2>/dev/null || true
mv ./*K.phononinfo . 2>/dev/null || true
cp $B/work_x$X/out/*.Tsweep . 2>/dev/null || true
sha256sum *.phononinfo $(ls *.Tsweep 2>/dev/null) > sha256sums
echo 0 > exit_status
