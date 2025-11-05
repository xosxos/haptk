#!/bin/bash

# Pipe in the list of files to read

IGV=$1
GENOME=$2
NAME=$3
COORDS=$4

if [[ -z $NAME || -z $COORDS || -z $IGV || -z $GENOME ]]; then
  echo 'no name, coordinates or directory given: check args from source code'
  exit 1
fi

cat << EOT > ${NAME}_batch.txt
new

genome $GENOME

preference SAM.SHOW_SOFT_CLIPPED true
preference SAM.MAX_VISIBLE_RANGE 700

$(while read line; do echo "load $line"; done;)

snapshotDirectory ${NAME}
goto ${COORDS}
sort strand
sort position
snapshot ${NAME}.png
exit
EOT

sh $IGV --batch=${NAME}_batch.txt;

rm ${NAME}_batch.txt

for file in ${NAME}/${NAME}.png; do convert -quality 100 -density 300 -trim $file ${file%%.*}.pdf; done
