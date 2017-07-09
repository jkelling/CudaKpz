#!/bin/bash

if [ -z "$TINYMT32DC" ]; then
	TINYMT32DC="tinymt32dc"
fi
TINYMT32DCFILE="tinymt32_param.h"
TMPFILE="/tmp/$TINYMT32DCFILE.tmp"

if [ -z "$1" ]; then
	N="$(( 32*2048*16 ))"
else
	N=$1
fi

echo -n \
"/* tinyMT parameters generated for CudaKmc
 * using $TINYMT32DC
 * " > $TINYMT32DCFILE
date >> $TINYMT32DCFILE
echo -e ' */\n' >> $TINYMT32DCFILE

echo -e '#ifndef KMC_TMT_PARAM_H\n#define KMC_TMT_PARAM_H\n' >> $TINYMT32DCFILE
echo -e "const int KMC_N_TMT_PARAM_SETS = $N;\n" >> $TINYMT32DCFILE

I=0;
rm -f "$TMPFILE"
while [ "$I" -lt "$N" ]; do
	$TINYMT32DC $I | tail -n1 >> "$TMPFILE"
	I=$(( $I + 1 ))
done

echo "static const unsigned int tinyMTmat1[KMC_N_TMT_PARAM_SETS] = {" >> $TINYMT32DCFILE
cut -d , -f 4 $TMPFILE | sed 's/^/\t0x/;s/$/,/' >> $TINYMT32DCFILE
echo -e '};\n' >> $TINYMT32DCFILE
echo "static const unsigned int tinyMTmat2[KMC_N_TMT_PARAM_SETS] = {" >> $TINYMT32DCFILE
cut -d , -f 5 $TMPFILE | sed 's/^/\t0x/;s/$/,/' >> $TINYMT32DCFILE
echo -e '};\n' >> $TINYMT32DCFILE
echo "static const unsigned int tinyMTtmat[KMC_N_TMT_PARAM_SETS] = {" >> $TINYMT32DCFILE
cut -d , -f 6 $TMPFILE | sed 's/^/\t0x/;s/$/,/' >> $TINYMT32DCFILE
echo -e '};\n' >> $TINYMT32DCFILE

echo '#endif' >> $TINYMT32DCFILE

rm "$TMPFILE"
