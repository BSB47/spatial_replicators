#! /bin/sh
# $Id: shuffle.sh 6392 2007-05-03 13:37:57Z ckarney $
#
# Shuffles the lines of a file (or standard input).  Requires that
#   RandomPermutation n
# produces a random permutation of the integers [0,n).
#
# Written by by Charles Karney <charles@karney.com> and licensed under
# the GPL.  For more information, see http://charles.karney.info/random/

usage="$0 [file]"

while getopts h c; do
    case $c in
	h ) echo "usage: $usage" 1>&2; exit 0;;
	* ) echo "usage: $usage" 1>&2; exit 1;;
    esac
done
shift `expr $OPTIND - 1`

case $# in
  0 ) FILE=;;
  1 ) FILE=$1;;
  * ) echo usage: $0 [file] 1>&2; exit 1;;
esac

TEMP=
trap 'trap "" 0; test "$TEMP" && rm -rf "$TEMP"; exit 1' 1 2 3 9 15
trap            'test "$TEMP" && rm -rf "$TEMP"'            0
TEMP=`mktemp -d ${TMPDIR:-/tmp}/shufXXXXXXXX`

if [ $? -ne 0 ]; then
   echo "$0: Can't create temp directory, exiting..." 1>&2
   exit 1
fi

if test -z "$FILE"; then
  FILE=$TEMP/in
  cat > "$FILE"
fi

l=`wc -l < "$FILE"`
# If file doesn't end with a newline, need to increment line count
b="`tail -1c "$FILE"`"
test "$b" = "`echo`" -o "$b" = "" || l=`expr $l + 1`
RandomPermutation $l > $TEMP/shuf
paste $TEMP/shuf "$FILE" | sort | cut -f2-
