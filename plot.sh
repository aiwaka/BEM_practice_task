#!/bin/sh

SCRIPTS_DIR="./gnuplot_scripts"

while (( $# > 0 ))
do
  case $1 in
    -l)
      ls $SCRIPTS_DIR/
      ;;
    -*)
      echo "invalid option"
      exit 1
      ;;
    *)
      SCRIPT_NAME=$1

      gnuplot -persist "$SCRIPTS_DIR/$1"
      ;;
  esac
  shift
done
