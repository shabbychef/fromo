#!/bin/bash -e
#
# Created: 2016.03.27
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
# Comments: Steven E. Pav


COMPRESS_LEVEL=0
VERBOSE=""

HELP_STRING="$0 [-v] [-Z zlev]";
while getopts "Z:vh" flag
do
	case "$flag" in
		(h) echo "$HELP_STRING foo";
			  exit 0;;
		(Z) COMPRESS_LEVEL=$OPTARG;;
		(v) VERBOSE='--verbose';;
	esac
done

shift $((OPTIND-1))

# c.f. http://kvz.io/blog/2013/11/21/bash-best-practices/
# and  http://fahdshariff.blogspot.com/2013/10/shell-scripting-best-practices.html
# and  http://redsymbol.net/articles/bash-exit-traps/

function finish() {
	git stash pop -q
}

# http://codeinthehole.com/writing/tips-for-using-a-Git-pre-commit-hook/
git stash -q --keep-index

trap finish EXIT
 
# blame Hadley: http://r-pkgs.had.co.nz/release.html
if [[ nodist/README.Rmd -nt README.md ]]; then
  echo "README.md is out of date; please re-knit README.Rmd"
  exit 1
fi 

# make attributes
r -l Rcpp -e 'compileAttributes(".")'
exit $?

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=sh:ft=sh:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
