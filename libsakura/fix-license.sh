#!/bin/sh

# @SAKURA_LICENSE_HEADER_START @
# @SAKURA_LICENSE_HEADER_END @

header_pre="@SAKURA_LICENSE_HEADER"
header_start="_START@"
header_end="_END@"

replace() {
	file="$1"
	echo -n "Processing $file: "
	(echo | cat - $file) | sed "1,/$header_pre$header_start/d" | grep --quiet "$header_pre$header_end" || { echo "Missing end marker. Skipped." ; return; }
	prefix="`grep "$header_pre$header_start" $file | sed -n \"s/\\(.*\\)${header_pre}${header_start}.*/\\1/;1p\" `"
	(echo | cat - $file) | sed -n "1,/$header_pre$header_start/p" | sed '1d'  > $tmp && sed "s/^/$prefix/" $notice >> $tmp && sed -n "/$header_pre$header_end/,\$p" $file >> $tmp && mv $tmp $file && echo "Done."
}

usage() {
	echo "Usage: $0 license_description_file root_of_source_tree"
	exit 1
}

[ $# = 2 ] || usage

notice="$1"
root="$2"
tmp="$root/tmp.$$"

find . -type f | grep -v '/.svn/' | while read file; do
	file $file | grep --quiet text && grep --quiet "$header_pre$header_start" $file && replace $file
done
