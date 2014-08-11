#!/bin/sh

# @SAKURA_LICENSE_HEADER_START @
# @SAKURA_LICENSE_HEADER_END @

header_pre="@SAKURA_LICENSE_HEADER"
header_start="_START@"
header_end="_END@"

replace() {
	file="$1"
	echo "Processing $file"
	prefix="`grep "$header_pre$header_start" $file | sed -n \"s/\\(.*\\)${header_pre}${header_start}.*/\\1/;1p\" `"
	sed -n "1,/$header_pre$header_start/p" $file > $tmp
	sed "s/^/$prefix/" $notice >> $tmp
	sed -n "/$header_pre$header_end/,\$p" $file >> $tmp
	mv $tmp $file
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