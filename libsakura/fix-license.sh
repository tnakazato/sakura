#!/bin/sh

# @SAKURA_LICENSE_HEADER_START@
# Copyright (C) 2013-2015
# National Astronomical Observatory of Japan
# 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
# 
# This file is part of Sakura.
# 
# Sakura is free software: you can redistribute it and/or modify it under 
# the terms of the GNU Lesser General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your 
# option) any later version.
# 
# Sakura is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License 
# along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
# @SAKURA_LICENSE_HEADER_END@

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
