#!/bin/bash 

usage() {
	echo "Usage: $0 root_of_source_tree tarball_file"
	exit 1
}

[ $# = 2 ] || usage
src_dir=$1
tarball_file=$2

work_dir=`pwd`

# Normalized path of source root directory
cd "$src_dir" || { echo "Error: cannot access root of source tree" ; exit 1 ; }
normalized_src_dir=`pwd`
src_basename=`basename $normalized_src_dir`
src_parent=`dirname $normalized_src_dir`

# Tarball file contents selection 
sed_script=`mktemp` { echo "Error: unable to create temporary file" ; exit 1 ; }
cat > $sed_script <<\EOF
# Remove first . directory returned by find
1 d
# Exclude anything under ./build directory, but keep build directory
\:^\./build/: d
# Exclude any path containing a hidden file or directory: .svn or .project or ...
\:/\.: d
# Exclude any path starting with ./gtest
\:^\./gtest: d
# Replace leading . with basename of source root directory
EOF
echo 's/^\./'$src_basename'/g' >> $sed_script

tarball_contents=`mktemp` || { echo "Error: unable to create temporary file" ; exit 1 ; }
find | sed -f "$sed_script" > "$tarball_contents"

# Tarball file creation
cd "$work_dir"
tar --no-recursion --directory="$src_parent" --files-from="$tarball_contents" --create --verbose --gzip --file "$tarball_file" 
    





