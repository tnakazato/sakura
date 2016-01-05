#!/bin/bash 

usage() {
	echo "Usage: $0 [options...]"
	exit 1
}


# TODO: Command line parsing
#[ $# = 2 ] || usage
#src_dir=$1
#tarball_file=$2

# Get the sources
#src_url='https://alma-dms.mtk.nao.ac.jp/svn/sakura/trunk/libsakura'
#src_url='https://dms.alma.nao.ac.jp/svn/sakura/trunk/libsakura'
src_url=$(cd $(dirname $0)/.. ; env LANG=en_US.UTF-8 svn info | grep -E "^URL:" | cut -d " " -f 2) 
project_name=`basename "$src_url"`

tmp_dir=`mktemp -d /tmp/sakura.XXXX`
src_dir="$tmp_dir/$project_name"

mkdir "$src_dir"
svn co "$src_url" "$src_dir" 1>/dev/null

# Script working directory
work_dir=`pwd`

cd "$src_dir" || { echo "Error: cannot access source directory:" ; echo "$1"  ; exit 1 ; }
src_dir_abs=`pwd`
src_parent_abs=`dirname $src_dir_abs`

# Tarball file
version_major=`egrep 'set.+libsakura_VERSION_MAJOR' CMakeLists.txt | head -n 1 | egrep --only-matching '[0-9]+'`
version_minor=`egrep 'set.+libsakura_VERSION_MAJOR' CMakeLists.txt | head -n 1 | egrep --only-matching '[0-9]+'`
svn_revision=`env LANG=en_US.UTF-8 svn info | grep 'Revision:' | egrep --only-matching '[0-9]+'`
release_version="${version_major}.${version_minor}.${svn_revision}"

tarball_dir="$work_dir"
tarball_name="libsakura-${release_version}.tar.gz"
tarball_file="${work_dir}/${tarball_name}"

# Tarball file contents selection 
sed_script="${src_parent_abs}/select_files.sed"
cat > $sed_script <<\EOF
# Remove first . directory returned by find
1 d
# Exclude anything under ./build directory, but keep build directory
\:^\./build/: d
# Exclude whole bin directory
\:^\./bin: d
# Exclude any path containing a hidden file or directory: .svn or .project or ...
\:/\.: d
# Exclude any path starting with ./gtest
\:^\./gtest: d
# Replace leading . with basename of source root directory
EOF
echo 's/^\./'$project_name'/g' >> $sed_script

release_contents="${src_parent_abs}/release_contents.txt"
find . | sed -f "$sed_script" > "$release_contents"

# Tarball file creation
cd "$work_dir"
rm -f "$tarball_file"
tar --no-recursion --directory="$src_parent_abs" --files-from="$release_contents" --create --verbose --gzip --file "$tarball_file" 1>/dev/null

echo 'Created tarball file:'
echo "$tarball_file"

# Tarball check
n_mismatches=`tar tfz "$tarball_file" | sed 's:/$::g' | diff - ${release_contents} | wc -l`

# TODO: Rpm file creation
    





