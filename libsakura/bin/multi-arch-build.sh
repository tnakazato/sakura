#!/bin/sh

ProjDirPath="`dirname $0`"
ProjDirPath="`(cd \"$ProjDirPath\"/..; pwd)`"
BuildDirPath="$PWD" # "${ProjDirPath}/build"
prefix_base="${BuildDirPath}/installed"

CMAKE_BUILD_TYPE=RelWithDebInfo
BUILD_DOC=ON	# Use Doxygen to create the HTML based API documentation
SCALAR=OFF	# Disable auto-vectorization by compiler
DOC_LANG=English	# "Language for the document generation
VECT_VERB=""	# vectorization verbose compile option, such as -fopt-info-vec-optimized -fopt-info-vec-missed -ftree-vectorizer-verbose=2
PROF_MODE=""	# -fprofile-generate or -fprofile-use -fprofile-correction for GCC

test -z "$CC" && CC=gcc
test -z "$CXX" && CXX=g++

case $# in
	0)
		;;
	*)
		prefix_base="$1"
		shift
		;;
esac

suffix_map () {
	case $1 in
		NATIVE)
			echo "default"
			;;
		*)
			echo "$1" | tr '[:upper:]' '[:lower:]'
			;;
	esac
	return 0
}

archs="SSE4 AVX AVX2"

for_each_arch () {
	for arch in $archs; do
		ARCH="$arch"
		PREFIX="${prefix_base}-`suffix_map $arch`"
		"$@" || return 1
	done
}

run_mkdir () {
	mkdir -p "${BuildDirPath}/$ARCH"
}

run_cmake () {
	cd "${BuildDirPath}/$ARCH"
	env CC="$CC" CXX="$CXX" cmake \
		-D CMAKE_MODULE_PATH="${ProjDirPath}/cmake-modules" \
		chdir "$BuildDirPath" \
		cmake -D CMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
		-D CMAKE_INSTALL_PREFIX:PATH="$PREFIX" \
		-D SIMD_ARCH:STRING=$ARCH \
		-D SCALAR:BOOL=$SCALAR \
		-D BUILD_DOC:BOOL=$BUILD_DOC \
		-D DOC_LANG:STRING="$DOC_LANG" \
		-D VECT_VERB:STRING="$VECT_VERB" \
		-D PROF_MODE:STRING="$PROF_MODE" \
		"$ProjDirPath"
}

run_command () {
	cd "${BuildDirPath}/$ARCH"
	"$@"
}

for_each_arch run_mkdir || exit 1
for_each_arch run_cmake || exit 1
for_each_arch run_command make || exit 1
for_each_arch run_command make test || exit 1
for_each_arch run_command make apidoc || exit 1
for_each_arch run_command make install || exit 1

