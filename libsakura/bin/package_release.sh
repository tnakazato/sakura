#!/bin/bash -e
# -e  Exit immediately if a command exits with a non-zero status.

# TODO: enable parallel rpm build on ana02 and an03
# -> need to sort out / simplify tmp dirs and tmp files :
#    only 1 top tmp dir in current dir, output files must also go there  

set -o pipefail
#set -x

export PATH=/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin

# Script exit status message
status() {
    if [[ $? -ne 0 ]] ; then
        echo "Status: FAIL"
    fi
}

trap status EXIT

# Constants
FALSE=0
TRUE=1

# Options default values
# ----------------------
# Sakura sources / svn repository url
src_url_default='https://dms.alma.nao.ac.jp/svn/sakura/trunk/libsakura'
src_url=${src_url_default}
# Sakura sources / svn revision
svn_revision_default='HEAD'
svn_revision=${svn_revision_default}
# Sakura sources / svn branch
svn_branch_default='trunk'
svn_branch=${svn_branch_default}
svn_branches_prefix='branches'
# Binary RPM / enable/disable generation
do_rpm=${TRUE}
# External resources directory
rpm_resources_dir_default="/nfsstore/sakura_casa/rpm"
rpm_resources_dir=${rpm_resources_dir_default}
# Binary RPM / package version
rpm_package_version_default='1'
rpm_package_version=${rpm_package_version_default}
# Binary RPM / file name / include or not svn revision 
rpm_short_name=${TRUE} # svn revision not included in rpm file name 
# Binary RPM / sakura version / include or not svn revision
rpm_short_version=${TRUE}
# Binary RPM / include legacy libraries for CASA
rpm_legacy_libs=${FALSE}

# Usage
usage() {
	echo "Usage: $(basename $0) [options...]"
cat <<USAGE_END

Synopsys:
    Package a release of Sakura into tar source archive and rpm binary archive. 

Options:
    -h|--help ..................: display this message and exit

    Release selection
    -u|--src_url <url> .........: URL of Sakura repository
         default ...............: ${src_url_default}
    -b|--svn_branch <branch>....: trunk or branch name under 
                                  ${src_url_default/trunk/libsakura/${svn_branches_prefix}}                        
         default ...............: ${svn_branch_default}
    -r|--svn_revision <rev>.....: svn revision number of release to package

    RPM-related options
    --no_rpm|--norpm ...........: do not create rpm archive
    --rpm_resources_dir ........: path to directory storing resources required
                                  to create the rpm file.
         default ...............: ${rpm_resources_dir_default}
    -p|--rpm_package_version <p>: package version
         default ...............: ${rpm_package_version_default}
    -ln|--rpm_long_name ........: include svn revision number in rpm file name
    -lv|--rpm_long_version .....: include svn revision number in rpm archive
    --rpm_legacy_libs ..........: include Sakura legacy libraries in rpm archive
                                  applies only to rhel5 and rhel6 rpms targeting CASA version < 5
                                  rhel6 rpm targeting CASA version >= 5 does not include Sakura legacy libraries
                                  rhel7 rpm does not include Sakura legacy libraries

Examples:
---- Package HEAD revision of trunk, tar + rpm:
$(basename $0)
---- Package HEAD revision of trunk, tar only:
$(basename $0)  --norpm
---- Package revision 1880 of trunk, tar + rpm:
$(basename $0) -r 1880
---- Package revision 1770 of branch 'release-casa-4_5', tar + rpm:
$(basename $0) -b release-casa-4_5 -r 1770
USAGE_END
}


# Command line parsing
while [[ $# > 0 ]] ; do
    key="$1"
    case $key in
        -h|--help)
        usage
        exit 0
        ;;
        -u|--src_url)
        shift
        src_url=$1
        ;;
        -b|--svn_branch)
        shift
        svn_branch=$1
        ;;
        -r|--svn_revision)
        shift
        svn_revision=$1
        ;;
        -p|--rpm_package_version)
        shift
        rpm_package_version=$1
        ;;
        --rpm_resources_dir)
        shift
        rpm_resources_dir=$1
        ;;
        -ln|--rpm_long_name)
        rpm_short_name=${FALSE}
        ;;
        -lv|--rpm_long_version)
        rpm_short_version=${FALSE}
        ;;
        --rpm_legacy_libs)
        rpm_legacy_libs=${TRUE}
        ;;
        --no_rpm|--norpm)
        do_rpm=${FALSE}
        ;;
        *) # Unknown option
        usage
        echo "Error: illegal option: ${key}"
        exit 1
        ;;
    esac
    shift
done

# Script working directory
work_dir=`pwd`

# Checkout sakura sources
project_name="libsakura"

# svn branch URL
if [[ ${svn_branch} != ${svn_branch_default} ]]; then
    src_url=${src_url/${svn_branch_default}/"${svn_branches_prefix}/${svn_branch}"} 
fi

tmp_dir=$(mktemp -d -p /tmp ${project_name}_src.XXXX)
src_dir="$tmp_dir/$project_name"
mkdir "$src_dir"

# Note: svn co command must be run from a clean, 'svn free' directory
cd "$src_dir" || { echo "Error: cannot access source directory:" ; echo "$src_dir"  ; exit 1 ; }
svn co --quiet --revision ${svn_revision} "$src_url" "$src_dir" 
svn_info_log="$tmp_dir/svn_info.log"
svn info $PWD > ${svn_info_log}

src_dir_abs=`pwd`
src_parent_abs=`dirname $src_dir_abs`

# 1. -------- Tarball file
version_major=`grep -E 'set.+libsakura_VERSION_MAJOR' CMakeLists.txt | head --lines=1 | grep -E --only-matching '[0-9]+' | head --lines=1`
version_minor=`grep -E 'set.+libsakura_VERSION_MINOR' CMakeLists.txt | head --lines=1 | grep -E --only-matching '[0-9]+' | head --lines=1`
# By default, svn_revision='HEAD', must replace with actual numeric value 
svn_revision=`env LANG=en_US.UTF-8 svn info | grep 'Revision:' | grep -E --only-matching '[0-9]+' | head --lines=1`
release_version="${version_major}.${version_minor}.${svn_revision}"

tarball_dir="$work_dir"
tarball_name="libsakura-${release_version}.tar.gz"
tarball_file="${tarball_dir}/${tarball_name}"

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
find | sed -f "$sed_script" > "$release_contents"

# Tarball file creation
cd "$work_dir"
rm -f "$tarball_file"
tar --no-recursion --directory="$src_parent_abs" --files-from="$release_contents" --create --verbose --gzip --file "$tarball_file" 1>/dev/null

if [[ ${do_rpm} -eq ${FALSE} ]]; then
    echo 'The following release ...'
    cat ${svn_info_log}
    echo
    echo '... has been packaged into:'
    echo "$tarball_file"
    exit 0
fi

# 2. -------- Binary RPM file creation
rpmbuild_dir=$(mktemp -d -p /tmp ${project_name}_rpmbuild.XXXX)

# Google test sources
gtest_url="https://github.com/google/googletest/archive/release-1.7.0.zip"
gtest_name="googletest-release-1.7.0.zip"
gtest_file="${rpm_resources_dir}/${gtest_name}"
gtest_expansion_name=$(basename ${gtest_name} ".zip") 

if [[ ! -f ${gtest_file} ]] ; then
   echo "File not found:"
   echo ${gtest_file}
   echo "==> 1. Download from:"
   echo ${gtest_url}
   echo "==> 2. Copy to:"
   echo ${gtest_file}
   echo "ERROR: Google Testing Framework zip file missing"
   exit 1
fi

# RHEL version
rhel_version_major=$(lsb_release -r | grep -E --only-matching '[0-9]+' | head --lines=1)

# Legacy libsakura libraries for CASA
# -- RHEL7: never include legacy libsakura libraries, enforce rpm_legacy_libs=false
if [[ "${rhel_version_major}" = "7" ]] ; then
    rpm_legacy_libs=${FALSE}
fi
# -- RHEL5 and RHEL6: include legacy libsakura libraries based on ${rpm_legacy_libs}
if [[ ${rpm_legacy_libs} -eq ${TRUE} ]] ; then
	export SAKURA_LEGACY_LIBS=$(find ${rpm_resources_dir}/libsakura-{0.1.1352,1.1.1690,2.0.1882}/el${rhel_version_major} -type f)
fi

# RPM Macros
my_rpm_macros="${rpmbuild_dir}/rpmmacros"
rpm_define_options=""
sign_option=""
# -- Default install directory
case "${rhel_version_major}" in
    5)
    # CASA version < 5
    casa_install_prefix=/usr/lib64/casa/01
    ;;
    6)
    # CASA version <5: casa_install_prefix=/usr/lib64/casa/01
    # CASA version >=  5
    casa_install_prefix=/opt/casa/02
    ;;
    7)
    # CASA version <5: casa_install_prefix=/opt/casa/01
    # CASA version >=5:
    casa_install_prefix=/opt/casa/02
    ;;
    *) # Unknown RHEL distribution
    echo "Error: RHEL version major must be 5,6 or 7 but is set to:" 
    echo "${rhel_version_major}"
    exit 1
    ;; 
esac

rpmbuild_macros() {
cat > ${my_rpm_macros} <<EOF_RPM_MACROS
%_topdir ${rpmbuild_dir}
%_prefix ${casa_install_prefix}
EOF_RPM_MACROS
    # %dist tag and gpg signature on RHEL5   
    if [ ${rhel_version_major} = "5" ] ; then
        # %dist tag is not defined on RHEL5
        echo >> ${my_rpm_macros} "%dist .el${rhel_version_major}" 
        # gpg signature
        echo >> ${my_rpm_macros} "%_signature gpg" 
        echo >> ${my_rpm_macros} "%_gpg_name NAOJ Sakura" 
        # rpmbuild command / sign option
        sign_option="--sign"
        # PGP private key directory
        export GNUPGHOME=${rpm_resources_dir}/gnupg
    fi
    # Equivalent command line options 
    rpm_define_options=$(sed "s:%\(.*\):--define='\1':g" ${my_rpm_macros})
}

# Directory structure required by rpmbuild command
rpmbuild_prepare(){
    rm -rf ${rpmbuild_dir}
    mkdir -p  ${rpmbuild_dir}/{BUILD,SOURCES,SPECS,BUILDROOT,RPMS}
    cp ${tarball_file} ${rpmbuild_dir}/SOURCES
    cp ${gtest_file}   ${rpmbuild_dir}/SOURCES
}

# RPM Spec file generation
spec_file="${project_name}.el${rhel_version_major}.spec"
# -- RPATH control
case "${rhel_version_major}" in
    5)
    need_rpath=${FALSE}
    ;;
    6|7)
    need_rpath=${TRUE}
    ;;
    *) # Unknown RHEL distribution
    echo "Error: RHEL version major must be 5,6 or 7 but is set to:" 
    echo "${rhel_version_major}"
    exit 1
    ;; 
esac
create_spec_file() {
rpm_release_version=${release_version}
if [[ ${rpm_short_version} -eq ${TRUE} ]]; then
    rpm_release_version=${version_major}.${version_minor} # hide .${svn_revision}
fi
cat > $spec_file <<PREAMBLE_VERSION
Name:           libsakura 
Version:        ${rpm_release_version}
PREAMBLE_VERSION
echo >> $spec_file \
"Release:        ${rpm_package_version}"'%{?dist}'
cat >> $spec_file <<'PREAMBLE_URL'
Summary:        High performance library for astronomical data analysis.
Vendor:         National Astronomical Observatory of Japan
Group:          Development/Libraries
License:        LGPLv3 or later
URL:            http://alma-intweb.mtk.nao.ac.jp/~sakura/
PREAMBLE_URL
cat >> $spec_file <<PREAMBLE_SOURCES
Source0:        ${tarball_name}
Source1:        ${gtest_name}
PREAMBLE_SOURCES
cat >> $spec_file <<'PREAMBLE_PREFIX'
BuildRoot:      %(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXX)
Prefix:         %{_prefix}
PREAMBLE_PREFIX

cat >> $spec_file <<'GTEST_BEGIN'
# Fix for rpm generation error on RHEL5: 
# disable debuginfo package generation
%define debug_package %{nil}

%description
Sakura library is optimized to use SIMD instructions to maximize utilization of CPU performances. 
%prep
%setup -q -n %{name}
%setup -q -T -D -a 1 -n %{name}
GTEST_BEGIN
cat >> $spec_file <<GTEST_END
ln -s ${gtest_expansion_name} gtest
GTEST_END
cat >> $spec_file <<'CMAKE_COMMAND_BEGIN'
%build
cd build && rm -rf *
prefix_root=%{_prefix}
prefix_no_root=${prefix_root:1}
sse4_install_prefix=${RPM_BUILD_ROOT}/${prefix_no_root}/lib/%{name}/default
CMAKE_COMMAND_BEGIN
# RPATH / LDFLAGS
ld_flags=""
if [[ ${need_rpath} -eq ${TRUE} ]] ; then
    ld_flags='LDFLAGS=-Wl,--enable-new-dtags'
fi
# CMake command / common options
# Note: need to escape last backslash when here-document is parsed    
cat >> $spec_file <<CMAKE_COMMAND
${ld_flags} cmake \\
CMAKE_COMMAND
cat >> $spec_file <<'CMAKE_OPTION_DOC'
  -D CMAKE_MODULE_PATH=$(dirname $PWD)/cmake-modules \
  -D CMAKE_INSTALL_PREFIX=${sse4_install_prefix} \
  -D CMAKE_BUILD_TYPE=Release \
  -D SIMD_ARCH=SSE4 \
  -D BUILD_DOC:BOOL=ON \
CMAKE_OPTION_DOC
# CMake command / rpath options
if [[ ${need_rpath} -eq ${TRUE} ]] ; then
	# Careful: editor must not insert spaces instead of tabs for this to work
	# First character of lines below must be tabs
	cat >> $spec_file <<-'CMAKE_RPATH_OPTIONS'
	  -D CMAKE_SKIP_BUILD_RPATH=FALSE \
	  -D CMAKE_BUILD_WITH_INSTALL_RPATH=FALSE \
	  -D CMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE \
	  -D CMAKE_INSTALL_RPATH='/cas-8362/some/non/existent/path/for/later/rpath/replacement' \
	CMAKE_RPATH_OPTIONS
fi
# CMake command / common options / path to sources   
cat >> $spec_file <<'CMAKE_COMMAND_END'
  ..
CMAKE_COMMAND_END
# Make + install
cat >> $spec_file <<'SAKURA_LEGACY_BEGIN'
make
make apidoc
%install
rm -rf ${RPM_BUILD_ROOT}
cd build
make install
# Dirty: move share doc directory to the right place: just under casa_prefix
prefix_root=%{_prefix}
prefix_no_root=${prefix_root:1}
mv ${RPM_BUILD_ROOT}/${prefix_no_root}/lib/%{name}/default/share ${RPM_BUILD_ROOT}/${prefix_no_root}
# Recursively delete empty directories (Doxygen: CREATE_SUBDIRS=YES)
# from the installation tree
find ${RPM_BUILD_ROOT}/${prefix_no_root}/share -type d -empty -delete
SAKURA_LEGACY_BEGIN
# Sakura legacy libraries
if [[ ${rpm_legacy_libs} -eq ${TRUE} ]] ; then
	cat >> $spec_file <<-'PREAMBLE_END'
		cp ${SAKURA_LEGACY_LIBS} ${RPM_BUILD_ROOT}/${prefix_no_root}/lib/%{name}/default/lib
	PREAMBLE_END
fi
# Installed files to be packaged
cat >> $spec_file <<'INSTALLED_FILES_SELECTION'
#
# -------- Spec file %files stage:
# select which installed files must be packaged
%files
%defattr(-,root,root)
INSTALLED_FILES_SELECTION
#%{_prefix}: does not work, so write the install prefix from bash
cat >> $spec_file <<EOF_SPEC_FILE
${casa_install_prefix}
EOF_SPEC_FILE
}

# RPM build environment
set_build_env() {
    build_host=$(hostname --short)
    case ${build_host} in
        ana01)
	source /opt/rh/devtoolset-4/enable
        ;;
        ana02)
	source /opt/rh/devtoolset-4/enable
        ;;
        ana03)
        source /opt/rh/devtoolset-2/enable
        ;;
        *) # Unknown build host
        echo "Error: unknown build host: ${build_host}"
        exit 1
        ;; 
    esac
}

# RPM generation
rpmbuild_prepare
rpmbuild_macros
set_build_env
create_spec_file

rm -f ${HOME}/.rpmmacros
echo 'rpm build: start ...'
echo "  . build log file: rpm_build.$(hostname --short).log"
# TODO: if rpmbuild fails FAIL message is not displayed on standard output
# since it is redirected to rpm_build.log
eval rpmbuild -v ${sign_option} ${rpm_define_options} -bb ${spec_file} 1>rpm_build.$(hostname --short).log 2>&1 
echo 'rpm build: done'

# Move rpm file to working directory and rename if needed
rpm_file=$(find ${rpmbuild_dir}/RPMS -type f)
rpm_name=$(basename ${rpm_file})
if [[ ${rpm_short_name} -eq ${TRUE} ]]; then
    # Remove .<svn_revision> from rpm_name
    rpm_name=${rpm_name/\.${svn_revision}/} 
fi
mv ${rpm_file} ${work_dir}/${rpm_name}

echo
echo 'The following release: ...'
echo '-------------------------------------------------------------------------'
sed '/^$/ d' ${svn_info_log}
echo '-------------------------------------------------------------------------'
echo '... has been packaged into:'
echo "$tarball_file"
echo "${work_dir}/${rpm_name}"







    





