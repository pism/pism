#!/bin/bash
#
# This script automates some release steps.
#
# - Updates lists of configuration parameters and diagnostics in the manual.
#
# - Updates lists of funding sources (in the manual and ACKNOWLEDGE.rst).
#
# - Sets Pism_BRANCH to "stable" in CMakeLists.txt. This variable is used to distinguish
#   PISM releases ("stable") from the code in development ("development").
#
# - Updates PISM version in CMake/PISM_CMake_macros.cmake. This variable is used to
#   determine the PISM version reported by "pismr -version" for PISM binaries built from a
#   tarball (instead of a working copy obtained using "git clone").
#
# - Updates PISM version in the manual.
#
# - Tags the new release.
#
# - Prints commands needed to complete the release.
#
# Run this script *after* merging new code into main.

set -u
set -e
set -x

version=${1:?Usage: \"$0 X.Y.Z\" where X.Y.Z is the new PISM version.}

# update the list of diagnostics and configuration parameters
make -C doc/sphinx

# update funding sources
make -C doc

# set Pism_BRANCH to "stable"
sed -Ei 's/(Pism_BRANCH)[^)]+/\1 "stable"/' CMakeLists.txt
git add CMakeLists.txt

# update Pism_VERSION
sed -Ei "s/(set +\(Pism_VERSION)[^)]+/\1 \"v${version}\"/" CMake/PISM_CMake_macros.cmake
git add CMake/PISM_CMake_macros.cmake

# update version in the manual
sed -Ei "s/^(version|release).+/\1 = '${version}'/" doc/sphinx/conf.py
git add doc/sphinx/conf.py

# commit these changes
git commit -m "Update PISM version to ${version}"

# print changes made by this script
git --no-pager diff HEAD^

# tag this commit
git tag -a v${version} -m "The v${version} release. Please see CHANGES.rst for the list of changes."

set +x
echo "Now check if everything looks good and run:"
echo ""
echo "git push -u origin HEAD"
echo "git push --tags"
echo ""
echo "Then create a release on GitHub."
