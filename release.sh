#!/bin/bash
#
# This script automates some release steps.
#
# - Updates lists of configuration parameters and diagnostics in the manual.
#
# - Updates lists of funding sources (in the manual and ACKNOWLEDGE.rst).
#
# - Updates PISM version in VERSION. This file is used to determine the PISM version
#   reported by "pismr -version" for PISM binaries built from a tarball (instead of a
#   working copy obtained using "git clone") and in documentation.
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

# update VERSION
echo ${version} > VERSION
git add VERSION

# update the list of diagnostics and configuration parameters
make -C doc/sphinx

# update funding sources
make -C doc
git add ACKNOWLEDGE.rst
git add doc/sphinx/funding.txt

# commit these changes
git commit -m "Update PISM version to ${version}"

# print changes made by this script
git --no-pager diff HEAD^

# tag this commit
git tag -a v${version} -m "The v${version} release. Please see CHANGES.rst for the list of changes."

set +x
echo "Now check if everything looks good and run:"
echo ""
echo "git push -u origin HEAD && git push origin tag v${version}"
echo ""
echo "Then create a release on GitHub."
