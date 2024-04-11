#!/usr/bin/env bash
#
# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
#
# Copyright 2024- The Mem3DG Authors
# and the project initiators Cuncheng Zhu, Christopher T. Lee, and
# Padmini Rangamani.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Please help us support Mem3DG development, by citing the research
# papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/.

# This script runs copyright header checks on modified files and
# reports/applies the necessary changes. It is shamelessly borrowed
# from the GROMACS project.
#
# See `copyright.sh -h` for a brief usage description.

# Parse command-line arguments
function usage() {
    echo "usage: copyright.sh [-f|--force] [--rev=REV]"
    echo "           [--copyright=<cmode>]"
    echo "           [--warnings=<file>] [<action>]"
    echo "<action>: (check*|diff|update)[-(index|workdir*)] (*=default)"
    echo "<cmode>:  off|add|update*|replace|full"
}


### SHORTCIRCUIT LOGIC
# Switch to the root of the source tree and check the config file
srcdir=`git rev-parse --show-toplevel`
pushd $srcdir >/dev/null
admin_dir=$srcdir/admin

# Actual processing starts: create a temporary directory
tmpdir=`mktemp -d -t mem3dg_copyright.XXXXXX`

git ls-tree --full-tree --name-only -r HEAD > $tmpdir/files
cat $tmpdir/files | git check-attr --stdin filter | sed -e 's/.*: filter: //' | paste $tmpdir/files - | grep -E 'copyright$' > $tmpdir/filtered; cat $tmpdir/filtered
grep -E 'copyright$' < $tmpdir/filtered| cut -f1 > $tmpdir/filelist_copyright

if ! $admin_dir/copyright.py -F $tmpdir/filelist_copyright --add-missing --update-header --replace-header
    then
        echo "Copyright checking failed!"
        rm -rf $tmpdir
        exit 2
    fi

rm -rf $tmpdir
exit 0


#### THE REST IS IGNORED
action="check-workdir"
declare -a diffargs
baserev="origin/main"
force=
copyright_mode=update
warning_file=
for arg in "$@" ; do
    if [[ "$arg" == "check-index" || "$arg" == "check-workdir" || \
          "$arg" == "diff-index" || "$arg" == "diff-workdir" || \
          "$arg" == "update-index" || "$arg" == "update-workdir" ]]
    then
        action=$arg
    elif [[ "$arg" == "check" || "$arg" == "diff" || "$arg" == "update" ]] ; then
        action=$arg-workdir
    elif [[ "$action" == diff-* ]] ; then
        diffargs+=("$arg")
    elif [[ "$arg" == --copyright=* ]] ; then
        copyright_mode=${arg#--copyright=}
    elif [[ "$arg" == "-f" || "$arg" == "--force" ]] ; then
        force=1
    elif [[ "$arg" == --rev=* ]] ; then
        baserev=${arg#--rev=}
    elif [[ "$arg" == --warnings=* ]] ; then
        warning_file=${arg#--warnings=}
    elif [[ "$arg" == "-h" || "$arg" == "--help" ]] ; then
        usage
        exit 0
    else
        echo "Unknown option: $arg"
        echo
        usage
        exit 2
    fi
done

# Switch to the root of the source tree and check the config file
srcdir=`git rev-parse --show-toplevel`
pushd $srcdir >/dev/null
admin_dir=$srcdir/admin

# Actual processing starts: create a temporary directory
tmpdir=`mktemp -d -t mem3dg_copyright.XXXXXX`

# Produce a list of changed files
# Only include files that have proper filter set in .gitattributes
internal_diff_args=
if [[ $action == *-index ]]
then
    internal_diff_args="--cached"
fi
git diff-index $internal_diff_args --diff-filter=ACMR $baserev >$tmpdir/difflist
cut -f2 <$tmpdir/difflist | \
    git check-attr --stdin filter | \
    sed -e 's/.*: filter: //' | \
    paste $tmpdir/difflist - | \
    grep -E 'copyright$' >$tmpdir/filtered
cut -f2 <$tmpdir/filtered >$tmpdir/filelist_all
grep -E 'copyright$' <$tmpdir/filtered | \
    cut -f2 >$tmpdir/filelist_copyright
git diff-files --name-only | grep -Ff $tmpdir/filelist_all >$tmpdir/localmods

# Extract changed files to a temporary directory
mkdir $tmpdir/org
if [[ $action == *-index ]] ; then
    git checkout-index --prefix=$tmpdir/org/ --stdin <$tmpdir/filelist_all
else
    rsync --files-from=$tmpdir/filelist_all $srcdir $tmpdir/org
fi
# Duplicate the original files to a separate directory, where all changes will
# be made.
cp -r $tmpdir/org $tmpdir/new

# Create output file for what was done (in case no messages get written)
touch $tmpdir/messages

# Run code formatting on the temporary directory
cd $tmpdir/new

# Update the copyright headers using the requested mode
if [[ $copyright_mode != "off" ]] ; then
    cpscript_args="--update-year"
    case "$copyright_mode" in
        year)
            ;;
        add)
            cpscript_args+=" --add-missing"
            ;;
        update)
            cpscript_args+=" --add-missing --update-header"
            ;;
        replace)
            cpscript_args+=" --replace-header"
            ;;
        full)
            cpscript_args+=" --add-missing --update-header --replace-header"
            ;;
        *)
            echo "Unknown copyright mode: $copyright_mode"
            exit 2
    esac
    if [[ $action == check-* ]] ; then
        cpscript_args+=" --check"
    fi
    # TODO: Probably better to invoke python explicitly through a customizable
    # variable.
    if ! $admin_dir/copyright.py -F $tmpdir/filelist_copyright $cpscript_args >>$tmpdir/messages
    then
        echo "Copyright checking failed!"
        rm -rf $tmpdir
        exit 2
    fi
fi

cd $tmpdir

# If a diff was requested, show it and we are done
if [[ $action == diff-* ]] ; then
    git diff --no-index --no-prefix "${diffargs[@]}" org/ new/
    rm -rf $tmpdir
    exit 0
fi

# Find the changed files
git diff --no-index --name-only --exit-code org/ new/ | \
    sed -e 's#new/##' > $tmpdir/changed
changes=
if [[ -s $tmpdir/changed ]]
then
    changes=1
fi

# Check if changed files have changed outside the index
if grep -Ff $tmpdir/localmods $tmpdir/changed > $tmpdir/conflicts
then
    awk '{print $0 ": has changes in work tree"}' $tmpdir/conflicts \
        >> $tmpdir/messages
    if [[ ! $force && $action == update-* ]] ; then
        echo "Modified files found in work tree, skipping update. Use -f to override."
        echo "The following would have been done:"
        sort $tmpdir/messages
        rm -rf $tmpdir
        exit 2
    fi
fi

# Update the index/work tree if requested
if [[ $action == update-index ]] ; then
    grep -Ff $tmpdir/changed $tmpdir/filtered > $tmpdir/tohash
    cd $tmpdir/new
    IFS='
'
    for change in `cut -f2 $tmpdir/tohash | \
                   git --git-dir=$srcdir/.git hash-object -w --stdin-paths --no-filters | \
                   paste - $tmpdir/tohash`
    do
        # NOTE: the patterns below contain literal tabs
        sha1=${change%%	*}
        rest=${change#*	}
        mode=${rest:8:6}
        path=${rest#*	}
        path=${path%%	*}
        # Contains a literal tab
        echo "$mode $sha1	$path" >> $tmpdir/toindex
    done
    unset IFS
    git --git-dir=$srcdir/.git update-index --index-info < $tmpdir/toindex
elif [[ $action == update-workdir ]] ; then
    rsync --files-from=$tmpdir/changed $tmpdir/new/ $srcdir/
fi

# Get back to the original directory
popd >/dev/null

# Report what was done
sort $tmpdir/messages | tee $warning_file

rm -rf $tmpdir
exit $changes
