# update_version_h.sh
# usage: bash update_version_h.sh source_dir dest_dir/polymec_version.h

# Search source_dir/CMakeLists.txt for version information.
major_version=`grep '(POLYMEC_MAJOR_VERSION' $1/CMakeLists.txt | sed -e 's/set (//' -e 's/)$//' -e 's/POLYMEC_MAJOR_VERSION //'`
minor_version=`grep '(POLYMEC_MINOR_VERSION' $1/CMakeLists.txt | sed -e 's/set (//' -e 's/)$//' -e 's/POLYMEC_MINOR_VERSION //'`
patch_version=`grep '(POLYMEC_PATCH_VERSION' $1/CMakeLists.txt | sed -e 's/set (//' -e 's/)$//' -e 's/POLYMEC_PATCH_VERSION //'`

# Do we have git?
which_git_wc=`which git | wc -l`

if [ $which_git_wc -gt 0 ]; then
  # Ask git for the revision and for diff information.
  git_revision=`git log -1 --format=format:%h`
  git_diffs=`git diff | wc -l` 
else
  git_revision=""
  git_diffs=0
fi

# Write a temporary file with the new contents.
cat > $2.tmp <<- END
// This file is automagically generated by update_version_h.py.
#ifndef POLYMEC_VERSION_H
#define POLYMEC_VERSION_H

END

echo "#define POLYMEC_MAJOR_VERSION $major_version" >> $2.tmp
echo "#define POLYMEC_MINOR_VERSION $minor_version" >> $2.tmp
echo "#define POLYMEC_PATCH_VERSION $patch_version" >> $2.tmp
echo "static const char* POLYMEC_VERSION = \"$major_version.$minor_version.$patch_version\";" >> $2.tmp
echo "static const char* POLYMEC_REVISION = \"$git_revision\";" >> $2.tmp
if [ $git_diffs -gt 0 ]; then
  echo "static const bool POLYMEC_HAS_DIFFS = true;" >> $2.tmp
else
  echo "static const bool POLYMEC_HAS_DIFFS = false;" >> $2.tmp
fi
echo "#endif" >> $2.tmp

if [ -f $2 ]; then
  # Diff the temporary file against what's there.
  tmp_diff_wc=`diff $2 $2.tmp | wc -l`

  # If the files differ, move the temporary file into place. Otherwise
  # discard it.
  if [ $tmp_diff_wc -gt 0 ]; then
    mv $2.tmp $2
  else
    rm $2.tmp
  fi
else
  mv $2.tmp $2
fi
