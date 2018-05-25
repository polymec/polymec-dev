# generate_headers.sh
# usage: sh generate_headers.sh source_dir dest_dir lib1 lib2 ... libN
# This script generates polymec.h and friends at the top level of the 
# destination directory, including header files found within the source 
# directory tree.

# Silent pushd.
pushd () 
{
  command pushd "$@" > /dev/null
}

# Silent popd.
popd () 
{
  command popd "$@" > /dev/null
}

# This function generates a top-level header for a library.
gen_lib_header()
{
  lib=$1
  source_dir=$2
  dest_dir=$3

  # Change to our source directory.
  pushd $source_dir

  # Write a temporary file with the new contents.
  tmp=$dest_dir/polymec_$lib.h.tmp
  echo "// polymec_$lib.h -- automatically generated." > $tmp
  cat >> $tmp <<- END1
// This file is part of the polymec HPC library. See the license
// in the actual source files for details of distribution.

END1
  echo "#ifndef POLYMEC_${lib}_LIBRARY_H" >> $tmp
  echo "#define POLYMEC_${lib}_LIBRARY_H" >> $tmp
  cat >> $tmp <<- END2
#ifdef __cplusplus
extern "C" {
#endif

END2
  if [ "$lib" = "core" ]; then
    echo "#include \"core/polymec.h\"" >> $tmp
  fi
  find $lib -regex ".*/[a-zA-Z0-9_]*.h" -exec echo "#include \"{}\"" >> $tmp \;
  cat >> $tmp <<- END3
#ifdef __cplusplus
}
#endif

#endif
END3

  header=$dest_dir/polymec_$lib.h
  if [ -f $header ]; then
    # Diff the temporary file against what's there.
    tmp_diff_wc=`diff $header $tmp | wc -l`

    # If the files differ, move the temporary file into place. Otherwise
    # discard it.
    if [ $tmp_diff_wc -gt 0 ]; then
      mv $tmp $header
    else
      rm $tmp
    fi
  else
    mv $tmp $header
  fi

  # Change back to the original directory.
  popd 
}

# This function generates the highest-level header for polymec.
gen_top_header()
{
  source_dir=$1
  dest_dir=$2
  shift; shift
  libs=$@;

  # Change to our source directory.
  pushd $source_dir

  # Write a temporary file with the new contents.
  tmp=$dest_dir/polymec.h.tmp
  echo "// polymec.h -- automatically generated." > $tmp
  cat >> $tmp <<- END1
// This file is part of the polymec HPC library. See the license
// in the actual source files for details of distribution.
#ifndef POLYMEC_LIBRARY_H
#define POLYMEC_LIBRARY_H

END1
  for lib in $libs
  do
    echo "#include \"polymec_$lib.h\"" >> $tmp
  done
  cat >> $tmp <<- END2

#endif
END2

  header=$dest_dir/polymec.h
  if [ -f $header ]; then
    # Diff the temporary file against what's there.
    tmp_diff_wc=`diff $header $tmp | wc -l`

    # If the files differ, move the temporary file into place. Otherwise
    # discard it.
    if [ $tmp_diff_wc -gt 0 ]; then
      mv $tmp $header
    else
      rm $tmp
    fi
  else
    echo "$tmp -> $header"
    mv $tmp $header
  fi

  # Change back to our original directory.
  popd
}

source_dir=$1
dest_dir=$2
shift; shift
libs=$@

for lib in $libs
do
  gen_lib_header $lib $source_dir $dest_dir
done

gen_top_header $source_dir $dest_dir $libs
