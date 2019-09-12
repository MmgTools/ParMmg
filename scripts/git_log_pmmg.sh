#!/bin/bash

# get arguments
CMAKE_SOURCE_DIR=$1
PMMG_SOURCE_DIR=$2
HEADER_FILE="git_log_pmmg.h"

# test if $HEADER_FILE is present
if [ ! -f "$PMMG_SOURCE_DIR/$HEADER_FILE" ]; then
   \touch $PMMG_SOURCE_DIR/$HEADER_FILE
   echo "#define PMMG_GIT_BRANCH \"No git branch found\"" >> "$PMMG_SOURCE_DIR/$HEADER_FILE"
   echo "#define PMMG_GIT_COMMIT \"No git commit found\"" >> "$PMMG_SOURCE_DIR/$HEADER_FILE"
   echo "#define PMMG_GIT_DATE   \"No git commit found\"" >> "$PMMG_SOURCE_DIR/$HEADER_FILE"
fi

# parameter
git_is_present=1

# test if .git is present
if [ ! -d $CMAKE_SOURCE_DIR/.git ]; then
   git_is_present=0
fi

git_commit=`git rev-parse HEAD 2> /dev/null`
if [ "$git_commit" == "" ]; then
   git_is_present=0
fi

if [ $git_is_present -eq 1 ]; then
   git_branch=`git rev-parse --abbrev-ref HEAD 2> /dev/null`
   echo "   > Found a git branch: $git_branch"
   echo "   > Found a git commit: $git_commit"
   git_date=`git show -s --format="%ci" $git_commit`
   echo "   > Found a git date: $git_date"

   # erase any previous version file
   \rm -f $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp"
   \touch $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp"
   echo "#define PMMG_GIT_BRANCH \"$git_branch\"" >> $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp"
   echo "#define PMMG_GIT_COMMIT \"$git_commit\"" >> $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp"
   echo "#define PMMG_GIT_DATE   \"$git_date\""   >> $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp"

   # diff
   diff=`diff $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp" "$PMMG_SOURCE_DIR/$HEADER_FILE"`
   if [ "$diff" != "" ]; then
      \cp $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp" "$PMMG_SOURCE_DIR/$HEADER_FILE"
   fi

   # clean
   \rm -f $PMMG_SOURCE_DIR/$HEADER_FILE"_tmp"
fi
