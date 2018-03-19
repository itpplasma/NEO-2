#!/bin/bash
GITHASH=`git rev-parse HEAD`
GITSHORTHASH=`git rev-parse --short HEAD`
GITCHANGEDFILES=`git diff-index --stat HEAD`
GITHASHDATE=$(git log -n 1 --date=format:'%A %Y-%m-%d %H:%M:%S' --pretty=format:"%cd" $GITHASH)

export GITSHORTHASH=$GITSHORTHASH

echo "character(len=*), parameter :: Neo2_Version = '${GITHASH}'" > ./version.f90
echo "character(len=*), parameter :: Neo2_Version_Date = '${GITHASHDATE}'" >> ./version.f90
echo ""
echo "Git Hash is ${GITHASH}"
echo "was produced on ${GITHASHDATE}"

if [ -n "$GITCHANGEDFILES" ]; then
    echo 'character(len=*), parameter :: Neo2_Version_Additional = "WARNING, &' >> ./version.f90
    echo "&THERE ARE UNCOMMITTED CHANGES. Run may not be reproduceable. &" >> ./version.f90
    echo "&Changed Files to Git Revision:&" >> ./version.f90
    echo "THERE ARE UNCOMMITTED CHANGES!"
    echo ""
    while read -r line; do
        echo '&"//char(13)//char(10)//"&' >> ./version.f90
        echo "& ${line} &" >> ./version.f90
    done <<< "$GITCHANGEDFILES"
    echo '&"' >> ./version.f90

else

    echo 'character(len=*), parameter :: Neo2_Version_Additional = ""' >> ./version.f90

fi
