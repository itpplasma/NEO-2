#!/bin/bash
GITHASH=`git rev-parse HEAD`
GITSHORTHASH=`git rev-parse --short HEAD`
GITCHANGEDFILES=`git diff-index --name-only HEAD`

export GITSHORTHASH=$GITSHORTHASH

echo "character(len=*), parameter :: Neo2_Version = '${GITHASH}'" > ./version.f90
echo ""
echo "Git Hash is ${GITHASH}"

if [ -n "$GITCHANGEDFILES" ]; then
    echo 'character(len=*), parameter :: Neo2_Version_Additional = "WARNING, &' >> ./version.f90
    echo "&THERE ARE UNCOMMITTED CHANGES. Run may not be reproduceable: &" >> ./version.f90
    echo "THERE ARE UNCOMMITTED CHANGES!"
    echo ""
    while read -r line; do
        echo "&${line} &" >> ./version.f90
    done <<< "$GITCHANGEDFILES"
    echo '&"' >> ./version.f90

else

    echo 'character(len=*), parameter :: Neo2_Version_Additional = ""' >> ./version.f90

fi
