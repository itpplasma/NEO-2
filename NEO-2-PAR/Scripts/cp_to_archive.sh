#!/bin/bash

arch=`uname -m`
branch=`git rev-parse --abbrev-ref HEAD`
shorthash=`git rev-parse --short HEAD`
gituncommitted=`git diff --exit-code --quiet; echo $?`
gituncommittedstr=""
if [ $gituncommitted -eq "1" ]; then
	gituncommittedstr="-INCL_UNCOMMITTED"
fi

cp ${1} ${2}/neo_2.x-${3}-${4}-${branch}-${shorthash}${gituncommittedstr}
echo "Copy file to ${2}/neo_2.x-${3}-${4}-${branch}-${shorthash}${gituncommittedstr}"
