#!/bin/bash
version=$(date +%F).$1
echo $version

sed -i "s#Version.*#Version ${version})#" ProjectVersion.cmake.in
