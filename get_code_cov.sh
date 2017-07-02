#!/bin/bash
for filename in `find 'libvipster' | egrep '\.cpp'`; 
do 
    gcov-5 -n -o . $filename > /dev/null; 
done
