#!/bin/bash

if [ -d build2/ ]
then
    cd build2

    modes=(union intersection difference difference2)

    for i in $(seq 0 6); do
        for j in 0 1 2 3; do
            ./tests $i ${modes[j]} > deb_$i.txt

            if [ ! -s deb_$i.txt ]
            then
                rm deb_$i.txt
            fi
        done
    done

fi
