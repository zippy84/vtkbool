#!/bin/bash

[ ! -d coverage_report ] && mkdir coverage_report
docker build -t vtkbool_coverage -f Dockerfile.archlinux .
id=$(docker run -t -d vtkbool_coverage)
docker exec $id sh -c 'tar -cf - /master/build/coverage.*' | tar --strip-components=2 -xf - -C coverage_report
docker stop -t 0 $id
