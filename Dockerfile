FROM archlinux:latest
RUN pacman -Syyu --noconfirm
RUN pacman -S cmake vtk make gcc libjpeg-turbo libpng libtiff python-pytest openmpi fmt python gcovr --noconfirm
WORKDIR /refactoring
COPY . .
RUN ls -alh
WORKDIR /refactoring/build
RUN cmake -DCMAKE_BUILD_TYPE=PROFILE -DVTKBOOL_COVERAGE=ON ..
RUN make VERBOSE=1
RUN ctest --output-on-failure
RUN gcovr -r .. . --cobertura-pretty -o coverage.xml

FROM ubuntu:latest
RUN apt-get update
RUN apt-get -y upgrade
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata
RUN apt-get -y install libvtk9-dev python3-vtk9 cmake python3-pytest
WORKDIR /refactoring
COPY . .
RUN ls -alh
WORKDIR /refactoring/build
RUN cmake ..
RUN make
RUN ctest --output-on-failure
