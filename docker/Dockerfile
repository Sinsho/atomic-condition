FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install -y --no-install-recommends \
  ca-certificates \
  git \
  build-essential \
  clang-3.9 \
  llvm-3.9 \
  python3 \
  python3-setuptools \
  python3-pip

#install dependencies for libcmaes
RUN apt-get update
RUN apt-get install -y \
  vim \
  autoconf \
  automake \
  libtool \
  libgoogle-glog-dev \
  libgflags-dev \
  libeigen3-dev

#install dependencies for nomad
RUN apt-get update
RUN apt-get install -y \
    unzip \
    wget \
    libssl-dev \
    gcc-8 \
    g++-8

RUN pip3 install wheel
RUN pip3 install mpmath

RUN git clone https://github.com/Sinsho/atomic-condition atom && \
  cd atom

# Build necessary objects.
WORKDIR /atom
RUN make prep
RUN tar xf gsl-2.5.tar.gz

# Build GSL
WORKDIR /atom/gsl-2.5
RUN ./configure --prefix=/atom/gslbuild CC=clang-3.9
RUN make -j4
RUN ln -s /atom/misc/compile-gsl-sf.sh /atom/gsl-2.5/specfunc
WORKDIR /atom/gsl-2.5/specfunc
RUN ./compile-gsl-sf.sh

# Build libcmaes
WORKDIR /atom
RUN git clone https://github.com/beniz/libcmaes.git
WORKDIR /atom/libcmaes
RUN ./autogen.sh
RUN echo "#define CMAES_EXPORT" > include/libcmaes/cmaes_export.h
RUN ./configure --prefix /usr
RUN make && make install
RUN ldconfig

# Build cmake since nomad needs a newer cmake version than the one available in the package repo from Ubuntu 18.04
WORKDIR /tmp
RUN wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1.tar.gz
RUN tar -xzvf cmake-3.23.1.tar.gz
WORKDIR /tmp/cmake-3.23.1
RUN ./bootstrap
RUN make -j$(nproc)
RUN make install

# Swap gcc version for nomad
RUN rm /usr/bin/gcc
RUN ln -s /usr/bin/gcc-8 /usr/bin/gcc
RUN rm /usr/bin/g++
RUN ln -s /usr/bin/g++-8 /usr/bin/g++

# Build nomads
WORKDIR /atom
ENV NOMAD_HOME=/atom/nomad
RUN git clone https://github.com/bbopt/nomad.git
WORKDIR /atom/nomad
RUN cmake -S . -B build/release
RUN cmake --build build/release
RUN cmake --install build/release
RUN cp -R /atom/nomad/build/release/lib/* /usr/lib/


# Ready for building solver.
WORKDIR /atom

