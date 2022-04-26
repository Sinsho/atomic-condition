FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils

RUN apt-get install -y --no-install-recommends \
  ca-certificates \
  git \
  vim \
  build-essential \
  clang-3.9 \
  llvm-3.9 \
  python3 \
  python3-setuptools \
  python3-pip \
  autoconf \
  automake \
  libtool \
  libgoogle-glog-dev \
  libgflags-dev \
  libeigen3-dev

#libeigen is needed for libcmaes

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

# Ready for building solver.
WORKDIR /atom

