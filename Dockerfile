FROM debian:11 AS build_image

RUN groupadd -g 901 superSTR_grp && useradd -u 901 -g superSTR_grp -ms /bin/sh superSTR_user

RUN apt-get update && apt-get install -y \
 autoconf \
 build-essential \
 cmake \
 git \
 libbz2-dev \
 libcurl4-openssl-dev \
 liblzma-dev \
 libncurses-dev \
 zlib1g-dev \
 wget

WORKDIR /tmp/htslib/
RUN git clone https://github.com/samtools/htslib
WORKDIR /tmp/htslib/htslib/
RUN git checkout tags/1.12
RUN git submodule update --init --recursive
RUN autoheader && autoconf && ./configure --prefix=/usr/local/superstr/htslib
RUN make && make install
ENV HTSLIB_ROOT /usr/local/superstr/htslib

WORKDIR /usr/local/sratoolkit/
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz 
RUN tar --strip-components=1 -zxvf sratoolkit.2.9.6-ubuntu64.tar.gz
RUN rm sratoolkit.2.9.6-ubuntu64.tar.gz

USER superSTR_user
ENV HOME /home/superSTR_user/
WORKDIR /home/superSTR_user/
RUN git clone https://github.com/bahlolab/superSTR
WORKDIR /home/superSTR_user/superSTR/C/
ENV PATH="/usr/local/sratoolkit/bin:/home/superSTR_user/superSTR/C:${PATH}"
RUN cmake .
RUN make
RUN mkdir /home/superSTR_user/.ncbi/
RUN mkdir /home/superSTR_user/sra-cache/
ENV HTSLIB_ROOT /usr/local/superstr/htslib
ENV PATH="/usr/local/sratoolkit/bin:${PATH}"
RUN vdb-config --cfg-dir /home/superSTR_user/.ncbi/
RUN echo '/config/default = "false"' > /home/superSTR_user/.ncbi/user-settings.mkfg
RUN echo '/repository/remote/disabled = "false"' >> /home/superSTR_user/.ncbi/user-settings.mkfg
RUN echo '/repository/user/cache-disabled = "false"' >> /home/superSTR_user/.ncbi/user-settings.mkfg
RUN echo '/repository/user/default-path = "/usr/local/sratoolkit/"' >> /home/superSTR_user/.ncbi/user-settings.mkfg
RUN echo '/repository/user/main/public/root = "/sra-cache/"' >> /home/superSTR_user/.ncbi/user-settings.mkfg
RUN vdb-config --root -s /repository/site/main/public/superSTR_user/sra-cache=/home/superSTR_user/sra-cache

USER root
RUN apt-get clean
