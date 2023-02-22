FROM debian:stable-slim
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda init bash

RUN conda install -c conda-forge mamba

COPY . /PitViper

RUN cd /PitViper/ \
    && bash install_PitViper_env.sh

EXPOSE 5000
EXPOSE 8888

CMD cd /PitViper/PitViper && bash run.sh
