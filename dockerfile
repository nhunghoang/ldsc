# syntax=docker/dockerfile:1

FROM ubuntu:22.04

WORKDIR /app

RUN  apt-get update \
  && apt-get install -y wget python3.10 bedtools python3-pip \
  && rm -rf /var/lib/apt/lists/* \
  && pip install --upgrade pip

# creating a folder for ldsc
RUN mkdir ldsc

COPY make_annot.py /app/ldsc/
COPY munge_sumstats.py /app/ldsc/
COPY ldsc.py /app/ldsc/
COPY ./ldscore/ /app/ldsc/ldscore/

# install all of the python dependencies
COPY requirements.txt /app/

RUN pip install -r requirements.txt



