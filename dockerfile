# syntax=docker/dockerfile:1

FROM debian:bullseye-slim AS build-container

WORKDIR /app

RUN  apt-get update \
  && apt install -y python3-pip git \
  && rm -rf /var/lib/apt/lists/* 

COPY pyproject.toml README.md pdm.lock LICENSE simulate.py /app/
COPY ./src /app/src
COPY ./tests /app/tests
# Install pdm
RUN pip install -U pdm

# disable update check
ENV PDM_CHECK_UPDATE=false
# Install all the dependencies
RUN pdm install --check --prod --no-editable

RUN pdm simulate
RUN pdm test

FROM debian:bullseye-slim AS app-container

RUN  apt-get update \
  && apt install -y bedtools python3-pip \
  && rm -rf /var/lib/apt/lists/* 

# retrieve packages from build stage
COPY --from=build-container /app/.venv/ /app/.venv
ENV PATH="/app/.venv/bin:$PATH"