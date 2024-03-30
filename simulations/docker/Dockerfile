FROM python:3.11-slim-bullseye

RUN \
	set -ex; \
	apt-get update; \
	apt-get install -y libgsl-dev

RUN \
	set -e; \
	# create virtual environment
	python -m venv ve; \
	# update pip, setuptools
	ve/bin/pip install -U pip setuptools; \
	# install requirements
	ve/bin/pip install fwdpy11

ADD run_replicates.py /opt
ADD demo.yaml /opt

WORKDIR /opt

ENTRYPOINT ["/ve/bin/python", "/opt/run_replicates.py"]
