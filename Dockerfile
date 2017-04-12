FROM ubuntu:16.04

MAINTAINER Ralf Weber, r.j.weber@bham.ac.uk

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF && \
    echo "deb http://download.mono-project.com/repo/debian wheezy/snapshots/4.8.1.0 main" | sudo tee /etc/apt/sources.list.d/mono-xamarin.list && \
    apt-get -y update && apt-get -y install --no-install-recommends libglib2.0-dev mono-complete python-dev python-pip && \
    pip install --upgrade pip && pip install -U setuptools && \
    pip install dimspy==0.1.0 && \
    pip uninstall -y pip && \
    apt-get purge -y python-pip && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENTRYPOINT ["dimspy"]
