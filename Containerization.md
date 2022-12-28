layout: page
title: "Containerization"
permalink: /Containerization

## <a name='TOC'>Containerization with Docker</a>

1. [Building containers](#Build)
2. [Running containers](#Run)
3. [Remove all stopped containers](#Prune)
4. [Reattach to a stopped container](#Reattach)


**<a name="Build">Building containers</a>**

First, use a text editor to make a Dockerfile. Here is an example:
```bash
FROM debian:latest

# Update packages in base image, avoid caching issues by combining statements, install build software and deps
RUN	apt-get update && apt-get install -y build-essential git pkg-config libssl-dev bzip2 wget zlib1g-dev libswscale-dev gettext nettle-dev libgmp-dev libssh2-1-dev libgnutls28-dev libc-ares-dev libxml2-dev libsqlite3-dev autoconf libtool libcppunit-dev automake autotools-dev autopoint && \
    #Install aria2 from git, cleaning up and removing all build footprint	
    git clone https://github.com/tatsuhiro-t/aria2.git /opt/aria2 && \
    cd /opt/aria2 && \
    autoreconf -i && ./configure && \
    make && make install 

CMD ["/usr/local/bin/aria2c","--conf-path=/config/aria2.conf"]
```

To build an image named `image_name`, use the provided Dockerfile. Note that as part of the `build` command, *all files in the same directory as the Dockerfile are included in the build*. So, the Dockerfile should probably be placed in an empty directory.

```bash
docker build -t image_name path/to/Dockerfile
```

**<a name="Run">Running containers</a>**

Run a docker image, interactively, named 'downloader', mounting the current working directory to `mnt` (note that the `-v` command expects the *absolute path* to the directories):
```bash
docker run -it -v "$(pwd):/mnt" downloader
```

**<a name="Prune">Remove all stopped containers</a>**

The Docker [Prune](https://docs.docker.com/engine/reference/commandline/container_prune/) command is a useful way to remove a lot of stopped containers at the same time
```bash
docker container prune
```

**<a name="Reattach">Reattach to a stopped container</a>**

To re-attach (-a) to a stopped docker container, interactively (-i)
```bash
docker start -a -i {container_id}
```
