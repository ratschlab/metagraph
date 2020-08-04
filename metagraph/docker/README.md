# Dockerized Metagraph

To easily and flexibly deploy the metagraph server the `metagraph` docker image can be used.

This is still somewhat work in progress, but should be already usable now.

# Makefile usage

The Makefile in the `docker` directory serves the purpose to simplify the building of metagraph binaries and especially docker images.

By default, commands are run in a docker container, but by passing `docker=false` commands are run on the host itself.

# Docker Images

There are several docker images involved, the main reason for this, that the current main deployment environment doesn't support multi-stage builds.

The images are as follows:
1. `metagraph_dev_env`: contains all dependencies to build metagraph. Can also be used during development by mounting the code base and build dir on the host (this is done in the `build-metagraph` make target)
2. `metagraph_bin`: based on `docker_dev_env` but contains the `metagraph` binary. It is more of an intermediary image and not really used per se.
3. `metagraph`: "the final output", the image used in production. It contains a basic runtime environment for metagraph (no build tools) along with the metagraph binary. The binary is copied out of the `metagraph_bin` image. It also contains the python API code. This image could also be published on dockerhub later on.

## Example Usage

### Build metagraph image

* `make build-docker-all`: takes care of rebuilding all 3 docker images
* `make build-docker`: just recreates the `metagraph` image

### Run metagraph in docker

To run the metagraph binary in a container, do
```
docker run --rm metagraph
```

In order to do something useful, you'd want to mount some data into the container, e.g.:

```
docker run --rm -v /some_dir_on_host:/mnt metagraph stats /mnt/mygraph.dbg
```

(see for more details see e.g. https://docs.docker.com/storage/bind-mounts/)

To start a bash in a metagraph container (e.g. to explore what is in the image), do
```
docker run -it --rm --entrypoint bash metagraph
```

### Run local test server in docker container

For local testing only, production deployment currenlty uses a docker-compose setup (not included here)

```
make build-server
make start-server
```

The latter command may needs some tweaking in the `Makefile` on your system. The data directory on the host can be set using an additional file `metagraph_local.mk` in the same directory as the `Makefile`. The directory in `DATA_DIR_HOST` gets mapped to `/data` in the container. The variables `SERVER_GRAPH_FILE` and `SERVER_ANNOT_FILE` are paths in the container. In the example below, the graph file would be in `/local/metagraph/my_graph.dbg` on the host filesystem.

```
DATA_DIR_HOST := /local/metagraph
SERVER_GRAPH_FILE := /data/my_graph.dbg
SERVER_ANNOT_FILE := /data/my_graph.relaxed.brwt.annodbg

```
