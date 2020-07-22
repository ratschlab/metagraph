IMG_NAME := metagraph # final image name
IMG_NAME_DEV := metagraph_dev_env # development environment image name

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CODE_BASE_HOST := $(abspath $(MKFILE_PATH)/../../..)
BUILD_DIR_HOST := $(CODE_BASE_HOST)/metagraph/build
BUILD_DIR_STATIC_HOST := $(CODE_BASE_HOST)/metagraph/build_static

# build dir on host when used with the build environment provided by metagraph_dev_env
BUILD_DIR_HOST_DOCKER := $(CODE_BASE_HOST)/metagraph/build_docker
BUILD_DIR_STATIC_HOST_DOCKER := $(CODE_BASE_HOST)/metagraph/build_docker_static

OS := $(shell uname -s)

ifeq ($(OS), Darwin)
DOCKER_GRP = 0
else
DOCKER_GRP = `stat -c '%g' /var/run/docker.sock`
endif

ifeq ($(docker), false)
	CODE_BASE := $(CODE_BASE_HOST)
	BUILD_DIR := $(BUILD_DIR_HOST)
	BUILD_DIR_STATIC := $(BUILD_DIR_STATIC_HOST)
else
	CODE_BASE := /opt/metagraph
	BUILD_DIR := /opt/metagraph_build
  BUILD_DIR_STATIC := /opt/metagraph_build_static
endif

DATA_DIR := /data

CONFIG_PATH := metagraph_local.mk

ifneq ("$(wildcard $(CONFIG_PATH))","")
        include $(CONFIG_PATH)
endif

DOCKER_OPTS := -it -u `id -u ${USER}`:$(DOCKER_GRP) -v $(BUILD_DIR_HOST_DOCKER):${BUILD_DIR} -v $(BUILD_DIR_STATIC_HOST_DOCKER):${BUILD_DIR_STATIC} -v  $(CODE_BASE_HOST):$(CODE_BASE) -v $(DATA_DIR_HOST):/data
DOCKER_BASE_CMD := docker run --rm $(DOCKER_OPTS) $(IMG_NAME_DEV)

ifeq ($(docker), false)
  EXEC_CMD := bash -c
else
  EXEC_CMD := $(DOCKER_BASE_CMD) bash -c
endif

sync:
	cd $(CODE_BASE_HOST) && git submodule sync && git submodule update --init --recursive

ifneq ($(docker), false)
build-docker-dev-env:
	cd $(CODE_BASE_HOST)/metagraph && docker build . -f docker/Dockerfile_dev_env -t metagraph_dev_env
else
build-docker-dev-env:
	echo "docker deactivated"
endif


build-docker-bin:
	cd $(CODE_BASE_HOST)/metagraph && docker build . -f docker/Dockerfile_bin -t metagraph_bin

BIN_TMP := $(CODE_BASE_HOST)/metagraph/docker/bin_tmp

build-docker:
	[ -d $(BIN_TMP) ] || mkdir -p $(BIN_TMP)
	docker create --name metagraph_bin_tmp metagraph_bin
	docker cp metagraph_bin_tmp:/opt/metagraph/metagraph/build/metagraph_DNA $(BIN_TMP)/metagraph_docker
	docker cp metagraph_bin_tmp:/usr/local/lib/libhts.so.1.10.2 $(BIN_TMP)/libhts_docker.so
	docker rm metagraph_bin_tmp
	cd $(CODE_BASE_HOST)/metagraph && docker build . -f docker/Dockerfile -t metagraph

build-docker-all: build-docker-dev-env build-docker-bin build-docker

build-sdsl-lite:
	$(EXEC_CMD) 'cd $(CODE_BASE)/metagraph/external-libraries/sdsl-lite && ./install.sh $${PWD}'

build-libmaus:
	$(EXEC_CMD) 'cd $(CODE_BASE)/metagraph/external-libraries/libmaus2; rm -f CMakeCache.txt; cmake -DCMAKE_INSTALL_PREFIX:PATH=$${PWD} . && make -j $$(($$(getconf _NPROCESSORS_ONLN) - 1)) && make install'

build-metagraph:
	[ -d $(BUILD_DIR_HOST_DOCKER) ] || mkdir -p $(BUILD_DIR_HOST_DOCKER)
	$(EXEC_CMD) 'mkdir -p $(BUILD_DIR) && cd $(BUILD_DIR) && cmake $(CODE_BASE)/metagraph && make -j $$(($$(getconf _NPROCESSORS_ONLN) - 1))'

build-metagraph-static:
	[ -d $(BUILD_DIR_STATIC_HOST_DOCKER) ] || mkdir -p $(BUILD_DIR_STATIC_HOST_DOCKER)
	$(EXEC_CMD) 'mkdir -p $(BUILD_DIR_STATIC) && cd $(BUILD_DIR_STATIC) && cmake -DBUILD_STATIC=ON $(CODE_BASE)/metagraph && make metagraph -j $$(($$(getconf _NPROCESSORS_ONLN) - 1)) metagraph'

api-integration-tests:
	$(EXEC_CMD) 'cd $(BUILD_DIR) && ./integration_tests --test_filter="test_api*"'

integration-tests:
	$(EXEC_CMD) 'cd $(BUILD_DIR) && ./integration_tests'

build-server: sync build-docker-dev-env build-sdsl-lite build-libmaus build-metagraph api-integration-tests

start-server:
	docker run -d --restart=always -it $(DOCKER_OPTS) -p 5555:5555 metagraph server_query -i $(DATA_DIR)/bigsi/graph_subset_3000.dbg -a $(DATA_DIR)/bigsi/annotation_subset_3000.relaxed.brwt.annodbg -v -p 8

run-interactive:
	docker run -it $(DOCKER_OPTS) --entrypoint bash $(IMG_NAME_DEV)
