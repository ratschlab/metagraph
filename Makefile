IMG_NAME := metagraph # final image name
IMG_NAME_DEV := metagraph_dev_env # development environment image name

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CODE_BASE_HOST := $(abspath $(MKFILE_PATH)/..)
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

alphabet := $(or $(alphabet), DNA)

additional_cmake_args := $(or $(cmake_args), -DBUILD_KMC=OFF)

ifeq ($(env), docker)
	CODE_BASE := /opt/metagraph
	BUILD_DIR := /opt/metagraph_build
	BUILD_DIR_STATIC := /opt/metagraph_build_static
else
	CODE_BASE := $(CODE_BASE_HOST)
	BUILD_DIR := $(BUILD_DIR_HOST)
	BUILD_DIR_STATIC := $(BUILD_DIR_STATIC_HOST)
endif

DATA_DIR := /data

# TODO: keep mechanism?
CONFIG_PATH := metagraph_local.mk

ifneq ("$(wildcard $(CONFIG_PATH))","")
        include $(CONFIG_PATH)
endif

DOCKER_OPTS := -it -u `id -u ${USER}`:$(DOCKER_GRP) -v $(BUILD_DIR_HOST_DOCKER):${BUILD_DIR} -v $(BUILD_DIR_STATIC_HOST_DOCKER):${BUILD_DIR_STATIC} -v  $(CODE_BASE_HOST):$(CODE_BASE) -v $(DATA_DIR_HOST):/data
DOCKER_BASE_CMD := docker run --rm $(DOCKER_OPTS) $(IMG_NAME_DEV)

ifeq ($(env), docker)
 	$(shell mkdir -p $(BUILD_DIR_HOST_DOCKER) $(BUILD_DIR_STATIC_HOST_DOCKER) )
	EXEC_CMD := $(DOCKER_BASE_CMD) bash -c
else
	EXEC_CMD := bash -c
endif

all: update build test

ifeq ($(env), docker)
update: sync build-docker-dev-env
else
update: sync
endif

sync:
	cd $(CODE_BASE_HOST) && git submodule sync && git submodule update --init --recursive

## BUILDING BINARIES AND IMAGES

build: build-sdsl-lite build-metagraph

build-sdsl-lite:
	$(EXEC_CMD) 'cd $(CODE_BASE)/metagraph/external-libraries/sdsl-lite && ./install.sh $${PWD}'

build-metagraph:
	[ -d $(BUILD_DIR_HOST_DOCKER) ] || mkdir -p $(BUILD_DIR_HOST_DOCKER)
	$(EXEC_CMD) 'mkdir -p $(BUILD_DIR) && cd $(BUILD_DIR) && cmake -DCMAKE_DBG_ALPHABET=$(alphabet) $(additional_cmake_args) $(CODE_BASE)/metagraph && make metagraph -j $$(($$(getconf _NPROCESSORS_ONLN) - 1))'

build-metagraph-static:
	[ -d $(BUILD_DIR_STATIC_HOST_DOCKER) ] || mkdir -p $(BUILD_DIR_STATIC_HOST_DOCKER)
	$(EXEC_CMD) 'mkdir -p $(BUILD_DIR_STATIC) && cd $(BUILD_DIR_STATIC) && cmake -DCMAKE_DBG_ALPHABET=$(alphabet) -DBUILD_STATIC=ON $(additional_cmake_args) $(CODE_BASE)/metagraph && make metagraph -j $$(($$(getconf _NPROCESSORS_ONLN) - 1))'

build-docker:
	docker build -t metagraph $(CODE_BASE_HOST)

# explicitly generating intermediate stages
build-docker-dev-env:
	docker build --target metagraph_dev_env -t metagraph_dev_env $(CODE_BASE_HOST)

build-docker-bin:
	docker build --target metagraph_bin -t metagraph_bin $(CODE_BASE_HOST)


# TESTING

test: integration-tests

integration-tests:
	$(EXEC_CMD) 'cd $(BUILD_DIR) && ./integration_tests'

integration-tests-api:
	$(EXEC_CMD) 'cd $(BUILD_DIR) && ./integration_tests --test_filter="test_api*"'




