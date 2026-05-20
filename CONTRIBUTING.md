# Contributing to MetaGraph

## Building from source

The full installation guide (custom alphabets, debug builds, dependencies) lives in the [online docs](https://metagraph.ethz.ch/static/docs/installation.html#install-from-source).

### Build a docker container

```bash
docker build .
```

### Makefile shortcuts

The top-level `Makefile` wraps the common build / test invocations. Useful arguments:

- `env`: `""` (host) or `docker`
- `alphabet`: e.g. `DNA`, `DNA5`, `Protein` (default `DNA`)
- `additional_cmake_args`: extra flags forwarded to CMake

Example:

```bash
# compile in a docker container for the DNA alphabet
make build-metagraph env=docker alphabet=DNA
```

## Releases

1. Bump the version in `package.json`.
2. Tag the commit with the new version.
3. Create a GitHub release.
