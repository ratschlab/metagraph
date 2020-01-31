#!/bin/bash

# exit on error
set -e

# ask the client to shutdown
curl localhost:8001/quit

echo "Shutdown script finished."
