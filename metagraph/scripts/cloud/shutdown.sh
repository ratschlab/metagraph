#!/bin/bash

# exit on error
set -e

# ask the client to shutdown
echo "Shutdown script invoked, calling /quit"
curl localhost:8001/quit

if [ -f "/tmp/shutdown.sh" ]; then
  echo "Calling /tmp/shutdown.sh"
  bash /tmp/shutdown.sh
fi

echo "Shutdown script finished."
