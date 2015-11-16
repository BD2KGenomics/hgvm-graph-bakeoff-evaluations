#!/usr/bin/env bash

# toilAzureCleanup.sh: delete all the tables and containers that Toil makes,
# assuming you named all your Toil Azure job stores with the prefix "tree". DO
# NOT RUN THIS IF YOU HAVE ANY IMPORTANT DATA in anything on Azure named with
# the prefix "atree" or in the "toilRegistry" table.

# Requires AZURE_CONNECTION_STRING to be set.

set -e

# Safety interlock
read -p "Deleting toilRegistry and atree* from Azure. Type 'delete' to continue: " DELETE

if [ "$DELETE" == "delete" ]
then
    # Actually do it
    echo "Deleting..."
    azure storage table list --json | jq -r '.[].name' | grep "^atree" | xargs -n1 -r azure storage table delete -q
    azure storage table list --json | jq -r '.[].name' | grep '^toilRegistry$' | xargs -n1 -r azure storage table delete -q
    azure storage container list --json | jq -r '.[].name' | grep "^atree" | xargs -n1 -r azure storage container delete -q
else
    echo "Not deleting anything."
fi

