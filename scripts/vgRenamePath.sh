#!/bin/bash
# rename a path.  this is not a clever script, and relies on the
# given path name being the only name in the graph.  new vg
# printed to stdout

set -e

VG=$1
VG_PATH=$2
VG_NEW_PATH=$3

vg view $1 -j | jq . | sed -e "s/\"name\"\:\ \"${2}\",/\"name\"\:\ \"${3}\",/" | vg view -Jv -

