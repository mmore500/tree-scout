#!/bin/sh -ex

black tree_scout tests scripts
ruff check tree_scout tests scripts --fix
