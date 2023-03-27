#!/bin/sh -ex

mypy tree_scout tests
black tree_scout tests --check
ruff tree_scout tests scripts
