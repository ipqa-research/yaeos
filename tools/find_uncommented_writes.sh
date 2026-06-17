#!/bin/sh
grep -R "write([0-9]" ../test | grep -v '!.*write'
