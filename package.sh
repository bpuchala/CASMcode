#!/bin/bash

./version.sh

version=$(cat build-aux/casm_version.txt)
reversion="0.3a0-${version#*-}"

echo $reversion > build-aux/casm_version.txt

autoreconf -ifv

cd builds
../configure
make dist
mv casm-*-*.tar.gz $HOME/www/public_html/casm
