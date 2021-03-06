# for running "make"
# - uses CONDA_PREFIX as CASM_PREFIX default value

### initialization - shouldn't need to touch
set -e
export CASM_BUILD_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
. $CASM_BUILD_DIR/build_scripts/install-functions.sh
detect_os
check_var "CASM_PREFIX" "CASM install location" "$CONDA_PREFIX"

### end initialization ###

### variables - Control how CASM is built  ###

check_var "CASM_CXXFLAGS" "Compiler flags" ""
check_var "CASM_NCPU" "Compiler -j option" 2

# set OS-dependent variable defaults
. $CASM_BUILD_DIR/build_scripts/variables-$CASM_OS_NAME.sh

### end variables ###


bash $CASM_BUILD_DIR/build_scripts/make-cpp.sh
