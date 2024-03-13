#!/bin/bash
DESIRED_COVERAGE=90

DID_TEST=0

COVER=0

echoerr() { echo -e "$@" 1>&2; }

install_fpm() {
    apt install pipx
    pipx install fpm
}

green() {
    echo -e "\e[1;32m$@\e[m"
}

red() {
    echo -e "\e[1;31m$@\e[m"
}

run_test() {
    echo y | fpm clean
    DID_TEST=1
    echo Checking tests files names...
    NAMING_ERRORS=0
    for file in $(find test/*f90); do
        filename="$(basename $file)"
        prefix="$(echo $filename | cut -d '_' -f1)"

        if [ "$prefix" = "test" ]; then
            green "$file [OK]"
        else
            NAMING_ERRORS=$((NAMING_ERRORS + 1))
            red "$file [X]"
        fi
    done

    [ $NAMING_ERRORS -ge 1 ] && 
        echoerr "There are wrongly named files in the test directory"

    echo Running tests...
    fpm test --flag "--coverage"
}

run_coverage() {
    gcovr \
        --exclude "build" \
        --exclude "test/test_runner.f90" \
        --exclude "test/fixtures/taperobinson.f90" \
        --exclude "src/adiff/hyperdual.f90" \
        --exclude "example" \
        --exclude "src/legacy/*" \
        --exclude "app"\
        --exclude "tools" \
        --fail-under-line 90 \
        --jacoco coverage.xml
    
    gcovr \
        --exclude "build" \
        --exclude "test/test_runner.f90" \
        --exclude "test/fixtures/taperobinson.f90" \
        --exclude "src/adiff/hyperdual.f90" \
        --exclude "example" \
        --exclude "src/legacy/*" \
        --exclude "app"\
        --exclude "tools" \
        --fail-under-line 90
}

resumee() {
    [ $DID_TEST = 1 ] &&
        echo There has been $NAMING_ERRORS test naming errors

    if [ ${COVER} -le 90 ]; then
        echo "COVERAGE: " $(red $COVER)
    else
        echo "COVERAGE: " $(green $COVER)
    fi
}

case $1 in
    "install")  install_fpm;;
    "test") run_test;;
    "coverage") run_coverage;;
    *)
        run_test
        run_coverage
        resumee;;
esac

