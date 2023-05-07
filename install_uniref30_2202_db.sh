#!/bin/bash -ex
# shell script adapted from https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh to install only the uniref30_2202_db
ARIA_NUM_CONN=8
WORKDIR="${1:-$(pwd)}"

cd "${WORKDIR}"

hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
	    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && set -e && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && set -e && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && set -e && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

if [ ! -f UNIREF30_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2202.tar.gz" "uniref30_2202.tar.gz"
  tar xzvf "uniref30_2202.tar.gz"
  mmseqs tsv2exprofiledb "uniref30_2202" "uniref30_2202_db"
  mmseqs createindex "uniref30_2202_db" tmp1 --remove-tmp-files 1
  if [ -e uniref30_2202_db_mapping ]; then
    ln -sf uniref30_2202_db_mapping uniref30_2202_db.idx_mapping
  fi
  if [ -e uniref30_2202_db_taxonomy ]; then
    ln -sf uniref30_2202_db_taxonomy uniref30_2202_db.idx_taxonomy
  fi
  touch UNIREF30_READY
fi