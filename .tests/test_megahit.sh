# Fix for #167:
# Check that if megahit gives a nonzero exit code it is handled appropriately.
# The two main cases are 255 (empty contigs) and anything else nonzero
# (presumed to be memory-related in the assembly rules).
# Checking for successful behavior is already handled in test_all.
function test_assembly_failures {
    # Up to just before the assembly rules, things should work fine.
    sunbeam run --profile $TEMPDIR/ all_decontam --configfile=$TEMPDIR/tmp_config.yml
    # Remove previous assembly files, if they exist.
    rm -rf $TEMPDIR/sunbeam_output/assembly

    # If megahit exits with 255, it implies no contigs were built.
#    mkdir -p "$TEMPDIR/megahit_255"
#    echo -e '#!/usr/bin/env bash\nexit 255' > $TEMPDIR/megahit_255/megahit
#    chmod +x $TEMPDIR/megahit_255/megahit
#    (
#    export PATH="$TEMPDIR/megahit_255:$PATH"
#    txt=$(sunbeam run -- --configfile=$TEMPDIR/tmp_config.yml -p all_assembly)
#    echo "$txt" > /mnt/d/Penn/sunbeam/log.txt
#    echo "$txt" | grep "Empty contigs"
#    )

    # If megahit gives an exit code != 0 and != 255 it is an error.
    mkdir -p "$TEMPDIR/megahit_137"
    echo -e '#!/usr/bin/env bash\nexit 137' > $TEMPDIR/megahit_137/megahit
    chmod +x $TEMPDIR/megahit_137/megahit
    (
    export PATH="$TEMPDIR/megahit_137:$PATH"
    # (This command should *not* exit successfully.)
    ! txt=$(sunbeam run --profile $TEMPDIR/ all_assembly --configfile=$TEMPDIR/tmp_config.yml)
    echo "$txt" | grep "Check your memory"
    )
}