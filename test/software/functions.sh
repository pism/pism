if [ -z $SUMMARY ];
then
    SUMMARY=summary.txt
fi

if [ -z $PASSCOUNT ];
then
    PASSCOUNT=0
fi

if [ -z $FAILCOUNT ];
then
    FAILCOUNT=0
fi

# Print both to the standard out and in the summary file
print () {
    echo "$1"
    echo "$1" >> $SUMMARY
}

# run a command and throw away the output if $SILENT is set

# This function is used to be able to suppress the output from a command when
# the test is run in the batch mode, but not when it is run manually.
run () {
    if [ -z $SILENT ]; then
	echo Running $@
	"$@"
    else
	"$@" &> /dev/null
    fi
}

# print a message and count a success
pass () {
    print "PASS: $test"
    PASSCOUNT=$(( $PASSCOUNT + 1 ))
    cleanup
}

# print a message and count a failure
fail () {
    print "FAIL: $test (Reason: $1)"
    FAILCOUNT=$(( $FAILCOUNT + 1 ))
    files=""			# so that the next test script does not delete them
}

cleanup () {
    cd "$dir"
    rm -f *~ $files
}

