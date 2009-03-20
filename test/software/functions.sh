# Use mpiexec if the MPIDO environment variable if empty:
if [ -z $MPIDO ];
then
    export MPIDO=mpiexec
fi

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
    echo "$@"
    echo "$@" >> $SUMMARY
}

# run a command and throw away the output if $SILENT is set

# This function is used to be able to suppress the output from a command when
# the test is run in the batch mode, but not when it is run manually.

# Note that it uses the $MPIDO variable; use "run -n NN command" with NN =
# number of processors if a command has to run under MPI and "run command" if
# you need to run it without MPI ("command" can not be "-n" in the latter
# case).
run () {

    # Figure out if we need to dump the output:
    if [ -z $SILENT ]; then

	# Print a message explaining what will happen:
	if [ $1 == "-n" ]; then
	    echo Running $MPIDO $@
	    $MPIDO $@
	else
	    echo Running $@
	    $@
	fi
    else
	if [ -z $VERBOSE ]; then
	    if [ $1 == "-n" ];
	    then
		$MPIDO $@ &> /dev/null
	    else
		$@ &> /dev/null
	    fi
	else
	    if [ $1 == "-n" ];
	    then
		print "Running $MPIDO $@ &> /dev/null"
		$MPIDO $@ &> /dev/null
	    else
		print "Running $@ &> /dev/null"
		$@ &> /dev/null
	    fi
	fi
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

