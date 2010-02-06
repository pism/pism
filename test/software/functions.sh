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

# This function is used to allow suppressing the output from a command when
# the test is run in the batch mode, but not when it is run manually.

# Note that it uses the $MPIDO variable; use "run -n NN command" with NN =
# number of processors if a command has to run under MPI and "run command" if
# you need to run it without MPI ("command" can not be "-n" in the latter
# case).
run () {
    COMMAND="echo ERROR: function run: no arguments!"
    PROGRAM="echo"
    # Find out if we need to use $MPIDO:
    if [ $1 == "-n" ]; then
	COMMAND="$MPIDO $@"
	PROGRAM=$3
    else
	COMMAND=$@
	PROGRAM=$1
    fi

    if [ -z `which $PROGRAM` ]; then
	echo ERROR: $PROGRAM not found!
    fi

    # Figure out if we need to dump the output:
    if [ $SILENT ]; then
	$COMMAND &> /dev/null
    else 
	echo "Running $COMMAND"
	$COMMAND
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

