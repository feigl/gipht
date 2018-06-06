# splitp.awk --- do split in awk based on pattern
#
# based on split.awk function from
# http://www.gnu.org/manual/gawk/gawk.html#Split-Program
# usage: split [-num] [file] [outname]

BEGIN {
    outfile = "pair."    # default
	count = 1000
	if (ARGC > 4)
	    usage()
		
		i = 1
		if (ARGV[i] ~ /^-[0-9]+$/) {
		    count = -ARGV[i]
		    ARGV[i] = ""
		    i++
		}
# test argv in case reading from stdin instead of file
    if (i in ARGV)
	i++    # skip data file name
	    if (i in ARGV) {
		outfile = ARGV[i]
		    ARGV[i] = ""
		    }
    
#s1 = s2 = s3 = "0"
#    s1 = 1
# 2012-JUL-05 Kurt
    s1 = 0
	out = sprintf("pair%03d.lst",s1)
#out = (outfile s1 s2)
	}


#The next rule does most of the work. 
#tcount (temporary count) tracks how many lines have been printed to the output file so far. 
#If it is greater than count, it is time to close the current file and start a new one. 
#s1 and s2 track the current suffixes for the file name. If they are both ¡z¢, the file is just too big. 
#Otherwise, s1 moves to the next letter in the alphabet and s2 starts over again at ¡a¢:


{
    if ($0 ~ /a/ ) {
#      if (++tcount > count) {
	close(out)
	    s1 = s1+1
	    out = sprintf("pair%03d.lst",s1)
	    }
    print > out
	}


#The usage function simply prints an error message and exits:


function usage(   e)
{
    e = "usage: split [-num] [file] [outname]"
	print e > "/dev/stderr"
	exit 1
	}


#The variable e is used so that the function fits nicely on the page.
#This program is a bit sloppy; it relies on awk to automatically 
#close the last file instead of doing it in an END rule. 
#It also assumes that letters are contiguous in the character set, which isn't true for EBCDIC systems.

