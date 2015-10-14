#!/usr/bin/python -u
from string import *
import sys
import commands 
import os

help_msg="""Simple program to test selenoprofiles. A couple of fasta files included in the test package are necessary for this program to be working.

### Options:
-b                        path to the Selenoprofiles executable. If not provided, it is determined with a "which" command
-print_opt                print currently active options
-h OR --help              print this help and exit
"""


#########################################################
###### start main program function

import traceback 
class notracebackException(Exception):   """ just a puppet class to show a cleaner output in case of errors"""
is_file=os.path.isfile
def bash(command, print_it=0):
  """Utility to run bash commands. a tuple (exit_status, message) is returned, where message includes both std input and stderr. If argument print_it==1  the command is printed before execution."""
  if print_it:    print command
  return commands.getstatusoutput(command)
def printerr(msg, newline):
  if newline: msg=str(msg)+'\n'
  print >> sys.stderr, str(msg)

def main(args={}):
#########################################################
############ loading options

  if '-h' in sys.argv or '-help' in sys.argv or '--help' in sys.argv:
    print help_msg
    sys.exit()

  print "Looking for Selenoprofiles executable..."

  # determining path to executable
  if '-b' in sys.argv:    
    selenoprofiles_bin =    sys.argv[  sys.argv.index('-b') +1 ]
  else: 
    b=bash('which Selenoprofiles')
    if not b[0]:  selenoprofiles_bin = b[1]
    else:     
      b=bash('which selenoprofiles_3.py')
      if not b[0]:  selenoprofiles_bin = b[1]
      else:         raise notracebackException, "ERROR selenoprofiles was not found! Please provide the executable with option -b "
    if not selenoprofiles_bin or not is_file(selenoprofiles_bin):
      raise notracebackException, "ERROR selenoprofiles was not found: "+selenoprofiles_bin+"    Please provide the executable with option -b "

  results_folder='test_results'
  sequence_file='test_sequences.fa'

  if not is_file(sequence_file):
    raise notracebackException, "ERROR cannot run tests: test_sequences.fa was not found! Please run this script in selenoprofiles_3 installation directory"

  print "Starting tests..."

  ## running first test 
  print '\n## Test 1:   SelR profile (should take ~1 min or less)' 
  b=  bash (selenoprofiles_bin+' '+  results_folder +' '+sequence_file+' -s "Homo sapiens" -P SelR -B -output_p2g  > log_test1 2>&1', 1)
  if b[0]:  raise notracebackException, "Selenoprofiles ERROR: test failed. see error below, or in log_test1\n\n"+b[1]
  
  b=bash('grep ERROR log_test1')
  if not b[0]:  raise notracebackException,  "-- ERROR!\n"+b[1]

  b=bash('grep WARNING log_test1')
  if not b[0]:  print  "-- WARNING!\n"+b[1]
  b=bash('grep P2G log_test1 | grep cysteine ')
  if b[0]:   raise notracebackException, "-- ERROR! Although selenoprofiles didn't crash, it didn't find what was supposed to."
 
  print "##### OK!"  


  ## running second test 
  print '\n## Test 2:   DI profile (should take ~1 min or less)' 
  b=  bash (selenoprofiles_bin+' '+  results_folder +' '+sequence_file+' -s "Homo sapiens" -P DI -B -output_p2g > log_test2 2>&1', 1)
  if b[0]:  raise notracebackException, "Selenoprofiles ERROR: test failed. see error below, or in log_test2\n\n"+b[1]

  b=bash('grep ERROR log_test2')
  if not b[0]:  raise notracebackException,  "-- ERROR!\n"+b[1]

  b=bash('grep WARNING log_test2')
  if not b[0]:  print  "-- WARNING!\n"+b[1]
  b=bash('grep P2G log_test2 | grep selenocysteine ')
  if b[0]:   raise notracebackException, "-- ERROR! Although selenoprofiles didn't crash, it didn't find what was supposed to."
  b=bash('grep SECISearch log_test1')
  if not b[0]:  print  "-- WARNING! SECISearch did not work correctly. Ignore this if you're using selenoprofiles for non-selenoprotein families."
  print "##### OK!"  

  ## running third test 
  print '\n## Test 3:   SPS profile (should take ~3 min or less)' 
  b=  bash (selenoprofiles_bin+' '+  results_folder +' '+sequence_file+' -s "Homo sapiens" -P SPS -B -output_p2g > log_test3 2>&1', 1)
  gb= bash('grep ERROR log_test3')
  if not gb[0]:     raise notracebackException, "ERROR the tag_score system is not installed: you need the NR protein database in your system. You can safely ignore this if you're not going to use tag_score based filtering, and neither the built-in profiles (it is necessary for selenoprotein search). See complete error below:\n\n"+gb[1]
  elif b[0]:         raise notracebackException, 'unknown ERROR: please inspect log_test3'    
  print "##### OK!"  


  ## running fourth test 
  print '\n## Test 4:   GPx profile (should take ~3 min or less)' 
  b=  bash (selenoprofiles_bin+' '+  results_folder +' '+sequence_file+' -s "Homo sapiens" -P GPx -B -output_p2g > log_test4 2>&1', 1)
  gb= bash('grep ERROR log_test4')
  if not gb[0]:     raise notracebackException, "ERROR gene_ontology utilities are not installed. You can safely ignore this if you're not going to use gene_ontology based filtering, and neither the built-in profiles (it is necessary for selenoprotein search). See complete error below:\n\n"+gb[1]
  elif b[0]:         raise notracebackException, 'unknown ERROR: please inspect log_test4'
  print "##### OK!"  

  


  ###############



#######################################################################################################################################


def close_program():

  if sys.exc_info()[0]: #an exception was raised. let's print the information using printerr, which puts it in the logfile as well, if any.   
    if issubclass(sys.exc_info()[0], notracebackException):      printerr( sys.exc_info()[1], 1)
    elif issubclass(sys.exc_info()[0], SystemExit):      pass
    else:                                                                  printerr('ERROR '+ traceback.format_exc( sys.exc_info()[2]) , 1)
    
  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except:
    close_program()

