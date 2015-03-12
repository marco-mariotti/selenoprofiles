#!/soft/bin/python -u
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
sys.path.append('/users/rg/mmariotti/selenoprofiles/trunk')
from MMlib import *
from selenoprofiles_3 import selenoprofiles_db, load, set_selenoprofiles_var

help_msg="""Utility to manipulate the sqlite results database created by selenoprofiles_3.
Mostly, its function is to fix entries resulting from executions that went wrong, but it can be used for fast output as well

Usage:    selenoprofiles_database.py      some_path/results.sqlite   [options]

### Operations: inspect or modify the database
-list                     list the profiles annotated for this target in the database and their state
-clean                    tell the database that no profile is currently under computation
-remove_overlaps          remove all overlaps recorded in the database
-remove    +              remove all results for a profile given as argument
-remove_list +            same as -remove, but provide a list of profile arguments in a file, one per line

## Fast output: specify any of these option to get fast output
-fasta                    print peptide sequence in fasta of all results (with label "kept") to standard output
-cds                      print coding sequence in fasta of all results (with label "kept") to standard output
-gff                      print genomic coordinates in gff of all results (label "kept"). NOTE: features associated to this result will not be present in the gff. If you're using custom features, you should output through the standard procedure, with selenoprofiles_3.py

-filter                   get fast output only for results passing these sql evaluation. This is directly concanenated to the sql query, so be careful. Check the results table structure to formulate a filter. Default is ' state=="kept" '
-p                        prompt: open the python interactive environment

-print_opt                print currently active options
-h OR --help              print this help and exit"""

command_line_synonyms={}

def_opt= { 'temp':'/home/mmariotti/temp', 
'i':0, 'p':0,
'clean':0, 'list':0, 'remove_overlaps':0, 'remove':None, 'remove_list':0, 
'sp_config': '/users/rg/mmariotti/scripts/selenoprofiles_3.config',
'fasta':0, 'cds':0, 'gff':0,  
'filter': 'state=="kept"',
'v':0,
}


#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input
  global input_file;   input_file=opt['i'];   check_file_presence(input_file, 'results.sqlite file')
  config_filename=opt['sp_config']

  db=selenoprofiles_db(input_file)

  
  try:   load(config_filename, {'puppet_option':1}) #putting a fake option just to give a non-null hash as args; this will avoid selenoprofiles trying to attempt loading options form commandline
  except: pass #the load function crashes since it doesn't find the file necessary for run, for example the target file. nonetheless the variables are in place.

  set_selenoprofiles_var('max_attempts_database', 0)

  if opt['clean']:
    """ to set the property "has_results_for_profile" as  "UNCHECKED-2" for entries with "UNCHECKED" or "WAITING", and "NO" for results with "ONGOING".  """
    write('Cleaning database... ', 1)
    db_cursor=db.cursor()
    db_cursor.execute('UPDATE has_results_for_profile SET  has_output="UNCHECKED-2" WHERE has_output == "UNCHECKED" OR  has_output == "WAITING" ')
    db_cursor.execute('UPDATE has_results_for_profile SET  has_output="NO" WHERE has_output == "ONGOING"')
    db.save()
    write('--finished cleaning.', 1)

  if opt['remove_overlaps']:
    """ change the state attribute of all results labelled as "overlapping"  """
    db.remove_overlapping_states(silent=False)

  profile_list_to_remove=None
  if opt['remove_list']:           profile_list_to_remove=[ line.strip() for line in open(opt['remove_list']) if line.strip() ] 
  if not opt['remove'] is None:    profile_list_to_remove=[ opt['remove'] ]
  
  if profile_list_to_remove:
    """ Removing entries from the table results and also has_results_for_profile"""
    write('Removing entries for '+str(len(profile_list_to_remove))+' profile'+'s'*int( len(profile_list_to_remove)>1   )+'... ', 1 )
    n_profile=0    
    for profile_name in profile_list_to_remove:
      db_cursor=db.cursor()
      db_cursor.execute( 'SELECT id FROM results WHERE profile=="'+profile_name+'"' ) 
      for resultid in db_cursor.fetchall():         db_cursor.execute('DELETE FROM features  WHERE resultid=="'+str(resultid)+'"')
      db_cursor.execute('DELETE FROM has_results_for_profile  WHERE profile=="'+profile_name+'"')
      db_cursor.execute('DELETE FROM results  WHERE profile=="'+profile_name+'"')
      n_profile+=1
      if len(profile_list_to_remove)>1: 
        service(str(n_profile)+' profiles removed')        
    db.save()        
    write('-- finished removing entries.', 1)

  if opt['list']:
    """ print table has_results_for_profile"""
    db_cursor=db.cursor()    
    db_cursor.execute( 'SELECT * FROM has_results_for_profile' ) 
    for obj in db_cursor.fetchall():
      write(str(obj[1]).ljust(40)+' '+str(obj[2]), 1)

  if opt['cds'] or opt['fasta'] or opt['gff']:
    info_target_file= directory_name(input_file)+'/info_target.txt'
    info_target_hash= configuration_file(info_target_file)
    target_file= info_target_hash['Target file']
    species= info_target_hash['Species name']

    db_cursor=db.cursor()    
    db_cursor.execute( 'SELECT profile, program, target_header, chromosome, alignment_target, cds_sequence, label FROM results WHERE '+opt['filter'] ) 
    obj=db_cursor.fetchone()
    while obj:
      profile, program, target_header, chromosome, alignment_target, cds_sequence, label= obj
      pred_index, strand_field, pos_field = target_header.split()
      fasta_header=profile+'.'+pred_index+'.'+label+' chromosome:'+chromosome+' '+strand_field+' '+pos_field+' species:"'+species+'" target:'+target_file+' prediction_program:'+program
      if opt['fasta'] or opt['cds']:
        write(">"+fasta_header , 1)
        if opt['fasta']: write( fasta(nogap(alignment_target)), 1)
        if opt['cds']:   write( fasta(nogap( cds_sequence  )), 1)
      if opt['gff']: 
        g=gene();         g.load_from_header(  fasta_header )
        write(g.gff(program="selenoprofiles_"+program), 1 )  
            
      obj=db_cursor.fetchone()

  if opt['p']:
    db_cursor=db.cursor()        
    interactive_mode(message='Variables you may be interested into: \ndb         -- class: selenoprofiles_db\ndb_cursor  -- class: pysqlite cursor')()


    
  ###############



#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass

  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
