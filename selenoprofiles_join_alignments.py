#!/usr/bin/python -u
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
#from profiles_classes import *

help_msg="""This program reads working directories of selenoprofiles and join the output profile alignments for the various species/targets found in the directory.
Usage:

selenoprofiles_join_alignments.py  -d selenoprofiles_results_folder   -o output_folder    -p profiles

"profiles" is a list of  comma-separated profile names, e.g. SelK,SelP
 if not provided, the program searches for .ali files
 if provided, you can restrict the list of  species/targets
-t        list of comma separated targets, each one in the form species_name.target_name
-s        list of species name, comma separated

##### other options  
-f          do no transfer the alignments: generate an unaligned fasta file with all predictions and the profile sequences
-i          file with list of .ali files to parse. Overrides -d and -p
-ds         do not attempt to shrink alignment at the end by realigning columns with lots of gaps with mafft
-suffix     suffix of alignments searched: instead of family.ali, looks for family.suffix.ali; suffix is used also for output files
-debug      debug mode; print file attempts
-p_list  | -s_list  | -t_list       = like -p | -s | -t, but read the list from a file, one per line
-h OR --help    print this help and exit """

command_line_synonyms={'p':'fam', 'p_list':'fam_list', 'P':'fam', 'profile':'fam', 's':'species', 's_list':'species_list', 't':'target', 't_list':'target_list'}

def_opt= {'temp':'/home/mmariotti/temp', 
'i':0, 'f':0, 'ds':0,
'o':'',
'v':0, 'Q':0, 
'fam':'', 'fam_list':'', 
'species':'', 'species_list':'', 
'target':'', 'target_list':'', 'd':'', 'suffix':'' , 'debug':0,
}


#########################################################
###### start main program function

def is_selenoprofiles_output_title(title):
  return 'chromosome:' in title and 'target:' in title and 'positions:' in title and 'strand:' in title

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'do', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
 
  #checking input
  working_directory='./' 
  if opt['d']: 
    check_directory_presence(opt['d'], 'working directory')
    working_directory=Folder(opt['d'])
    
  output_folder='./'
  if opt['o']:
    output_folder=Folder(opt['o'])
    test_writeable_folder(output_folder, 'output folder')
    
  i_handler=None
  if not opt['i'] and not opt['fam'] and not opt['fam_list']:
    find_cmnd="""find $(find """+working_directory+""" -mindepth 2  -maxdepth 2 -type d -name "output" )  -name "*.ali" """
    print 'Profile list not specified; searching all those in: '+working_directory
    i_handler=bbash(find_cmnd).split('\n')
    if i_handler==['']: i_handler=[]
  elif opt['i']:
    i_handler=open(opt['i'])

  if not i_handler:
    #determining families
    families_list=[]
    if opt['fam_list']:
      check_file_presence(opt['fam_list'], 'fam_list file')
      for f in open(opt['fam_list'], 'r').readlines():
        if f[:-1]:      families_list.append(f[:-1])
    elif opt['fam']:     families_list=opt['fam'].split(',')
    
    #determining targets, either by specified targets name or species name, or just trying all subfolders
    target_list=[]; species_list=[]
    if opt['target_list']:
      check_file_presence(opt['target_list'], 'target_list file')
      for f in open(opt['target_list'], 'r').readlines():
        if f[:-1]:      target_list.append(f[:-1])
    elif opt['target']:     target_list=opt['target'].split(',')

    if not target_list:
      all_subdirectories= bbash('find '+working_directory+' -maxdepth 1 -type d ').split('\n')    

      if opt['species_list']:  
        check_file_presence(opt['species_list'], 'species_list file')
        for f in open(opt['species_list'], 'r').readlines():
          if f[:-1]:      species_list.append(f[:-1])
      elif opt['species']:    species_list=opt['species'].split(',')

      if not species_list: target_list=all_subdirectories
      else: 
        for d in all_subdirectories:
          if base_filename(d).split('.')[0] in species_list:  target_list.append(d)

    if not target_list:     raise Exception, "no target found ..."
  
    #correcting to have no final "/" in target_list
    a=[]
    for t in target_list:
      if t[-1]=='/': t=t[:-1]
      a.append(t)
    target_list=a

    def next_attempt(fam):
      for target in target_list:
        file_attempt=target+'/output/'+fam+'.ali'
        if opt['suffix']:    file_attempt=target+'/output/'+fam+'.'+opt['suffix']+'.ali' #using suffix  
        yield file_attempt

  else:
    fams_list={}
    for line in i_handler:
      afile=line.rstrip()
      if is_file(afile):
        fams_list.setdefault(base_filename(afile).split('.')[0], []).append(afile)
      else: print 'ERROR didn\'t find: '+afile

    families_list=fams_list.keys()
    def next_attempt(fam):
      for f in fams_list[fam]:
        yield f

  ################### file list (attempts) is now defined. loading alignments and transferring them to get a joined ali
  something_found=False
  for fam in families_list:
    joined_ali=False
    for file_attempt in next_attempt(fam):

      if opt['debug']: print 'ATTEMPT: '+file_attempt
      if is_file(file_attempt):
        print 'READING: '+file_attempt
        this_ali= alignment(file_attempt); corrected_ali=alignment()
        for title in this_ali.titles():
          if is_selenoprofiles_output_title( title ):
    
            g=gene(); g.load_from_header(title)

            if g.species:      
              new_title=g.id+'.'+replace_chars(mask_characters(g.species), ' ', '_') +'.'+join(base_filename(g.target).split('.')[:-1], '.')
            else: new_title=title
            dont_print_until_double_commas=False
            for i in title.split()[1:]:
              if dont_print_until_double_commas:                dont_print_until_double_commas=  i[-1]!='"'
              else:
                if not i.startswith('species:'):  new_title+=' '+i
                elif i.startswith('species:"'):   dont_print_until_double_commas=True
            
          else:        new_title=title
          corrected_ali.add(new_title, this_ali.seq_of(title))
        
        if opt['f']:
          if not joined_ali: joined_ali=alignment()
          for new_title in corrected_ali.titles():               
            if not joined_ali.has_title(new_title):             joined_ali.add(   new_title, nogap(corrected_ali.seq_of(new_title))   )
          
        else:   #### normal case: now transfering the alignment
          if joined_ali:          joined_ali= joined_ali.transfer_alignment(corrected_ali, dont_shrink=True)
          else:                   joined_ali= corrected_ali

    if joined_ali:
      outfile=output_folder+fam+'.ali'
      if opt['suffix']: outfile=output_folder+fam+'.'+opt['suffix']+'.ali'
      if not opt['f'] and not opt['ds']:   joined_ali.shrink()
      print 'WRITING ---------> '+outfile
      joined_ali.display(outfile)
      something_found=True

  if not something_found: print '... nothing was found. See -help'

#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
