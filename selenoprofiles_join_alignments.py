#!/usr/bin/python -u
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
#from profiles_classes import *

help_msg="""This program reads working directories of selenoprofiles and join the output profile alignments for the various species/targets found in the directory.

With no options, the program assumes we are in a working directory and tries to find alignments for all builtin families in all subdirectories (i.e. for all targets).
Options can be specified to decide which alignments to search for. All options require arguments except -D

-d        selenoprofiles result directory
-o        output directory

##### options to determine input profile alignments
-i        file with list of .ali files to parse, one per line. If this is active, none of the next six options described should be active

-p        list of comma separated profile names that the program will search for   #e.g. sps,SelK,SelP
-p_list   file with a list of profile names, one per line. Overrides -p

-t        list of comma separated targets. Each one must be in the form species_name.target_name, exactly as the folder name created by selenoprofiles.
-t_list   file with list of targets, one per line. Overrides -t

# If any of the 2 previous options is specified, the list of target is determined and these will be the only subfolders where the program looks for alignments. If instead any of the 2 next options are specified, the species list is determined, which means that any target associated to a desired species will be accepted.

-s        list of species name, comma separated. Careful! species names containing commas will not be read correctly: use -species_list instead.
-s_list   file with list of species, one per line. Overrides -species

-suffix   suffix of alignments searched. By default, this programs searches for alignments named like family.ali; with this option, it will look for files name like family.suffix.ali; The suffix will be used also for the files in output

##### other options  
-D              don't change the names of selenoprofiles predictions. Normally the species name is added to allow discriminating between results on different species.
-debug          debug mode; print file attempts
-f              do no transfer the alignments: generate an unaligned fasta file with all predictions and the profile sequences
-ds             do not attempt to shrink the alignment at the end. Normally, this procedure detects the columns which are misaligned due to the transfering alignment procedures, and attempts to refine these regions using mafft
-sp_config      path to selenoprofiles config file. This is used only to get the default list of profiles, and is used only if no option is active among -p, -p_list, -target, -target
-h OR --help    print this help and exit """

command_line_synonyms={'p':'fam', 'p_list':'fam_list', 'P':'fam', 'profile':'fam', 's':'species', 's_list':'species_list', 't':'target', 't_list':'target_list'}

def_opt= {'temp':'/home/mmariotti/temp', 
'i':0, 'f':0, 'ds':0,
'o':'',
'v':0, 'Q':0, 
'fam':'', 'fam_list':'', 
'species':'', 'species_list':'', 
'sp_config': '/users/rg/mmariotti/scripts/selenoprofiles_3.config',
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
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input
  
  try: 
    ### reading default list of profiles from selenoprofiles main configuration file
    sp_config_hash= configuration_file(opt['sp_config'])
    families_sets={}
    if sp_config_hash.has_key('families_set'):
      for k in sp_config_hash['families_set']:
        families_sets[ k ] = [del_white(word) for word in sp_config_hash['families_set'][k].split(',')]
    profile_names=sp_config_hash['profile'].split(',')
    index_profile=0
    while any ([families_sets.has_key(p) for p in profile_names]):
      this_profile=profile_names[index_profile]
      while families_sets.has_key(this_profile):
        profile_names[index_profile:index_profile+1]= families_sets[this_profile]
        this_profile=profile_names[index_profile]
      index_profile+=1
    # removing redundancy in profiles names
    profiles_names_length=len(profile_names)
    for i in range(profiles_names_length-1, -1, -1):
      profile_name=profile_names[i]
      if profile_name in profile_names[:i]:
        profile_names.pop(i)  
    default_families_list=profile_names
  except: 
    printerr('WARNING could not determine default list of profiles from configuration file: '+str(opt['sp_config'])+' ;  this is not a problem if the list of profiles is specified on the command line', 1)

  working_directory='./' 
  if opt['d']: 
    check_directory_presence(opt['d'], 'working directory')
    working_directory=Folder(opt['d'])
    
  output_folder='./'
  if opt['o']:
    output_folder=Folder(opt['o'])
    test_writeable_folder(output_folder, 'output folder')
    
  if not opt['i']:
    #determining families
    families_list=[]
    if opt['fam_list']:
      check_file_presence(opt['fam_list'], 'fam_list file')
      for f in open(opt['fam_list'], 'r').readlines():
        if f[:-1]:      families_list.append(f[:-1])
    elif opt['fam']:     families_list=opt['fam'].split(',')
    else:
      families_list= default_families_list
    
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
    for line in open(opt['i']):
      afile=line.rstrip()
      if is_file(afile):
        fams_list.setdefault(base_filename(afile).split('.')[0], []).append(afile)
      else: print 'ERROR didn\'t find: '+afile

    families_list=fams_list.keys()
    def next_attempt(fam):
      for f in fams_list[fam]:
        yield f

  ################### file list (attempts) is now defined. loading alignments and transferring them to get a joined ali
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
