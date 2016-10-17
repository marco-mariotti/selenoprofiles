#!/usr/bin/python
__author__  = "Marco Mariotti"
__email__   = "marco.mariotti@crg.eu"
__licence__ = "GPLv3"

help_msg="""This is the Selenoprofiles3 installation script. Please move the Selenoprofiles installation folder to its final destination, then cd into it and run this script with the option specified below. Before installing, take care to have installed and available in bash the following executables:
. python (2.6 <= Version < 3.0), gawk
. blastall, blastpgp, formatdb    (NCBI blastall package)
. exonerate, fastafetch, fastaindex, fastasubseq, fastalength, fastarevcomp   (exonerate package)
. genewise   (wise2 package)
. mafft

For help on installing these programs, check page: http://big.crg.cat/news/20110616/installing_programs_and_modules_needed_by_selenoprofiles
The installation script will check that the slave programs used by selenoprofiles are installed. The command "which" will be used determine which executables to use, and they will be linked inside a bin/ folder in the installation folder. If you want to use executables which are not those returned by a "which" command, the best option is to change the PATH bash variable temporarly so "which" will point to the desired executable, and run the installation. Afterwards, you may return to your original PATH settings.
This installation script will set up most variables in the selenoprofiles configuration file. One important option is the temporary folder used:
-temp               temporary folder later used by selenoprofiles (default: /tmp)

## Minimal installation:  
This is suitable for searching user-provided profiles. Tag-blast based filtering will not be available. To perform this installation run:
$ python  install_selenoprofiles.py -min

## Complete installation:  
This allows to use more advanced forms of filtering (tag-score and go-score), which are required to search for the selenoprotein and Sec machinery profiles that come with selenoprofiles. Run:
$ python  install_selenoprofiles.py -full
This installation requires the NCBI nr protein database, and two gene ontology related files. The installation script will attempt to fetch them through internet. If you have them already on your system, you can link them with the following options:
-nrdb               path to the "nr" fasta file as downloaded (and uncompressed) from ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz 
-formatted_nrdb     path to the formatted files produced by formatdb on the nr file. This option is useful only if these files are in a location different from the one provided with -nrdb. The common prefix should be provided; e.g. if formatted files are called like /path/nr.00.phd /path/nr.00.phi ... /path/nr.04.psq, you have to provide just /path/nr 
-godb               path to the gene_ontology_ext.obo file as downloaded from http://www.geneontology.org/ontology/obo_format_1_2/
-gomap              path to the idmapping_GI_GO.tab file as downloaded from http://genome.crg.es/~mmariotti/idmapping_GI_GO.tab.gz 

## Updating: 
If a prior Selenoprofiles installation is detected ("selenoprofiles_3" or "Selenoprofiles" are searched in your PATH variable), the script will ask whether you want to update it. If you want to avoid updating, add option -dont_update to the command line. To force updating without prompt, use -update. In any case, please use -min or -full depending on your type of installation. 

## Other options: 
-n_cpus             number of cpu cores used for blast searches; it is passed to the blastall command through its "-a" option 
-SS                 use this option if SECISearch3 is installed in your system, and it will be activated in Selenoprofiles. SECISearch3 is not distributed yet, but is available on request. However, for simplicity we recommend to run Selenoprofiles without installing SECISearch3, and then search the three prime sequences of candidates using the webserver at http://seblastian.crg.es/
-bSS                same for bSECISearch (bSeblastian), a tool under development for bacterial selenoproteins 

If you encounter any problem during installation, contact marco.mariotti at crg.es"""

class notracebackexception(Exception):
  """ to raise exceptions without printing the usual verbose -tree like report (see end of script)"""

try:  
  from MMlib import *
  def_opt={ 'temp':'/tmp/', 'no_nr':0, 'no_go':0,  'SS':0, 'bSS':0, 'nrdb':0, 'formatted_nrdb':0,  'dont_update':0, 'godb':0, 'gomap':0, 'min':0, 'n_cpus':1, 'full':0, 'update':0 }
  opt=command_line(def_opt, help_msg, '' )

  if opt['formatted_nrdb'] and not opt['nrdb']: 
    print "ERROR option -formatted_nrdb must be used together with option -nrdb"
    sys.exit(8)

  if opt['min']:   
    opt['no_nr']=1; opt['no_go']=1
    installation_type='min'
  elif opt['full']: 
    installation_type='full'
  else:
    raise notracebackexception,  "ERROR you must either provide -min or -full to perform a minimal or a full installation. Use --help for information."
  
  print "\n##########################  Selenoprofiles3 installation script  ##########################"
  installation_directory=Folder(os.path.abspath('.'))
  temp_directory=Folder(opt['temp'])
  updating=False

  ## looking for a previous installation of selenoprofiles. 
  if not opt['dont_update']:
    for p_name in ["Selenoprofiles", "selenoprofiles_3", "selenoprofiles_3.py"]:
      b=bash('which '+p_name)
      if not b[0]:  break # program found 
    if not b[0]: 
      #### installation found!
      previous_installation_dir= directory_name(   dereference(b[1])   )+'/'

      if  opt['update']: 
        updating=True
      else:
        possible_update_choices='YNQ'
        update_choice='uNmAtChAbLe'
        while (not update_choice) or  (not update_choice in possible_update_choices) :
          if not update_choice =='uNmAtChAbLe': print '#### Your choice was not recognized... please type it again!'
          update_choice= upper( raw_input('\nAn existing selenoprofiles installation found in '+previous_installation_dir+' ; do you want to update that installation ?\nType:  Y to update,   N to perform a new installation,   Q (or Ctrl-C) to exit\n>') )

        if   update_choice=='Y': updating=True
        elif update_choice=='N': updating=False
        elif update_choice=='Q': raise notracebackexception, 'Installation aborted'
        
      if updating:         
        print "--- Updating previous selenoprofiles installation found in "+previous_installation_dir+'...'
        test_writeable_folder(previous_installation_dir, 'ERROR with permissions: can\'t write in the installation directory' )      
        updating=True
        for file_or_folder in ['CHANGES.txt', 'MMlib.py', 'README' , 'blaster_parser.g' , 'profiles', 'selenoprofiles_3.py', 'selenoprofiles_build_profile.py', 'selenoprofiles_join_alignments.py', 'selenoprofiles_database.py',  'test_selenoprofiles.py', 'test_sequences.fa', 'selenoprofiles_tree_drawer.py', 'selenoprofiles_manual.pdf', 'version.txt', 'selenoprofiles_3.config']:
          print 'Updating file: '+file_or_folder
          if file_or_folder=='selenoprofiles_3.config':
            print 'NOTE: the configuration file has been updated, so your old settings are not currently used. You can recover your old settings in '+previous_installation_dir+'selenoprofiles_3.config.backup and paste them in the new configuration file: '+previous_installation_dir+'selenoprofiles_3.config'
            bash('cp '+previous_installation_dir+'selenoprofiles_3.config '+previous_installation_dir+'selenoprofiles_3.config.backup')
          bash('rm -fr '+previous_installation_dir+file_or_folder)
          bash('cp -fr '+installation_directory+file_or_folder+' '+previous_installation_dir)

        installation_directory=previous_installation_dir
    elif opt['update']: 
      raise notracebackexception, "ERROR option -update is active, but no existing Selenoprofiles executable was found! make sure it is available in terminal and rerun this installation script"
    
  
  confirm_choice='uNmAtChAbLe'; possible_confirm_choices='YQ'
  while (not confirm_choice) or (not confirm_choice in possible_confirm_choices):
    if confirm_choice != 'uNmAtChAbLe': print '#### Your choice was not recognized... please type it again!'
    confirm_choice= upper(raw_input('\nThis program will now install selenoprofiles_3 in this folder:    '+installation_directory+'\nIf you want to choose a different installation directory, quit the installation, move this folder to its final destination, then run again the installation script.\n\nCurrently you are using:\n installation type= '+installation_type+'\n -temp      '+temp_directory+'\n -n_cpus    '+str(opt['n_cpus'])+'      # for blast searches only \n\n Type Y and Enter to confirm,  Q or Ctrl-C to quit\n>'))
    if confirm_choice=='Q': raise notracebackexception, 'Installation aborted'

  test_writeable_folder(installation_directory, 'ERROR with permissions: can\'t write in the installation directory' )   
  wget_folder=Folder(installation_directory+'downloading_files')
  
  if not is_file(installation_directory+'selenoprofiles_3.py'): raise notracebackexception, "ERROR selenoprofiles_3.py file not found. install_selenoprofiles must be run inside the installation directory. "
  bin_folder=Folder(installation_directory+'bin')
  print "Selenoprofiles3 installation start!"

  #determining which program (interpreters) will be used. These will be written in the she bang lines of the selenoprofils scripts
  print 'Determining which python, perl and gawk executable to use ...'
  extension_to_program={'py':'python', 'g':'gawk'}
  program_to_executable={}
  for k in extension_to_program:
    program=extension_to_program[k]
    executable=bbash('which '+program)
    if program=='python': executable+=' -u'  
    if program=='gawk': executable+=' -f'
    program_to_executable[program]=executable
  #testin pysqlite2
  print 'Testing if module pysqlite is present for this python executable ...'
  b=bash(program_to_executable['python']+' -c "from sqlite3 import dbapi2 as sqlite"')
  if b[0]: raise notracebackexception, 'ERROR module sqlite3 not available in this python environment. Module sqlite3 is available in python 2.5 and newer. Please use python 2.6 with selenoprofiles3.'

  ## checking external programs
  programs_and_packages = [('blastall', "ncbi blast package 2.2.2x"), ('blastpgp', "ncbi blast package 2.2.2x"), ('formatdb', "ncbi blast package 2.2.2x"), ('exonerate', "exonerate version 2.0.0"),  ('fastafetch', "exonerate version 2.0.0"),  ('fastaindex', "exonerate version 2.0.0"),  ('fastasubseq', "exonerate version 2.0.0"),  ('fastalength', "exonerate version 2.0.0"),('fastarevcomp', "exonerate version 2.0.0"), ('genewise', "Wise2 package -- http://www.ebi.ac.uk/Tools/Wise2/"),('mafft', 'mafft from http://align.bmr.kyushu-u.ac.jp/mafft/software/')]
  error_msg=''
  for program, package_name  in  programs_and_packages:
    if not ( updating  and is_file(bin_folder+program) ): #happens when updating
      b=bash('which '+program)
      if b[0]:     error_msg+= 'ERROR can\'t find program '+program+' \tPlease install package '+package_name+' and make sure that program is available in your bash terminal.\n'
      else: 
        program_bin=b[1].split('\n')[0]
        if not os.path.isfile(installation_directory+program):
          print 'Linking '+program+' in '+bin_folder+' ... '
          bash('rm '+bin_folder+program) #deleting if necessary
          bbash('ln -fs '+program_bin+' '+bin_folder) #linking

  if error_msg:    raise notracebackexception, error_msg
      
  ### fetching files that will be used by selenoprofiles
  libraries_folder=Folder(installation_directory+'libraries')

  #checking/downloading nr
  if not opt['no_nr']:
    if not is_file(libraries_folder+'nr.fa'):
      if not opt['nrdb']:
        answer= raw_input('You did not specify the path to the nr database on your computer using option -nrdb. This will cause this script to try and download it through internet at the ncbi ftp site. The nr database is large (>3Gb) so this may take a long time. We recommend to avoid this if the nr database is already present on your computer.\nAre you sure to continue and start downloading nr? (Y/N)   ')
        if upper(answer) in ['YES', 'Y']:
          print 'Fetching nr database from ncbi (~3Gb compressed). This may take a long time ... '
          b=['','']
          try:
            b=bash('cd '+wget_folder+'; wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz && mv nr.gz '+libraries_folder);  assert not b[0]
          except: raise notracebackexception, "ERROR fetching nr database from ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz : "+b[1]+'\n\nPlease download it manually and provide the nr fasta file to this installation script with -nrdb ; if otherwise you want to skip the installation of the nr database use option -no_nr (or perform a minimal installation with -min)'
          print 'Uncompressing nr.gz ...'
          bbash('gunzip '+libraries_folder+'nr.gz')
          bash ('mv '+libraries_folder+'nr '+libraries_folder+'nr.fa ')
        else: raise notracebackexception, "Installation aborted"
      else:
        print 'Linking '+abspath(opt['nrdb'])+' in '+libraries_folder+' ...'
        bbash('ln -fs '+abspath(opt['nrdb'])+' '+libraries_folder+'nr.fa')
        #checking for the formatted files produced by formatdb. If they are present, they are linked so they won't be produced
        if not opt['formatted_nrdb']: opt['formatted_nrdb']=opt['nrdb']
        b=bash('ls -1 '+opt['formatted_nrdb']+'.*p[A-z][A-z] '+opt['formatted_nrdb']+'.*pal')
        if not b[0]: #found
          for ffile in b[1].split('\n'):
            #ffile name example:   /db/seq/databases/nr_uncompressed/nr.fa.02.psq
            link_name=  'nr.fa.'+   join( ffile.split('.')[-2:], '.')
            if ffile.endswith('pal'):               link_name=  'nr.fa.pal'
            bbash('ln -fs '+abspath(ffile)+' '+libraries_folder+link_name)
        #checking for the index file produced by fastaindex    ## if not found, it will be produced on the first use.
        for name_attempt in  [  opt['nrdb']+'.index', join(opt['nrdb'].split('.')[:-1], '.')+'.index']:
          if is_file(name_attempt):
            print "Index file found: "+name_attempt
            bbash('ln -fs '+abspath(name_attempt)+' '+libraries_folder+'nr.index')
            break
    else:
      print "File found: "+libraries_folder+'nr.fa'

    #check/produce formatted file on nr
    if not is_file(libraries_folder+'nr.index'): #index file is not already present
      print "Running fastaindex on the nr file in "+libraries_folder+'nr ...'
      bbash('fastaindex '+libraries_folder+'nr.fa '+libraries_folder+'nr.index')
    if bash('ls -1 '+libraries_folder+'nr.*p[A-z][A-z]')[0]: #formatted files are not already present
      print "Running formatdb (with options -p T -o T) on the nr file in "+libraries_folder+'nr ...'
      bbash('formatdb -i '+libraries_folder+'nr.fa -o T -p T')

  """ deprecated from version 3.2
  ## ncbi taxonomy database
  if not is_file(libraries_folder+'names.dmp'):
    if not opt['taxdb']:
      print 'Fetching taxonomy database from ncbi (~18Mb) ... '
      b=['', '']
      try:
        b=bash('cd '+wget_folder+'; wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && mv taxdump.tar.gz '+libraries_folder); assert not b[0]
      except:
        raise notracebackexception, "ERROR fetching taxonomy database from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz : "+b[1]+'\n\nPlease download it manually and provide the names.dmp file to this installation script with -taxdb'
      bbash('cd '+libraries_folder+'; tar -zxf taxdump.tar.gz')
      check_file_presence(libraries_folder+'names.dmp', 'names.dmp taxonomy file')
    else:
      check_file_presence(str(opt['taxdb']), 'taxdb provided')
      if not opt['taxdb'].split('/')[-1]=='names.dmp': raise notracebackexception, "ERROR argument of option -taxdb should be a names.dmp file as downloaded from ncbi taxonomy ftp site."
      print 'Linking '+abspath(opt['taxdb'])+' in '+libraries_folder+' ...'
      bbash('ln -fs '+abspath(opt['taxdb'])+' '+libraries_folder)
  else:
    print "File found: "+libraries_folder+'names.dmp' """

  # gi -> go codes file
  if not opt['no_go']:
    if not is_file(libraries_folder+'idmapping_GI_GO.tab'):
      if not opt['gomap']:
        print "Fetching file with GO annotations for protein gi codes from http://genome.crg.es/~mmariotti/idmapping_GI_GO.tab.gz ..."
        b=['', '']
        try:
          b=bash('cd '+wget_folder+'; wget http://genome.crg.es/~mmariotti/idmapping_GI_GO.tab.gz && mv idmapping_GI_GO.tab.gz '+libraries_folder); assert not b[0]
        except: 
          raise notracebackexception, "ERROR fetching http://genome.crg.es/~mmariotti/idmapping_GI_GO.tab : "+b[1]+'\n\nPlease download it manually and put it in '+libraries_folder+' ; if otherwise you want to skip the installation of the GO tools, run with option -no_go  (or use minimal installation with -min)'
        bbash('gunzip '+libraries_folder+'idmapping_GI_GO.tab.gz')
        check_file_presence(libraries_folder+'idmapping_GI_GO.tab', 'idmapping_GI_GO.tab file')
      else:
        check_file_presence(opt['gomap'], 'GI->GO mapping file provided with -gomap')
        print "Linking "+abspath(opt['gomap'])+' in '+libraries_folder+' ...'
        bbash('ln -s '+abspath(opt['gomap'])+' '+libraries_folder)
    else:
      print "File found: "+libraries_folder+'idmapping_GI_GO.tab'
    if not is_file(libraries_folder+'gene_ontology_ext.obo'):
      if not opt['godb']:
        print "Fetching the complete GO database as obo file from: http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo ..."
        b=['', '']
        try: 
          b=bash('cd '+wget_folder+'; wget http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo && mv gene_ontology_ext.obo '+libraries_folder); assert not b[0]
        except: raise notracebackexception, "ERROR fetching //www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo : "+b[1]+'\n\nPlease download it manually and put it in '+libraries_folder+' ; if otherwise you want to skip the installation of the GO tools, run with option -no_go  (or use minimal installation with -min)'
      else:
        check_file_presence(opt['godb'], 'GO obo file provided with -godb')      
        print "Linking "+abspath(opt['godb'])+' in '+libraries_folder+' ...'
        bbash('ln -s '+abspath(opt['godb'])+' '+libraries_folder)
    else:
      print "File found: "+libraries_folder+'gene_ontology_ext.obo'

  selenoprofiles_programs=["MMlib.py","selenoprofiles_3.py", "selenoprofiles_3.config", "test_selenoprofiles.py",  "blaster_parser.g", "selenoprofiles_build_profile.py", 'selenoprofiles_join_alignments.py', 'selenoprofiles_database.py', 'selenoprofiles_tree_drawer.py']
  words_to_replace={'/users/rg/mmariotti/bin':bin_folder[:-1], '/users/rg/mmariotti/scripts':installation_directory[:-1],  '/users/rg/mmariotti/libraries':libraries_folder[:-1], '/users/rg/mmariotti/temp':temp_directory[:-1], '/users/rg/mmariotti/selenoprofiles/trunk/profiles':installation_directory+'profiles', '/users/rg/mmariotti/Databases/nr.fa':libraries_folder+'nr.fa',  ' -a 7 ':' -a '+str(opt['n_cpus'])+' '}

  ### SECISearch3
  if opt['SS']:
    found_ss=False
    for ss_name in [ 'Seblastian', 'Seblastian.py' ]:
      b=bash('which '+ss_name)
      if not b[0]:  
        found_ss=True          #found
        ex_to_link= b[1]
        print 'SECISearch3/Seblastian found: '+ex_to_link+' ; linking it in '+bin_folder
        bash('ln -s '+ex_to_link+' '+bin_folder+'Seblastian.py')
        break # program found 
    if not found_ss: raise notracebackexception, "ERROR opt -SS active but SECISearch3/Seblastian was not found!"
    else:            print "SECISearch3 found. Activating SECISearch3 action in configuration file ... "
  else: 
    #turning off SECISearch3 action in configuration file'
    words_to_replace['ACTION.post_filtering.secisearch']='#ACTION.post_filtering.secisearch'

  ### bSECISearch3
  if opt['bSS']:
    found_ss=False
    for ss_name in [ 'bSeblastian', 'bSeblastian.py' ]:
      b=bash('which '+ss_name)
      if not b[0]:  
        found_ss=True          #found
        ex_to_link= b[1]
        print 'bSECISearch/bSeblastian found: '+ex_to_link+' ; linking it in '+bin_folder
        bash('ln -s '+ex_to_link+' '+bin_folder+'bSeblastian.py')
        break # program found 
    if not found_ss: raise notracebackexception, "ERROR opt -bSS active but bSeblastian was not found!"
    else:            print "bSeblastian found. Activating bSeblastian action in configuration file ... "
  else: 
    #turning off bSECISearch3 action in configuration file'
    words_to_replace['ACTION.post_filtering.bsecisearch']='#ACTION.post_filtering.bsecisearch'
      
  #reading selenoprofiles programs, making replacements and shebangs
  for program_file in selenoprofiles_programs:
    print "Preparing "+program_file+' ...'
    extension=program_file.split('.')[-1] ; new_file_text=''
    program_file_h=open(installation_directory+program_file, 'r');   its_lines=program_file_h.readlines();   program_file_h.close()
    cont_line=0
    for line in its_lines:
      if cont_line==0 and extension_to_program.has_key(extension) and line.startswith('#!') :     line='#!'+program_to_executable[extension_to_program[extension]]+'\n' #setting she bang for the program file
      else:
        for word in sorted(words_to_replace.keys()):
          replace_to_this=words_to_replace[word]
          line= replace(line, word, replace_to_this)
      new_file_text+= line
      cont_line+=1

    write_to_file(new_file_text, installation_directory+program_file)
    bbash('chmod u+x '+installation_directory+program_file)

  for program_file in [ "blaster_parser.g"]:
    print "Linking "+program_file+' in '+bin_folder+' ...'
    bash('ln -s '+installation_directory+program_file+' '+bin_folder)

  if not updating:
    print 'Linking selenoprofiles...'
    bbash('ln -fs '+os.path.abspath(installation_directory)+'/selenoprofiles_3.py  '+os.path.abspath(installation_directory)+'/Selenoprofiles' )
    print '-- Installation complete. Please link executable '+os.path.abspath('Selenoprofiles')+' to your bin folder '
  else: 
    print 'Linking selenoprofiles...'
    bash('rm '+os.path.abspath(installation_directory)+'/Selenoprofiles')
    bbash('ln -fs '+os.path.abspath(installation_directory)+'/selenoprofiles_3.py  '+os.path.abspath(installation_directory)+'/Selenoprofiles' )
    print '-- Update complete in folder '+installation_directory+' . You can delete this folder ('+bbash('pwd')+') . You can now run '+os.path.abspath(installation_directory)+'/Selenoprofiles' 
    
  print "\nYou may want to run test_selenoprofiles.py now. This programs needs the executable Selenoprofiles to be available in your path.\nIf the test fails, please check http://big.crg.cat/news/20110616/installing_programs_and_modules_needed_by_selenoprofiles or contact us"
  bash('rm -r '+wget_folder)  

except notracebackexception:
  if 'wget_folder' in globals():  bash('rm -r '+wget_folder)
  print sys.exc_info()[1]    
except: raise

#
#  install_selenoprofiles.py
#  
#  Created by Marco Mariotti
