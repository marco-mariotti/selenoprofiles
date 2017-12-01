#!/usr/bin/python -u
__author__  = "Marco Mariotti"
__email__   = "marco.mariotti@crg.eu"
__licence__ = "GPLv3"
__version__ = "3.4"
global temp_folder; global split_folder
from string import *
import sys
import traceback
from types import MethodType
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
try:        import annotations.GO.Parsers.oboparser as oboparser
except:     sys.exc_clear()

allowed_output_formats=['p2g', 'fasta', 'gff', 'gtf', 'cds',  'secis', 'three_prime', 'five_prime', 'dna', 'introns', 'aligned_cds']
allowed_output_formats_descriptions={'fasta':'fasta file with the predicted protein sequence',
'gff' : 'gff file with the genomic coordinates of the prediction',
'gtf':  'analog to gff, but with this syntax for the last field: gene_id "id_for_prediction"; transcript_id "id_for_prediction"; (...)',
'three_prime': 'fasta file with the sequence downstream of the prediction. Its length is defined by the three_prime_length parameter in the main configuration file',
'cds':'fasta file with the coding sequence of the prediction. If frameshifts are predicted, the frameshifts-causing nucleotides are excluded. If the prediction is complete at 3\', the stop codon is included',
'secis':'file with the SECISes found downstream of the prediction (for selenoprotein search - requires SECISearch3)',
'p2g':'selenoprofiles standard output. It contains all essential information for the prediction, including the query-target alignment, the genomic coordinates etc ',
'dna':'fasta file with the full predicted gene sequence, including the intronic sequences (and frameshifts if any)',
'introns':'fasta file with the sequence of the introns, split into different fasta headers',
'five_prime':'fasta file with the sequence upstream of the prediction. Its length can defined by the five_prime_length option or setting it as parameter in the main configuration file',
'aligned_cds': 'fasta file with a pairwise alignment  between the coding sequence and a fake back translation of the blast_query_master sequence  where gaps in the original alignment have been converted to X. useful to obtain alignments of cds between different predictions without realigning'}

help_msg="""Selenoprofiles version """+str(__version__)+""" for profile-based prediction of protein families in nucleotide databases. In a single run, selenoprofiles can search multiple profiles in a single target.
usage:          Selenoprofiles  results_folder  target_file  -s "species name"   -p  ARG   [other options]

As ARG of -p, you must provide one or more profiles (comma separated, if multiple).  Arguments can be filenames (aligned fasta files) or profile names. If they are profile names, files called like profile_name.fa or profile_name.fasta must be present in the profiles folder defined in your main configuration file. Profile alignments are normally built with default options: a profile_name.fa.config file containing the profile settings is created. See script selenoprofiles_build_profile.py to build a profile with non-default options.

The following routines are normally executed only if their output files are not found. You can force their execution specifying either the short or long option. Forcing a routine forces also the execution of all subsequent steps. Note that this may cause certain files to be overwritten, but none will be deleted.

-B || -blast          Run blast and filtering of its output. All next steps are for each blast hit (merged by colinearity)
-E || -exonerate      Run cyclic_exonerate
-G || -genewise       Run genewise. For the blast hits with empty exonerate outputs, genewise is run if the option genewise_to_be_sure is active.
-C || -choose         For each hit, choose the best prediction among available ones (blast, exonerate, genewise)
-F || -filter         Filter results according to p2g_filtering property of profile, then to p2g_refiltering property. 
-D || -database       Store filtered results in a sqlite database
-O || -output         Produce output files according to output options (see below or in the manual).

Before outputing, results are stored in a SQLite database inside the results_folder. If selenoprofiles finds in the database the results from a previous run, it will load them instead of running the pipeline, unless any routine prior to output is specified as above.

Output files will be created inside a folder called like:   results_folder/species_name.target_file_name/output/ 

These output formats are built-in in selenoprofiles: 
"""+join([ (format+':').ljust(2+max([len(i) for i in allowed_output_formats ]))+allowed_output_formats_descriptions.setdefault(format, 'no description available for this method')  for format in allowed_output_formats], '\n')+"""
Use option -output_FORMAT to activate the corresponding output for each prediction. For example option -output_gff will produce gff files.
To see the default active options, see your main configuration file.
Additionally, if a least one prediction is output, a fasta alignment called PROFILE.ali is created: this contains the sequences of the profile along with all predictions for this family in this target. 

A few other options:
-no_splice || -N   for use on RNA sequences or bacterial genomes. Genewise is deactivated in this mode
-genetic_code  +   use a non-standard genetic code; see NCBI codes; implies -tblastn (see -help full)
-test              prints the slave programs and modules available in this selenoprofiles installation, then quits
-print_opt         print currently active options
-h full            display full list of options and of accessory programs

Please refer to the manual, available at http://big.crg.cat/services/selenoprofiles
If selenoprofiles has been useful for your research, please cite:  
Mariotti M, Guigo R. Selenoprofiles: profile-based scanning of eukaryotic genome sequences for selenoprotein genes.
Bioinformatics. 2010 Nov 1;26(21):2656-63."""
full_help="""
####### Full list of other options #######
Options with + require an argument. Options can be specified either in command line with syntax -option value or in the configuration file as option = value. If no argument is required for an option, in the configuration file you must use value=1.  To deactivate an option which is active by default, use -option 0.

## system and global configuration
-bin_folder       +     folder where the executables run by selenoprofiles are searched
-profiles_folder  +     folder where the profile alignments are searched
-temp             +     temporary folder. A folder with random name is created here, used and deleted at the end of the computation
-save_chromosomes       active by default. Selenoprofiles tries to recycle the single-sequence fasta files extracted from the genome, to minimize computation. These files are kept in a subfolder of the folder specified with -temp. Turn off to save disk space
-no_colors              disable printing in colors to atty terminals. Put "no_colors=1" in your configuration file to set this as default
-GO_obo_file      +     path to the gene_ontology_ext.obo file used in GO tools-based filtering (see manual)

## prediction programs
-dont_exonerate         do not run exonerate. Not recommended. 
-dont_genewise          do not run genewise.  Use to reduce the time required for computation
-genewise_to_be_sure    active by default. When exonerate produce no output or its prediction does not overlap the seed blast hit, genewise is run, seeded using the blast hits coordinates. Turn off this option not to run genewise in these cases, to reduce the time required for computation
-no_blast               do not allow choosing a blast prediction (over a genewise or exonerate prediction). Use this if an accurate splice site prediction is crucial for you
-tblastn                use simple tblastn (single query) instead of psitblastn (profile based PSSM)
-exonerate_extension      +    nt lenght of extension used on both sides by the cyclic exonerate procedure (see paper or manual)
-genewise_extension       +    nt length of extension used on both sides when running genewise on gene boundaries defined by exonerate 
-genewise_tbs_extension   +    nt length of extension used on both sides when running genewise blast hits for which exonerate produced no output (only if option -genewise_to_be_sure is active)
-blast_opt/exonerate_opt/genewise_opt  "+"    command line options always used when running psitblastn/exonerate/genewise. Additional profile-specific options are also available (see profile attributes below)

## database 
-max_attempts_database  +   if you're running multiple parallel instances of selenoprofiles on the same target, sometimes selenoprofiles will find the database locked. This sets the maximum number of times that selenoprofiles will try again before crashing.
-sleep_time             +   when the database is found locked, this is the time in seconds separating two consecutive attempts
-full_db                    normally, selenoprofiles will store in the database only the results passing all filters. If this option is active, all results are written, allowing for fast retrieving (see also option -state)
-no_db                      do not store anything in the database, and exit after the filtering step. This option make sense only if you want to heavily parallelize the use of selenoprofiles, when thousands of accesses to the database may become a problem; this option allow to produce all intermediate text files, but will produce no output. Selenoprofiles must be run again without this option when files for all profiles are ready, so that results will be stored in the database, redundancy of results will be removed and output will be produced
-stop                       exit after storing results to the database, producing no output files. When all results are ready in the database, you should run selenoprofiles again with option -merge to remove redundancy and get output.
-merge                      force running the procedure to remove inter-profile redundancy of all results in the database, and proceed to output
-no_merge                   force skipping the procedure to remove inter-profile redundancy

## additional output options (see above for basics)
-five_prime_length  +   length in nucleotide used when outputing the five prime sequence (output option -output_five_prime must be active)
-three_prime_length +   length in nucleotide used when outputing the three prime sequence (output option -output_three_prime must be active)
-outfolder          +   use this as output folder, instead of results_folder/species_name.target_file_name/output/
-output_ali             active by default. A fasta alignment file is produced for each family comprising all filtered results. Turn off to speed up outputing
-output_FORMAT_file +   redirect all output of the specified format (see above for available formats) to the argument file
-output_filter     "+"  argument is a python expression evaluating a p2g result, named x. Only the results for which the evaluation is true are output. e.g:  -output_filter "x.chromosome=='chr1'" Note: if option -full_db is active, the filtered out results will also be considered unless specifically excluded by the output filter
-state              +   determines what kind of filtering label the results must have to be output. Possible values: kept, filtered, refiltered, redundant, overlapping. Multiple values can be provided comma-separated. Note: this option will work only if option -full_db was active when writing results to the database
-fasta_add         "+"  argument is a python expression manipulating a p2g result, named x. The expression is evaluated to a string which is added to the fasta header used for the result (for example in the fasta output, cds output, ali output etc). Useful to quickly add custom information to the output 

## miscellaneous
-fam_list         +       file with list of profiles, one per line. This is an alternative way to provide the list of profiles that overrides option -p
-name             +       specify manually the target file name. Normally it is derived from the file name
-log              +        write a copy of the output normally printed to screen into a file provided as argument
-add              +       provide a file with python code that is executed before the pipeline flow. It can be used to customize selenoprofiles (see manual)
-clean                    remove intermediate files (blast, exonerate, genewise etc), to save disk space. Use only for non-parallelized run of selenoprofiles, or making sure that results are stored in the db before cleaning
-filtered_blast_file  +   write an output file with all blast hits passing the blast filter. Mostly for refining filters
-blast_filtering_warning  active by default. When too many blast hits pass filtering for a profile, normally the program prints a warning and goes on with the first blast hits encountered. If you turn this option off, in these cases the program will crash instead
-debug                    if an error occurs, instead of crashing directly, it prints a preview of the error, then waits for keyboard input before exiting. Useful when debugging to inspect the files in the temporary folder (deleted when exiting)
-print_commands           print to screen every bash command before running it. Extremely verbose

#### Families set ####
If you often run the same list of families, you may want to define them with a keyword. This is done by adding a line in the configuration file with the syntax:
families_set.NAME  =  fam1,fam2,fam3
You can then use such defined keywords (NAME in the example) as argument(s) of the -p option. Families set definition can include other sets names. Examples:  families_set.first_set  =  fam1,fam2,fam3        families_set.extended_set  =  first_set,fam4,fam5

#### Actions ####
Actions are a fast way to customize the pipeline workflow for specific purposes. Each action is performed on all predictions active at a specific moment of the pipeline, depending on the action category. The possible categories of actions are post_blast_filter, post_blast, post_blast_merge, pre_choose, pre_filtering, post_filtering, pre_output.
Actions can be specified either in the command line or in the main configuration file, like other options. 
The syntax for command line actions  is:  -ACTION.category.name "+", where name is a user-defined name for the action and the argument is a piece of python code that will be executed, in which you can refer to the active p2g prediction as x. 
Actions can be used to modify/improve predictions (see for example the default actions in the main configuration file) or to add custom output to the log of selenoprofiles.
Example of command line action: (prints the sequence of results labelled as pseudogenes)
 -ACTION.post_filtering.print_pseudo       " if x.label=='pseudo': print x.output_id()+' '+ x.protein() "
See selenoprofiles manual for more information.

####### Profile attributes #######
These attributes are read from the .config file of each profile alignment, and can be manually added to such .config file, or chosen when running selenoprofiles_build_profile.py.
In the .config file, attributes are set with the syntax   attribute_name = value
All attributes not specified  for a certain profile are assumed to be defaults read from the main configuration file, there in the syntax attribute_name.DEFAULT = value 
If you want to use the same value for many profiles, then you should set attributes with keywords. In this case you need to define an attribute keyword name in the main configuration file with the syntax attribute_name.KEYWORD_NAME = value ; then, you can use the keyword in the .config file of any profile using the syntax attribute_name = KEYWORD_NAME ; as an example, see the SELENO keyword for attributes blast_options, exonerate_options and genewise_options in the main configuration file. This particular keyword is used for selenoprotein families.

### Full list of attributes (all attributes require an argument)
.name                  profile name, normally derived from the profile alignment filename [compulsory]
.blast_filtering       python expression evaluated to decide if a blast hit (named x) passes blast filtering; e.g.  x.evalue < 1e-2
.p2g_filtering         python expression evaluated to decide if a p2g hit (named x) passes filtering. This filter is applied after the prediction choice step; e.g.  x.coverage() >= 0.5
.p2g_refiltering       python expression evaluated to decide if a p2g hit (named x) passes refiltering. This filter is applied after the p2g_refiltering to allow two types of discarded hits, filtered and refiltered.
.max_blast_hits        maximum number of blast hits allowed for this profile. See also option blast_filtering_warning.
.queries               defines which sequences in the profile are elegible as queries for exonerate and genewise. There are several possible syntaxes: 
    all                 all sequences are queries (default)
    [...]               you provide directly the list of queries titles, or indices (starting from 0) in the ordered titles list
    last_query:TITLE    you provide the last title in the ordered titles list that will be labelled as query. This may be a incomplete title: it just have to be found in a title
    best:0.XX           you provide the proportion of the full list of queries in the ordered list that will be labelled
    NOTE for selenoprotein families: in the first, third and fourth cases, those titles for which at least one U position in the alignment contains a gap in the sequence are not considered as queries, unless the you prefix the value of the queries attribute with "not_only_U:"
  Examples of valid values for attribute queries:   best:0.5    last_query:XP_854030   not_only_U:all   not_only_U:best:0.8  
.clustering_seqid                    before running psitblastn, the profile alignment is split into clusters, and a search is performed for each of them, to ensure adequate sensitivity. This is the maximum sequence identity within each cluster
.max_column_gaps     within each cluster, a blast query is computed by building a consensus of sequences in the cluster. This sets the maximum proportion of gaps in a column of a blast cluster, above which the position will not be present in the blast query
.tags                  list of tags (strings in perl-like regexp syntax) used by the tag_score method, for tag blast filtering (see manual). e.g. ['deiodinase', '\\WDI\\d ']
.neutral_tags          list of tags scored neutrally in the tag_score method. This allows skipping non-informative titles without penalty
.tag_db                path to the fasta file used as database for the tag blast procedures, used when tag_score and go_score methods are invoked. Normally the same nr database can be used for all profiles but you may differentiate them to speed up or improve the process
.tag_blast_options     command line options used when running tag blast (blastp agains tag_db) for this profile
.go_terms              list of GO terms (as "GO:XXXXXXX" strings) used by the go_score method, for gene ontology/tag blast based filtering (see manual). e.g. ['GO:0004364']
.uniref2go_db          path to the file mapping the GO terms to the proteins in the tag_db. Used for go_score 
.blast_options/exonerate_options/genewise_options         command line options used when running psitblastn/exonerate/genewise for this profile

####### Accessory programs #######
This pipeline comes with a few accessory programs:
selenoprofiles_build_profile.py       to build a custom profile alignment, or inspect its features for filtering, such as AWSI distribution or associated GO terms
selenoprofiles_database.py            to manipulate a selenoprofiles sqlite database, to clean it, or for fast output of stored results 
selenoprofiles_join_alignments.py     useful to collect an alignment for profiles searched in multiple targets, typically different species
selenoprofiles_tree_drawer.py         for graphical output of the predicted genes along a species tree
"""

terminal_colors={'routine':'blue', 'database':'magenta', 'profile':'green'}
set_MMlib_var('colored_keywords', {'ERROR':'red,underscore', 'WARNING':'red'})

def load(config_filename='/users/rg/mmariotti/scripts/selenoprofiles_3.config', args={}):
  """Load all global variables later used, reading configuration file and command line options. """  
  ### initialising command line options and initialising file / objects lists . opt is a dictionary the all the command_line options. (0: False 1:True in case of boolean values). The -config FILE option allows to specify a different configuration file.  Every option in the configuration file is read, then eventually replaced by the correspondant command line option.
  for i in range(len(sys.argv)):
    if sys.argv[i] == "-config":     config_filename=sys.argv[i+1]
  command_line_synonyms={'B':'blast', 'E':'exonerate','G':'genewise','C':'choose', 'F':'filter',  'O':'output', 'P':'profile', 'S':'species', 's':'species', 'p':'profile', 'D':'database', 'no_splicing':'no_splice', 'N':'no_splice'}
  non_config_options=['t', 'v', 'name', 'blast', 'exonerate', 'genewise', 'filter', 'choose', 'output', 'fam_list', 'print_commands', 'species', 'dont_exonerate', 'dont_genewise', 'debug', 'log', 'state', 'filtered_blast_file', 'output_filter', 'add', 'clean', 'fasta_add', 'outfolder', 'five_prime_length', 'database', 'no_blast', 'stop', 'no_merge', 'no_db', 'merge', 'test', 'no_splice', 'genetic_code', 'tblastn']
  for keyword in allowed_output_formats: non_config_options.extend(['output_'+keyword+'_file', 'output_'+keyword])
  # reading configuration file
  def_opt= configuration_file(config_filename)
  # complete list of global variables
  global families_sets, keywords, opt, sleep_time, max_attempts_database, temp_folder, split_folder, bin_folder, profiles_folder, target_file, reference_genome_filename, three_prime_length, five_prime_length, profiles_names, profiles_hash, results_folder, target_name,  target_species, target_results_folder, results_db_file, results_db, actions, blast_folder, exonerate_folder,  genewise_folder, prediction_choice_folder, filtered_list_folder,  output_folder, blast_nr_folder, target_file_index, chromosome_length_file, exonerate_extension, genewise_extension, genewise_tbs_extension, blast_options, exonerate_options, genewise_options, genewise_tbs_options, output_file_handlers, max_chars_per_column
  # setting sets of families, defined by keywords starting with families_set. in the config file. e.g families_set.eukaryotic = sps,GPx,MsrA,DI,15-kDa,Fep15
  families_sets={}      # example: key_   machinery  ->   value_   ['sps', 'eEFsec', 'pstk', ... ] 
  if def_opt.has_key('families_set'):
    for k in def_opt['families_set']:
      families_sets[ k ] = [del_white(word) for word in def_opt['families_set'][k].split(',')]
  # setting filtering defaults or keywords functions. here, we define dynamically functions in this environment, reading their text representation (python code) from the configuration file
  
  for x in [x for x in def_opt   if '.' in x]:    del def_opt[x]  
  
  keywords={}; keywords_text={}  
  for category in profile_alignment.parameters:
    keywords[category]={}
    if not def_opt.has_key(category): raise notracebackException, "ERROR configuration file lacks default value for "+category+'  (example:  '+category+'.DEFAULT = value ) '
    keywords_text[category]=def_opt[category]
    #write((category, def_opt[category]), 1, how='green')

    for kword in def_opt[category]: # one is surely DEFAULT, then more can be defined. 
      function_text=def_opt[category][kword]
      if 'filtering' in category:       function_text='lambda x:'+function_text
      elif 'options' in category or '_db' in category:       function_text='"""'+function_text+'"""'
      try:        exec('keywords[category][kword]='+str(function_text))
      except:
        printerr("selenoprofiles ERROR can't assign function "+str(function_text)+' to keyword '+str(kword)+' in category '+str(category), 1)
        raise       
  

  #preparing to read command line options     
  for i in non_config_options: 
    if not i in def_opt:  def_opt[i]=0
  ## the "if args" construct is to allow running this "load" function from inside an external python script. Normally, options are read from command line
  if args:   
    opt=options()
    for k in def_opt: 
      if args.has_key(k):  opt[k]=args[k]
      else: opt[k]=def_opt[k]
  else:      opt=command_line(def_opt, help_msg, ['r','t'], synonyms=command_line_synonyms, tolerated_regexp=['ACTION.*']+[p+'.*' for p in profile_alignment.parameters], strict= notracebackException, advanced={'full':full_help} );   #### reading options from command line
  
  # allowing to specify profile parameters on the command line
  for x in opt:
    if '.' in x:
      category=x[:x.find('.')]
      kword=   x[x.find('.')+1:]
      value=opt[x]
      #write( (category, kword, value), 1, how='red')
      keywords_text[category][kword]=str(value)
      #write((category, kword, value), 1, how='magenta')
      #write(keywords_text[category], 1, how='green')
      function_text=value
      if 'filtering' in category:       function_text='lambda x:'+function_text
      elif 'options' in category or '_db' in category:       function_text='"""'+function_text+'"""'
      try:        exec('keywords[category][kword]='+str(function_text))
      except:
        printerr("selenoprofiles ERROR can't assign function "+str(function_text)+' to keyword '+str(kword)+' in category '+str(category), 1)
        raise            

  set_MMlib_var('opt', opt);
  sleep_time=opt['sleep_time']
  max_attempts_database=opt['max_attempts_database']
  
  if opt['log']:  
    if opt['log']==1: raise notracebackException, "selenoprofiles ERROR an log output must be provided with option -log"
    global log_file; log_file=open(opt['log'], 'w'); set_MMlib_var('log_file', log_file)
  # set options and parameters 
  temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'Please provide a different folder with option -temp'); set_MMlib_var('temp_folder', temp_folder);
  if opt['save_chromosomes']:  split_folder=Folder(opt['temp'])
  else:                        split_folder=temp_folder
  test_writeable_folder(split_folder, 'split_folder'); set_MMlib_var('split_folder', split_folder)
  bin_folder=Folder(opt['bin_folder']);     set_MMlib_var('bin_folder', bin_folder)
  profiles_folder=Folder(opt['profiles_folder']);      check_directory_presence(profiles_folder, 'profiles_folder', notracebackException)

  if opt['genetic_code']:
    set_genetic_code(opt['genetic_code'])
    opt['tblastn']=1
    if not opt['dont_exonerate']:
      exonerate_version=float(bash('exonerate --version')[1].split('\n')[0].split()[-1][:3])
      if exonerate_version<2.4: 
        raise notracebackException, "ERROR exonerate version detected: {}\nAlternative genetic codes are bugged in exonerate versions <2.4!\nPlease download and install exonerate version 2.4.0 if you want to use -genetic_code".format(exonerate_version)

    for  k  in  keywords['blast_options']:        keywords['blast_options'][k]+=    ' -D '+str(opt['genetic_code'])
    for  k  in  keywords['exonerate_options']:    keywords['exonerate_options'][k]+=' --geneticcode '+str(opt['genetic_code'])
    for  k  in  keywords['genewise_options']:     keywords['genewise_options'][k]=keywords['genewise_options'][k].format(GENETIC_CODE=opt['genetic_code'])
  else: 
    for  k  in  keywords['genewise_options']:     keywords['genewise_options'][k]=keywords['genewise_options'][k].format(GENETIC_CODE=1)

  ###############
  ### test routine 
  if opt['test']:
    write('\n'); write('TEST', how=terminal_colors['routine']); write(' routine', 1)
    write('Slave programs:', 1)
    write(' blastall    (ncbi blast)     : ')
    b=bash('blastall ')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! blastall is compulsory in the selenoprofiles pipeline, please install it!', 1)
    write(' blastpgp    (ncbi blast)     : ')
    b=bash('blastpgp --help')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! blastpgp is compulsory in the selenoprofiles pipeline, please install it!', 1)
    write(' formatdb    (ncbi blast)     : ')
    b=bash('formatdb --help')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! formatdb is compulsory in the selenoprofiles pipeline, please install it!', 1)
    write(' exonerate   (exonerate pkg)  : ')
    b=bash('exonerate ')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! you could still run using option -dont_exonerate' , 1)
    write(' fastaindex  (exonerate pkg)  : ')
    b=bash('fastaindex ')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! the fasta suite by exonerate is compulsory for the selenoprofiles pipeline, please install it! ' , 1)
    write(' fastafetch  (exonerate pkg)  : ')
    b=bash('fastafetch ')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! the fasta suite by exonerate is compulsory for the selenoprofiles pipeline, please install it! ' , 1)
    write(' fastasubseq (exonerate pkg)  : ')
    b=bash('fastasubseq ')
    if b[0]==256:    write('found', 1)
    else:            write('NOT FOUND! the fasta suite by exonerate is compulsory for the selenoprofiles pipeline, please install it! ' , 1)
    write(' genewise    (Wise2 package)  : ')
    b=bash('genewise ')
    if b[0]==16128:    write('found', 1)
    else:            write('NOT FOUND! you could still run using option -dont_genewise' , 1)
    write(' SECISearch3 (Seblastian)     : ')
    b=bash('Seblastian.py ')
    if b[0]==2048:    write('found', 1)
    else:            write('NOT FOUND! SECISearch is recommended to search for selenoproteins; you can ignore this if you want to search your custom profiles' , 1)
    ###     ##

    write('\nPython modules and miscellaneous:', 1)    
    write(' oboparser   (module for GO)  : ')
    try:        
      import annotations.GO.Parsers.oboparser as oboparser
      write('found', 1)
    except:     
      write('NOT FOUND! you cannot use method go_score() for filtering. This precludes the use of the built-in selenoprotein profiles', 1)
      sys.exc_clear()
    write(' obofile     (GO structure)   : ')
    if is_file(opt['GO_obo_file']):                 write('found', 1)
    else:                                           write('NOT FOUND! you cannot use method go_score() for filtering. This precludes the use of the built-in selenoprotein profiles', 1)
    write(' uniref2go db    (GO annotation)  : ')
    if  is_file( keywords ['uniref2go_db']['DEFAULT'] ):  write('found', 1)
    else:                                           write('NOT FOUND! You cannot use methods go_score(), unless you defined a different uniref2go_db attribute for that profile. This precludes the use of the built-in selenoprotein profiles', 1)
    write(' uniref database (tag/GO blast)   : ')
    if  is_file( keywords ['tag_db']['DEFAULT'] ):  write('found', 1)
    else:                                           write('NOT FOUND! You cannot use methods tag_score() and go_score(), unless you defined a different tag_db attribute for that profile. This precludes the use of the built-in selenoprotein profiles', 1)
    write(' Pylab  (interactive graphs)  : ')
    try:
      import Pylab
      write('found', 1)
    except:
      write('NOT FOUND! you cannot use the graphical tools ( -d or -D ) of selenoprofiles_build_profile.py', 1)
      sys.exc_clear()
    write(' ete2   (interactive trees)   : ')
    try:
      import ete2
      write('found', 1)
    except:
      write('NOT FOUND! you cannot use the results viewer program called selenoprofiles_tree_drawer.py', 1)
      sys.exc_clear()

    sys.exit()  
  ###############
  ##### setting target
  target_file=   opt['t'];                               check_file_presence(target_file, 'target file', notracebackException)
  target_file=abspath(target_file)
  reference_genome_filename= target_file; set_MMlib_var('reference_genome_filename', target_file) #to allow some [gene].fasta_sequence() calls without specifying the genome file. see [blasthit].place_selenocysteine() using [blasthit].cds() 
  three_prime_length=opt['three_prime_length'];   set_MMlib_var('three_prime_length', three_prime_length) 
  five_prime_length=opt['five_prime_length'];     set_MMlib_var('five_prime_length', five_prime_length)
  if opt['output_three_prime'] and not three_prime_length: raise notracebackException, "ERROR output_three_prime is active but three_prime_length is 0 or not defined!"
  if opt['output_five_prime'] and not five_prime_length:   raise notracebackException, "ERROR output_five_prime is active but five_prime_length is 0 or not defined!"
  # setting routine cascade
  next_must_be_active=0
  for o in ['blast', 'exonerate', 'genewise', 'choose', 'filter','database',  'output']:
    if opt[o]: next_must_be_active=1
    elif next_must_be_active: opt[o]=1
  if opt['print_commands']:   set_MMlib_var('print_commands', 1)
  ## determining families
  if opt['fam_list']:
    write('Reading families list from file: '+str(opt['fam_list']), 1)
    check_file_presence(opt['fam_list'], 'families list file', notracebackException)
    opt['profile']= join( [line.strip() for line in open(opt['fam_list'], 'r')], ',')
  # here I change each set of families with the list of its members, iteratively
  profile_files=opt['profile'].split(',')
  profiles_names=[]; profiles_hash={}
  index_profile=0
  while any ([families_sets.has_key(p) for p in profile_files]):
    this_profile=profile_files[index_profile]
    while families_sets.has_key(this_profile):
      profile_files[index_profile:index_profile+1]= families_sets[this_profile]
      this_profile=profile_files[index_profile]
    index_profile+=1
  # removing redundancy in profiles names
  profiles_files_length=len(profile_files)
  for i in range(profiles_files_length-1, -1, -1):
    profile_file=profile_files[i]
    if profile_file in profile_files[:i]:
      profile_files.pop(i)  
  # changing profile names with profile file
  write('/'+'-'*119, 1)
  for f in range(len(profile_files)):
    this_profile=profile_files[f]
    if this_profile == 'NONE':       continue
    this_profile_found=False
    for ext in ['', '.fa', '.fasta', '.ali']:
      for profiles_folder_attempt in ['', profiles_folder]:
        file_attempt=profiles_folder_attempt+this_profile+ext
        if not this_profile_found and is_file(file_attempt):        
          profile_files[f]=  file_attempt
          this_profile_found=True
    if not this_profile_found:                  raise notracebackException, "selenoprofiles ERROR can't find profile: "+this_profile
    else:
      write( '|-'+'-'*abs(f%8-4)+'  load profile  '+profile_files[f], 1)
      profile_ali=profile_alignment(profile_files[f])
      profiles_names.append(profile_ali.name)
      if profiles_hash.has_key(profile_ali.name):   raise notracebackException, "selenoprofiles ERROR profile names must be unique! This was found twice: "+profile_ali.name
      profiles_hash[profile_ali.name] = profile_ali
  write('\\'+'-'*119, 1)
  
  # filtering state
  if not opt['state']: opt['state']='kept'
  ## determining target name and output folders
  results_folder=Folder(opt['r']);                     
  if not opt['outfolder']:         test_writeable_folder(results_folder, 'selenoprofiles results folder')
  if opt['name']:   target_name=opt['name']
  else: 
    target_name=base_filename(target_file)
    if target_name.split('.')[-1] in ['fasta', 'fa']:      target_name=join(target_name.split('.')[:-1], '.')
  target_name=replace_chars(target_name, '.', '_')
  #species_library=opt['species_library']; check_file_presence(species_library, 'species_library'); set_MMlib_var('species_library', species_library)  ## this will be used for getting species name
  ## database for output,   linking to species 
  target_species=None; set_species_from_previous_run=False   
  if not target_name=='genome' and not opt['species']:
    ### no species information provided!
    # checking if there's a results subfolder with results on the same target, so we can use the same species    
    service('no species provided. searching existing folders for a previous run on the same target... ')
    target_links_found=bbash('find '+results_folder+' -maxdepth 2 -mindepth 2 -name link_target.fa -type l ').split('\n') 
    if target_links_found==['']: target_links_found=[]
    for possible_target_link in target_links_found:
      try: 
        possible_target_file=       dereference(possible_target_link)
        if possible_target_file == target_file: 
          masked_species_name  =  possible_target_link.split('/')[-2].split('.')[0]
          #taxid, library_species_name= get_species_from_library(    unmask_characters( replace(masked_species_name, '_', ' ') )     )
          #target_species=species( library_species_name );  target_species.taxid = taxid
          target_species=unmask_species( masked_species_name )
          write('Previous run on the same target found. Species name: '+str(target_species), 1)
          set_species_from_previous_run=True
          break
      except: raise
    if not set_species_from_previous_run:
    ## if we didn't find it, we set species to unidentified, printing a warning
      target_species='unidentified'
      printerr('WARNING no species name provided! Using: "unidentified"', 1)
  else:
    ### something provided; either the target is called genome, so we can derive the species name from the folder, or the option -species (-s) was provided. setting variable species_name, which may be converted to ncbi syntax
    if target_name=='genome' and not opt['species']:       target_species=abspath(target_file).split('/')[-2] 
    elif opt['species']:                                   target_species=opt['species']
    else:                                                  target_species='unidentified'    
    if '_' in target_species or '{ch' in target_species:  target_species=unmask_species(target_species)
#    try:   
#      taxid, library_species_name= get_species_from_library( species_name )
#      if not library_species_name in  [species_name, unmask_characters(species_name), unmask_characters(replace_chars(species_name, '_', ' '))]:  
#        printerr('WARNING species name '+species_name+' changed to ncbi-defined scientific name: '+library_species_name, 1)
#      target_species=species( library_species_name ) 
#      target_species.taxid=taxid 
#    except Exception:   
#      raise notracebackException, 'ERROR species not recognized! Only ncbi taxonomy names are accepted'

  ## if species name is unidentified, we make a uniq target name to be sure not to overwrite results
  if target_species=='unidentified':  target_name=uniq_id_for_file(target_file)
  target_results_folder=Folder(results_folder + mask_species(target_species) + '.'+target_name) #if you change this, change above as well, in the end of mv command

  results_db_file=target_results_folder+'results.sqlite'
  if not opt['no_db']:
    results_db=selenoprofiles_db(results_db_file)
    if not results_db.has_table('results'): 
      results_db.initialise_db() #creating database if necessary
      write('DATABASE', how=terminal_colors['database']); write(' initialising: '+results_db_file, 1)
    results_db.update_to_last_version() # attempt updating database from last versions of selenoprofiles if necessary

  link_target_file=target_results_folder+'link_target.fa'
  if not is_file(link_target_file):
    write('Creating link:  '+link_target_file+'   --> '+target_file, 1)
    bash('rm '+link_target_file) # if we entered here because link is broken, cleaning up and relinking
    bbash('ln -fs '+target_file+' '+link_target_file)
    
  # preparing actions hash. keys of actions are categories (pre_filtering, post_filtering), and their values are hashes. Keys of these hashes are ids (can be anything) which are evaluated in alphabetical order (the sorted() function is used). Values corresponding to such keys are strings which are executed in the python environment, denoting the chosen_hit object with "x" 
  actions={}       ; action_categories=['pre_blast_filter', 'post_blast', 'post_blast_merge', 'pre_choose', 'pre_filtering', 'post_filtering', 'pre_output']
  
  if not def_opt.has_key('ACTION'): printerr('WARNING no action is specified in the config file', 1)  
  else: actions=def_opt['ACTION']
  for category in action_categories:
    if not actions.has_key(category): actions[category]={}
  #expanding with action defined in command line
  for k in opt:
    if k.startswith('ACTION.'): 
      try:
        category=   k.split('.')[1]
        if not category in action_categories:    raise Exception
        action_id=  k.split('.')[2] 
      except: raise Exception, "selenoprofiles ERROR actions must have the form -ACTION.category.id \"action code\" (command line) or ACTION.category.id = \"action code\" (configuration file),   where category must be one of "+join(action_categories, ', ')+" and id is any non-null string."
      actions[ category ] [ action_id ]= opt[k]

  # setting folders depending on keep_ options
  blast_folder=Folder(target_results_folder+'blast')
  exonerate_folder=Folder(target_results_folder+'exonerate')
  genewise_folder=Folder(target_results_folder+'genewise')
  prediction_choice_folder=Folder(target_results_folder+'prediction_choice')
  filtered_list_folder=Folder(target_results_folder+'filtering')
  output_folder=Folder(target_results_folder+'output')
  if opt['outfolder']:     output_folder=Folder(opt['outfolder']); test_writeable_folder(output_folder, 'outfolder');
  blast_nr_folder=Folder(target_results_folder+'tag_blast')

  ## indexing and formatting if necessary
  # fastaindex
  target_file_index = join( target_file.split('.')[:-1],'.'   )+'.index';   index_in_progress_file=join( target_file.split('.')[:-1],'.'   )+'.index_in_progress'
  if not is_file(target_file_index) and not is_file(index_in_progress_file):     
    write_to_file('look in '+temp_folder, index_in_progress_file) ## creating file to indicate we're computing the index
    write('Indexing '+target_file+' with fastaindex ... ', 1)
    try:      
      target_file_index = fastaindex(target_file, silent=True)
      bbash('rm '+index_in_progress_file)
    except: 
      bash('rm '+index_in_progress_file)
      raise
  elif is_file(index_in_progress_file):
    index_waited=0
    while is_file(index_in_progress_file):
      write('Another instance of selenoprofiles is computing the index of '+target_file+' ; waiting ...', 1)
      index_waited+=1
      time.sleep(sleep_time)
      service('waited: '+str(sleep_time*index_waited)+' seconds ')
    if not is_file(target_file_index): raise notracebackException, "ERROR after waiting, the index was not found : "+target_file_index
  # formatdb
  formatdb_in_progress_file=  join( target_file.split('.')[:-1],'.'   )+'.formatdb_in_progress'
  if bash('ls '+target_file+'.*nin')[0]:
    write_to_file('look in '+temp_folder, formatdb_in_progress_file) ## creating file to indicate we're computing with formatdb
    write('Formatting '+target_file+' with formatdb ... ', 1)  
    try:      
      formatdb(target_file, is_protein=False, silent=True) 
      bbash('rm '+formatdb_in_progress_file)
    except: 
      bash('rm '+formatdb_in_progress_file)
      raise
  elif is_file(formatdb_in_progress_file):
    index_waited=0
    while is_file(formatdb_in_progress_file):
      write('Another instance of selenoprofiles is formatting '+target_file+' ; waiting ...', 1)
      index_waited+=1
      time.sleep(sleep_time)
      service('waited: '+str(sleep_time*index_waited)+' seconds ')
    if bash('ls '+target_file+'.*nin')[0]: raise notracebackException, "ERROR after waiting, the file produced by formatdb were not found : "+target_file+'.*nin'
  #chromosome lenghts      
  chromosome_length_file=join( target_file.split('.')[:-1],'.'   )+'.lengths';   chromosome_length_in_progress_file= join( target_file.split('.')[:-1],'.'   )+'.lengths_in_progress'
  if not is_file(chromosome_length_file) and not is_file(chromosome_length_in_progress_file):     
    write_to_file('look in '+temp_folder, chromosome_length_in_progress_file) ## creating file to indicate we're computing the index
    write('Computing length of all sequences in '+target_file+' --> '+chromosome_length_file, 1)    
    try:      
      bbash('fastalength '+target_file+ ' > '+temp_folder+'chrom_lengths && mv '+temp_folder+'chrom_lengths '+chromosome_length_file)
      bbash('rm '+chromosome_length_in_progress_file)
    except: 
      bbash('rm '+chromosome_length_in_progress_file)
      raise
  elif is_file(chromosome_length_in_progress_file):
    index_waited=0
    while is_file(chromosome_length_in_progress_file):
      write('Another instance of selenoprofiles is computing the length of the sequences in '+target_file+' ; waiting ...', 1)
      index_waited+=1
      time.sleep(sleep_time)
      service('waited: '+str(sleep_time*index_waited)+' seconds ')
    if not is_file(chromosome_length_file): raise notracebackException, "ERROR after waiting, the sequence lengths file was not found : "+chromosome_length_file
  
  # handling no_splice option
  if opt['no_splice']:     opt['dont_genewise']=1

  # determining parameter options of programs: blast, exonerate, genewise. These are the default options, but they can be overriden by the profile options.
  exonerate_extension=opt['exonerate_extension']
  genewise_extension=opt['genewise_extension']
  genewise_tbs_extension=opt['genewise_tbs_extension']
  blast_options={};  blast_options_string_split=str(opt['blast_opt']).split()
  while blast_options_string_split:
    try:
      option_name = blast_options_string_split.pop(0)[1:];        value = blast_options_string_split.pop(0) #######
      if value[0] in '\'"':  
        done=0
        while len(value)==1 or value[-1] != value[0] : 
          value+=blast_options_string_split.pop(0)+' '
      blast_options[option_name]=value
    except: raise notracebackException, 'selenoprofiles ERROR parsing blast options: '+blast_options_string_split
  exonerate_options={};  exonerate_options_string_split=str(opt['exonerate_opt']).split()
  while exonerate_options_string_split:
    if len(exonerate_options_string_split)==1:    raise notracebackException, "selenoprofiles exonerate_opt ERROR the words number is not even... can't find a value for option: "+str(exonerate_options_string_split[0])
    option_name = exonerate_options_string_split.pop(0)[1:];        value = exonerate_options_string_split.pop(0)
    exonerate_options[option_name]=value
  genewise_options={};  genewise_options_string_split=str(opt['genewise_opt']).split()
  while genewise_options_string_split:
    if len(genewise_options_string_split)==1:    raise notracebackException, "selenoprofiles genewise_opt ERROR the words number is not even... can't find a value for option: "+str(genewise_options_string_split[0])
    option_name = genewise_options_string_split.pop(0)[1:];        value = genewise_options_string_split.pop(0)
    genewise_options[option_name]=value
  genewise_tbs_options=genewise_options.copy()

  # preparing filehandlers for output files like fasta, gff, gtf.. if any has been specified in command_line ( e.g. -output_fasta_file any_file.fa )
  output_file_handlers={} # will host the filehandlers to write in specified output files, if any 
  for keyword in allowed_output_formats:
    if opt['output_'+keyword+'_file']:
      output_file_handlers[keyword]= open( opt['output_'+keyword+'_file'], 'w')
      try:     output_file_handlers[keyword]= open( opt['output_'+keyword+'_file'], 'w')
      except:  
        printerr('ERROR can\t open option -'+'output_'+keyword+'_file'+' file for writing: '+opt['output_'+keyword+'_file'])
        raise

  # for pretty printing
  max_chars_per_column={'id':18+4, 'chromosome':10, 'strand':1}

  ## computing summary of current options  
  summary=''
  summary+='       Options      '.center(120, '#')+'\n'
  #summary+='command line options: '
  output_options= []
  for k in sorted(opt.keys()):
    #if def_opt.has_key(k) and opt[k]!=def_opt[k]: summary+='-'+str(k)+' '+str(opt[k])+' '
    if k.startswith('output_') and opt[k]: output_options.append(k)  

  #summary+='\n'
  summary+='| output folder:       '+results_folder+'\n'
  summary+='| target file:         '+target_file+'\n'
  summary+='| target species:      '+str(target_species)+'\n' # (taxid:'+str(target_species.taxid)+')\n'
  summary+='| profiles list:       '+join(profiles_names, ' ')+'\n'
  summary+='| configuration file:  '+config_filename+'\n'
  summary+='| temporary folder:    '+temp_folder+'\n'      
  # program options/ filtering / db
  summary+='|\n##########      Filtering procedures, program options, profile attributes defined:\n'
  for category in  sorted(keywords_text.keys()):
    summary+="| "+category.ljust(19)
    for k_index, keyword in enumerate(keywords_text[category]):
      if k_index>0: summary+="\n| "+' '*19
      summary+='| '+keyword+': '+str(keywords_text[category][keyword]).ljust(25)+' '
    summary+='\n'
  # actions

  summary+='|\n##########      Active actions:\n'
  for category in actions:
    for action_id in actions[category]:
      summary+=('| '+category+'.'+str(action_id)).ljust(19)+' =   '+actions[category][action_id]+'\n'
  if not actions:  summary+='| None\n'
  other_options=''
  for k in opt: 
    if not k in {'__synonyms__':1, 't':1, 'r':1, 'species':1, 'temp':1, 'config':1, 'profile':1} and not k.startswith('ACTION') and (not k.startswith('output_') or k.endswith('_file')) and (k!='state' or opt[k]!='kept') and (not k in def_opt or opt[k]!=def_opt[k]):
      other_options+='\n| '+k+': '+str(opt[k])+''
  if other_options: summary+='|\n##########      Routines and other non-default options:'+other_options+'\n'
  summary+='|\n##########      Output options:   '
  if output_options:
    for o in output_options: summary+=o.split('output_')[1]+' '
  else:    summary+='None'
  summary+='\n'+'#'*120

  return summary

def mask_species(species_name):     return replace( mask_characters (   species_name   ), ' ', '_')
def unmask_species(species_name):   return unmask_characters (  replace( species_name , '_', ' ') )

def load_chromosome_lengths(chromosome_length_file, max_chars=0):
  """Utility to load chromosome lenghts from a fastalength output file and also set it as a MMlib variable; also performing controls on the file """
  global chromosome_lengths; chromosome_lengths={}
  for line in open(chromosome_length_file, 'r'):     
    fasta_identifier = line.split()[1]
    length=int(line.split()[0])
    if chromosome_lengths.has_key(fasta_identifier): 
      bash('rm '+chromosome_length_file)
      raise notracebackException, "ERROR the target file has a duplicate fasta identifier! ("+line.split()[1]+') Please modify it and rerun. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if length==0: 
      bash('rm '+chromosome_length_file)
      raise notracebackException, "ERROR the target file has a length zero entry! ("+line.split()[1]+') Please modify it and rerun. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if is_number(fasta_identifier) and fasta_identifier[0]=='0':
      bash('rm '+chromosome_length_file)
      raise notracebackException, "ERROR the target file has a numeric fasta identifier starting with zero!  ("+line.split()[1]+') This would cause an unexpected blast behavior. Please modify this or these ids and rerun. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if ':subseq(' in fasta_identifier: 
      bash('rm '+chromosome_length_file)
      raise  notracebackException, "ERROR with fasta header: "+fasta_identifier+' ; this was generated by fastasubseq and will cause unexpected behavior of this program, since it is using fastasubseq itself to cut sequences. Please clean the titles in your target file from ":subseq(" tags. Note: remove the .index and *.fa.n* blast formatting files after changing the target file '
    if max_chars and   len(fasta_identifier)>max_chars: 
      bash('rm '+chromosome_length_file)
      raise  notracebackException, "ERROR with fasta header: "+fasta_identifier+' is too long. The maximum length for a fasta identifier (first word of the title) is '+str(max_chars)+' characters. Please clean the titles in your target file. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if '/' in fasta_identifier:
      bash('rm '+chromosome_length_file)
      raise  notracebackException, "ERROR with fasta header: "+fasta_identifier+' has forbidden character: "/" \nPlease clean the titles in your target file. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'

    chromosome_lengths[fasta_identifier]=length
  set_MMlib_var('chromosome_lengths', chromosome_lengths)

############################################################################################
################################################3
### MAIN PROGRAM START
def main():
  write('', 1);
  write('|         '); write('Selenoprofiles v'+str(__version__), how='reverse'); write('          Host: '+bbash('echo $HOSTNAME')+'   Date: '+bbash('date'),  1)
  options_summary= load()    #loading variables and getting a nice summary to show
  write('\n'+options_summary+'\n', 1)

  if opt['add']:                  # with -add  options, the user can add code to be executed right here. This can include the definition of function that will be used in the customized steps
    added_file=opt['add']
    check_file_presence(added_file, 'added code file')
    all_lines_of_added_file=open(added_file, 'r').readlines()
    try:      exec(join(all_lines_of_added_file, ''))
    except:   printerr('add ERROR importing code from file : '+ added_file); raise
    
  global max_chars_per_column;   global blast_nr_folder_profile_subfolder
  # computing length of all chromosomes in the target, or just loading it from chromosome_length_file
  write('Loading length of all chromosomes from '+chromosome_length_file, 1)
  load_chromosome_lengths(chromosome_length_file)
  for k in chromosome_lengths:
    if len(k)> max_chars_per_column['chromosome']: max_chars_per_column['chromosome']=len(k)

  all_results_are_loaded_from_db = True 
  chromosomes_with_results={} #filled after filtering, when writing in database
  skipped_profiles={}         # taken off the profiles_names list
  for profile_index, family in enumerate(profiles_names):
    try:
      # load profile
      profile_ali=profiles_hash[family]
      if len(family)+18>max_chars_per_column['id']:        max_chars_per_column['id']=len(family)+18

      write('\n'+profile_ali.summary(), 1, how=terminal_colors['profile'])
      
      ### important data structures to keep links to all loaded data:        
      blast_hits_hash={}     #original_blast_id : -> (super)blasthit_object
      exonerate_hits_hash={} #original_blast_id : -> exonerate_object
      genewise_hits_hash={}  #original_blast_id : -> genewise_object
      blast_hits_of_empty_exonerates_hash={} #index collection
      considered_indexes=[]             #index collection 
      chosen_predictions={}  #original_blast_it: ->  p2g_hit object chosen (blast, exonerate or genewise)
          
      if profile_ali.nseq()>1: 
        if not profile_ali.has_title('BLAST_QUERY_MASTER'):  # it could have it already if the profile has been built now
          profile_ali.add( 'BLAST_QUERY_MASTER', profile_ali.blast_queries().seq_of('BLAST_QUERY_MASTER')  )
        clusters_relative_positions = profile_ali.clusters_relative_positions() # to convert quickly positions on any cluster query to the corresponding position on the blast master query, we precompute all of them and put them into a hash.

      blast_folder_profile_subfolder=Folder(blast_folder+family)
      #check if this target has already been scanned for this profile. 
      if not opt['no_db']:      db_has_results_for_this_profile=   results_db.has_results_for_profile(family)
      if not opt['no_db'] and db_has_results_for_this_profile in ['ONGOING', 'UNCHECKED', 'WAITING']:
        raise notracebackException, "ERROR the search for profile: "+family+' is already being computed by another instance of selenoprofiles! If this is not true, use the script selenoprofiles_database.py to clean the file results.sqlite'     
      elif   opt['no_db']  or db_has_results_for_this_profile=='NO' or any( [opt[o] for o in ['blast', 'exonerate', 'genewise', 'choose', 'filter', 'database']] ):
        all_results_are_loaded_from_db= False
        if not opt['no_db']:
          results_db.set_has_results_for_profile(family, 'ONGOING')
          results_db.save()
        #########
        ## BLAST
        blast_parsers=[]
        write('\n'); write('BLAST:', how=terminal_colors['routine']); write(' psitblastn of profile '+family+' --> '+blast_folder_profile_subfolder+' ', 1)
        #### change for tblastn

        for cluster_index in range( profile_ali.n_clusters() ): #index is 0 based
          
          blast_outfile=       blast_folder_profile_subfolder+family+'.psitblastn.'+str(cluster_index+1)
          cluster_profile_ali= profile_ali.clusters()[cluster_index]
          ## debug
          #cluster_profile_ali.display(temp_folder+'cluster_profile_ali.fa')
          #profile_ali.display(temp_folder+'profile_ali.fa')
          #raw_input(temp_folder+'profile_ali.fa  '+temp_folder+'cluster_profile_ali.fa')

          if profile_ali.nseq()==1:
            write(  ('cluster'+str(cluster_index+1)).ljust(15)+(' (1 seq) --tblastn--').ljust(40) )
            if opt['blast'] or not is_file(blast_outfile) or not is_valid_blast_output(blast_outfile):
              write('R> '+blast_outfile, 1)       
              blast_parsers.append( tblastn(profile_ali, target_file, outfile=blast_outfile, blast_options=blast_options) )
            else:
              write('L> '+blast_outfile, 1)
              blast_parsers.append( parse_blast(blast_outfile) )        
          else:
            write(  ('cluster'+str(cluster_index+1)).ljust(15)+(' ('+str(cluster_profile_ali.nseq()-1)+' seq'+ 's'*int(cluster_profile_ali.nseq()>2) +')').ljust(24) )
            if opt['blast'] or not is_file(blast_outfile) or not is_valid_blast_output(blast_outfile):
              write(' R> '+blast_outfile, 1)
              if not opt['tblastn']:
                blast_parsers.append( psitblastn(cluster_profile_ali, target_file, outfile=blast_outfile, blast_options=blast_options) )
              else: 
                blast_parsers.append( multi_tblastn(cluster_profile_ali, target_file, outfile=blast_outfile, blast_options=blast_options) )
            else:
              write(' L> '+blast_outfile, 1)
              blast_parsers.append( parse_blast(blast_outfile) )
              
        if opt['filtered_blast_file']: filtered_blast_file_h=open(opt['filtered_blast_file'], 'w')
        
        ## filtering blast hits depending on the options on the profile 
        blast_hits=[] ;  blast_hit_index=1
        for cluster_id,  blast_hits_parser  in enumerate(blast_parsers):
          cluster_id+=1 #to have it 1-based
          blast_queries_alignment= alignment() #alignment of blast_query_master and of the query for this cluster. Used to speed up the process of blast_h.set_query_to_master (see below)
          blast_queries_alignment.add('q',  profile_ali.blast_queries().seq_of('BLAST_QUERY_'+str(cluster_id)))
          blast_queries_alignment.add('BLAST_QUERY_MASTER', profile_ali.blast_queries().seq_of('BLAST_QUERY_MASTER'))
          blast_queries_alignment.remove_useless_gaps()

          for blast_h in blast_hits_parser:
            if not chromosome_lengths.has_key(blast_h.chromosome) and is_number(blast_h.chromosome.split('_')[0]) and base_filename(target_file) in blast_h.chromosome: 
              raise notracebackException, "selenoprofiles ERROR the chromosome of the blast hit is not recognized. This is a known bug of blast which happens when the fasta titles in the target have \"|\" characters, but are not gis codes from ncbi. Please, reformat your target file (taking care also that all chromosome names are unique) and rerun selenoprofiles."
            if blast_hit_index> profile_ali.max_blast_hits_number_value():
              ##too many blast hits. now a little trick to avoid a painful error message that gawk will print to screen since we wouldn't get to the end of the blaster_parser
              done_end_of_parser=False
              while not done_end_of_parser: 
                try: blast_hits_parser.next()
                except StopIteration: done_end_of_parser=True
              if opt['blast_filtering_warning']:
                printerr("profile "+family+" WARNING too many blast hits! you should set a more strict blast filtering for this family. Now keep running using only with the first "+str(profile_ali.max_blast_hits_number_value())+' blast hits', 1)
                break
              else:
                raise skipprofileException, 'profile '+family+' ERROR too many blast hits! you should change the blast filtering for this family or set a higher value for parameter max_blast_hits in its profile configuration file. See blast filtering chapter in the manual for details. To keep running with the first max_blast_hits blast hits, use option -blast_filtering_warning'
            blast_h.target=target_file;                blast_h.profile=profile_ali;    blast_h.species=target_species; blast_h.cluster_id=cluster_id
            if profile_ali.nseq()>1: blast_h.set_query_to_master(clusters_relative_positions, blast_queries_alignment)
            ### CHANGED; now called later unless necessary (sec_is_aligned); it was:  blast_h.place_selenocysteine() ## modifying in place the blast hit to have Us for selenocysteine, also in the target. After this, there are U both in query and target (here only in UGAs aligned to U in query)
            x=blast_h
            for action_id in sorted(actions['pre_blast_filter'].keys()):
              try: exec(actions['pre_blast_filter'][action_id])  ### running action     
              except:
                printerr('selenoprofiles ERROR trying to run pre_blast_filter action: '+str(actions['pre_blast_filter'][action_id]) +' on result:\n '+str(x), 1)
                raise
            #write('CTL '+blast_h.chromosome+' '+blast_h.positions_summary()+' '+str(round(blast_h.weighted_seq_identity_with_profile(), 4)), 1)              
            if profile_ali.blast_filtering_eval()(blast_h): #filtering blast hits
              blast_h.place_selenocysteine()
              blast_h.id = str(blast_hit_index);         blast_h.query.id = str(blast_hit_index)+'_query'
              blast_hits.append(blast_h)
              blast_hit_index+=1
              if opt['filtered_blast_file']: print >> filtered_blast_file_h,  'BLASTID:'+blast_h.id+'\n', blast_h
              x=blast_h
              for action_id in sorted(actions['post_blast'].keys()):
                #write('id: '+str(index_id)+' checking pre action '+str(action_id)+' : '+actions['pre_filtering'][action_id], 1)
                try: exec(actions['post_blast'][action_id])  ### running action
                except:   
                  printerr('selenoprofiles ERROR trying to run post_blast action: '+str(actions['post_blast'][action_id]) +' on result with index '+str(blast_hit_index), 1)
                  raise
            
            ## sometimes blast gives an output changing the chromosome names to titles like "1_genome.fa" , where 1 is the index in the file where the fasta title appears, and genome.fa is the target name. This happens presumably because of unwanted characters in the fasta titles. user should remove them by himself
        
        if opt['filtered_blast_file']: filtered_blast_file_h.close()
            
        if profile_ali.n_clusters()>1:
          write("-- Removing blast hits overlapping from different clusters. before: "+str(len(blast_hits))+ " ; after: ", )
          blast_hits= remove_overlapping_genes(blast_hits, cmp_fn=choose_among_overlapping_p2gs_intrafamily,  phase=True)
          write(str(len(blast_hits)), 1)
          
        if not opt['no_splice']:
          # merging by colinearity blast hits, meaning: if A is downstream of B, and also the query of A is downstream of the query of B, they're merged.
          write("-- Merging by colinearity blast hits. before: "+str(len(blast_hits))+ " ; after: ", )
          superblast_hits = merge_p2g_hits_by_colinearity(  blast_hits, post_function=merge_p2g_hits_by_colinearity_post_function_blast_hits, sequence_collection=profile_ali )
          write(str(len(superblast_hits)), 1)          
          del blast_hits
        else: 
          superblast_hits = blast_hits
          write("-- Number of filtered blast hits: "+str(len(blast_hits)), 1)

        # superblasthits is actually a collection of both blasthits and superblasthits. You can use any blasthit method nonetheless
        for superblast_hit in superblast_hits:
          blast_hits_hash[superblast_hit.id] = superblast_hit
          considered_indexes.append(int(superblast_hit.id))
          x=superblast_hit
          for action_id in sorted(actions['post_blast_merge'].keys()):
            #write('id: '+str(index_id)+' checking pre action '+str(action_id)+' : '+actions['pre_filtering'][action_id], 1)
            
            try: exec(actions['post_blast_merge'][action_id])  ### running action
            except:   
              printerr('selenoprofiles ERROR trying to run post_blast_merge action: '+str(actions['post_blast_merge'][action_id]) +' on result with index '+str(x.id), 1)
              raise

        del superblast_hits #so I am deleting from memory those which were merged, but not the other ones cause they will be linked in blast_hits_hash
        considered_indexes.sort()

        ############## 
        ## EXONERATE
        if opt['dont_exonerate']:
          for hit_index in considered_indexes: blast_hits_of_empty_exonerates_hash[hit_index]=1
        else:
          exonerate_folder_profile_subfolder=Folder(exonerate_folder+family)
          write('\n'); write('EXONERATE:', how=terminal_colors['routine']) ; write(' exonerate for profile '+family+' --> '+exonerate_folder+' # R> stands for run, L> for load', 1)
          exonerate_mode={False:'p2g', True:'p2d'}[opt['no_splice']]
          for hit_index in considered_indexes:
            exonerate_outfile=exonerate_folder_profile_subfolder+family+'.'+str(hit_index)+'.exonerate'
            superblast_hit=   blast_hits_hash[str(hit_index)]
            if opt['exonerate'] or not is_file(exonerate_outfile): #deciding if run exonerate or not
              ### running cyclic_exonerate, obtaining an (super)exonerate object. It can be empty (in this case its boolean evaluation will be False, and it will have a "error_message" attribute).
              # the id of the exonerate hit is set to the id of the original blas hit here below.
              # to output, a short description is printed after the output file, like this:
              #sps.1.exonerate -> on:scaffold_6 strand:+ positions:1818342-1819187,1819199-18192912
              exonerate_hit = cyclic_exonerate(profile_ali, target_file, outfile=exonerate_outfile, seed=superblast_hit, extension=exonerate_extension, exonerate_options=exonerate_options, mode=exonerate_mode, merge_multiple=not opt['no_splice'])
              run_or_load_code='R'        
            else:
              #loading file
              exonerate_hit = superexoneratehit(); error_message=exonerate_hit.load(exonerate_outfile, seed=superblast_hit, merge_multiple=not opt['no_splice'])
              if not exonerate_hit: exonerate_hit.error_message= error_message
              run_or_load_code='L'

            exonerate_hit.id=superblast_hit.id; exonerate_hit.query.id= exonerate_hit.id+'_query' #setting id property of exonerate object to its hit_index
            exonerate_hit.profile=profile_ali;   exonerate_hit.target= target_file; exonerate_hit.species=target_species
           
            if exonerate_hit:      short_description= exonerate_hit.header(no_species=True, no_id=True, compress=True) 
            else:                     short_description=exonerate_hit.error_message
            write(description_format(exonerate_outfile.split('/')[-1]+' '+run_or_load_code+'> '+short_description), 1 )
            exonerate_hits_hash[superblast_hit.id] = exonerate_hit   
           
            if not exonerate_hit:
              #insert exhaustive exonerate code here ############
              blast_hits_of_empty_exonerates_hash[hit_index]=1

          if not considered_indexes:   write(' -- no blast hit passed filtering --', 1)
          else:
            duplicated_exonerate_hits=[] #this will contain the list of exonerate identical or overlapping to some other exonerate hit
            remove_overlapping_genes(exonerate_hits_hash.values(), scoring=len, phase=False, out_removed_genes=duplicated_exonerate_hits) #dropping output.. we're interested in the removed ones
            duplicated_exonerate_hits_ids=sorted([int(e.id) for e in duplicated_exonerate_hits])
            write('-- Removing duplicated exonerate hit'+'s'*int(len(duplicated_exonerate_hits_ids)>1)+' with id'+'s'*int(len(duplicated_exonerate_hits_ids)>1)+' : ')
            if duplicated_exonerate_hits_ids:
              write(join([str(i) for i in duplicated_exonerate_hits_ids], ', '), 1)
              remove_items_from_list(considered_indexes, duplicated_exonerate_hits_ids, inplace=True) #removing duplicates from considered_indexes
            else: write('(None)', 1)
        ###################
        ## GENEWISE
        if opt['dont_genewise']:        pass
        else:
          genewise_folder_profile_subfolder = Folder(genewise_folder+family)
          write( '\n'); write('GENEWISE:', how=terminal_colors['routine']); write(' genewise for profile '+family+' --> '+genewise_folder+' # R> stands for run, L> for load', 1)
          if not considered_indexes: write(' -- no blast hit passed filtering --', 1)
          for i_i, hit_index in enumerate(considered_indexes): #i_i is the index of the index... it is not used.
            if not hit_index in blast_hits_of_empty_exonerates_hash: #it means we have a non-empty exonerate output for this hit_index
              genewise_outfile=genewise_folder_profile_subfolder+family+'.'+str(hit_index)+'.genewise' ; current_outfile=genewise_outfile
              exonerate_seed_hit= exonerate_hits_hash[str(hit_index)]
              
              if opt['genewise'] or not is_file(genewise_outfile): #### run genewise!
                genewise_hit = genewise(profile_ali, target_file, outfile=genewise_outfile, seed=exonerate_seed_hit, extension=genewise_extension, genewise_options=genewise_options)
                run_or_load_code='R'
              else: # or load genewise file
                genewise_hit = genewisehit(genewise_outfile)
                run_or_load_code='L'
              genewise_hit.id = str(hit_index); genewise_hit.query.id = genewise_hit.id+'_query' #setting id property of genewise object to its hit_index
              genewise_hits_hash[genewise_hit.id] = genewise_hit
              genewise_hit.profile=profile_ali
              genewise_hit.species=target_species
              genewise_hit.target=target_file
              if genewise_hit:           short_description= genewise_hit.header(no_species=True, no_id=True, compress=True) 
              else:                      short_description=genewise_hit.error_message
              
            else: # genewise to be sure. the correspondent exonerate output was empty
              if opt['genewise_to_be_sure']:
                genewise_tbs_outfile=genewise_folder_profile_subfolder+family+'.'+str(hit_index)+'.genewise_tbs' ; current_outfile=genewise_tbs_outfile 
                blast_seed_hit=blast_hits_hash[str(hit_index)]
                
                if opt['genewise'] or not is_file(genewise_tbs_outfile): #### run genewise tbs!
                  genewise_hit = genewise(profile_ali, target_file, outfile=genewise_tbs_outfile, seed=blast_seed_hit, extension=genewise_tbs_extension, genewise_options=genewise_tbs_options)
                  run_or_load_code='R'
                else:            #or load genewise tbs file
                  genewise_hit = genewisehit(genewise_tbs_outfile)
                  run_or_load_code='L'
                genewise_hit.id = blast_seed_hit.id; genewise_hit.query.id = genewise_hit.id+'_query' #setting id property of genewise object to its hit_index
                genewise_hits_hash[genewise_hit.id] = genewise_hit
                genewise_hit.profile=profile_ali
                genewise_hit.species=target_species
                genewise_hit.target=target_file
                
                if genewise_hit:           short_description= genewise_hit.header(no_species=True, no_id=True, compress=True) 
                else:                      short_description=genewise_hit.error_message

            #both in tbs and in normal routine: print summary
            if opt['genewise_to_be_sure'] or not hit_index in blast_hits_of_empty_exonerates_hash:
              genewise_hit.check_alignment  ()
              write(description_format(current_outfile.split('/')[-1]+' '+run_or_load_code+'> '+short_description), 1)               

        write('', 1)
        ###############
        ## choosing prediction and labelling 
        chosen_predictions_reasons_why={}
        chosen_predictions_outfile=prediction_choice_folder+family+'.tab'
        #running pre_choose actions. This must be here, and not inside the next if, to ensure they are executed even if the prediction choose file has already been produced
        for index_id in considered_indexes:
          exonerate_hit=False; genewise_hit=False; blast_hit=False 
          blast_hit=     blast_hits_hash[str(index_id)]
          if exonerate_hits_hash.has_key(str(index_id)):      exonerate_hit= exonerate_hits_hash[str(index_id)]
          if genewise_hits_hash.has_key(str(index_id)):       genewise_hit=  genewise_hits_hash[str(index_id)]
          candidates=[blast_hit]  
          if exonerate_hit:  candidates.append(exonerate_hit)
          if genewise_hit:   candidates.append(genewise_hit)
          for action_id in sorted(actions['pre_choose'].keys()):
            for x in candidates:
              try: exec(actions['pre_choose'][action_id])  ### running action
              except:   
                printerr('selenoprofiles ERROR trying to run pre_choose action: '+str(actions['pre_choose'][action_id]) +' on result with index '+str(index_id)+' and prediction program '+x.prediction_program(), 1)
                raise

        write('CHOOSE:', how=terminal_colors['routine']); write(' choosing among available predictions, assigning label --> '+chosen_predictions_outfile+' ')
        indexes_to_remove={} ## keep track of the predictions we want to drop immediately after the choose step. uniq example on creation date: when the only available prediction is from blast, and you don't want blast predictions according to options

        if opt['choose'] or not is_file(chosen_predictions_outfile):
          write('RUN', 1)
          chosen_predictions_outfile_text=''
          if not considered_indexes: write(' -- no blast hit passed filtering --', 1)
          for index_id in considered_indexes:
            #recovering the three predictions
            exonerate_hit=False; genewise_hit=False; blast_hit=False 
            blast_hit=     blast_hits_hash[str(index_id)]
            if exonerate_hits_hash.has_key(str(index_id)):      exonerate_hit= exonerate_hits_hash[str(index_id)]
            if genewise_hits_hash.has_key(str(index_id)):       genewise_hit=  genewise_hits_hash[str(index_id)]

            candidates=[blast_hit]   
            if opt['no_blast']:   candidates=[]
            if exonerate_hit:  candidates.append(exonerate_hit)
            if genewise_hit:   candidates.append(genewise_hit)
              
            chosen_hit, reason_why =  choose_prediction(candidates)  ## choosing prediction here. see choose_prediction function for details
            chosen_predictions[str(index_id)]=chosen_hit
            chosen_predictions_reasons_why[str(index_id)]=reason_why

            if not chosen_hit:               
              label='empty'
              indexes_to_remove[index_id]=True
            else:
              #assigning label
              if profile_ali.sec_pos():        label= assign_label_selenoprotein_families(chosen_hit)
              else:                            label= assign_label_non_selenoprotein_families(chosen_hit)
            chosen_hit.label=label ## assigning label property here
            chosen_predictions_outfile_text+=str(index_id)+'\t'+chosen_hit.prediction_program()+'\t'+reason_why+'\t'+label+'\n'
                    
          # writing file with results from this step
          if chosen_predictions_outfile_text: chosen_predictions_outfile_text=chosen_predictions_outfile_text[:-1] #removing last \n
          write_to_file(chosen_predictions_outfile_text, chosen_predictions_outfile)
            
        else:
          write('(just loading file)', 1)
          considered_indexes_check_hash={}
          if not considered_indexes: write(' -- no blast hit passed filtering --', 1)
          for index_id in considered_indexes: considered_indexes_check_hash[index_id]=False
          for line in open(chosen_predictions_outfile, 'r'):
            if line.split():
              try:
                index_id=int(line[:-1].split("\t")[0])
                program_chosen=line[:-1].split('\t')[1]
                reason_why=line[:-1].split('\t')[2]
                label=line[:-1].split('\t')[3]
                if program_chosen=='blast':           chosen_predictions[str(index_id)]= blast_hits_hash[str(index_id)]
                elif program_chosen=='exonerate':     chosen_predictions[str(index_id)]= exonerate_hits_hash[str(index_id)]
                elif program_chosen=='genewise':      chosen_predictions[str(index_id)]= genewise_hits_hash[str(index_id)]                        
                elif program_chosen=='None' and label=='empty':      
                  indexes_to_remove[index_id]=True 
                  chosen_predictions[str(index_id)]= empty_p2g() ## just for pretty printing, then it's removed
                else: raise notracebackException, 'ERROR loading '+chosen_predictions_outfile+' : program "'+program_chosen+'" not recognized! '
                chosen_predictions[str(index_id)].label=label ## assigning label property here
                chosen_predictions_reasons_why[str(index_id)]=reason_why
                considered_indexes_check_hash[index_id]=True
              except KeyError: raise notracebackException, "ERROR indexes loaded from chosen prediction file "+chosen_predictions_outfile+" are not consistent with predictions indexes loaded in memory! Run with option -C to repeat the prediction choice"  
          if not all( considered_indexes_check_hash.values() ):  
            indexes_not_found_in_file=[] #strings
            for index_id in considered_indexes_check_hash: 
              if not considered_indexes_check_hash[index_id]: indexes_not_found_in_file.append(str(index_id))
            raise notracebackException,  "ERROR loading "+chosen_predictions_outfile+' some hit indexes were not found: '+join(indexes_not_found_in_file, ' ')
        
        # printing information to screen. Also, adding profile and species property to chosen predictions
        for index_id in considered_indexes:
          chosen_hit=chosen_predictions[str(index_id)]
          reason_why=chosen_predictions_reasons_why[str(index_id)]
          write( (family+'.'+str(index_id)).ljust(max_chars_per_column['id'])+' '+chosen_hit.prediction_program().ljust(15)+reason_why.ljust(34)+chosen_hit.label, 1)
          
        if indexes_to_remove:          
          considered_indexes = [ index for index in considered_indexes if  not indexes_to_remove.has_key(index)]
          for i in indexes_to_remove: del chosen_predictions[str(i)]
        del chosen_predictions_reasons_why #useless afterwards

        blast_nr_folder_profile_subfolder=Folder(blast_nr_folder+family) #setting this to global variable so it can be read from the blast_against_tag_db method in the p2ghit class.
        
        #############
        ### prefiltering ACTIONS (none by default)
        for index_id in considered_indexes:
          x=chosen_predictions[str(index_id)]
          for action_id in sorted(actions['pre_filtering'].keys()):
            #write('id: '+str(index_id)+' checking pre action '+str(action_id)+' : '+actions['pre_filtering'][action_id], 1)
            try: exec(actions['pre_filtering'][action_id])  ### running action
            except:   
              printerr('selenoprofiles ERROR trying to run pre-filtering action: '+str(actions['pre_filtering'][action_id]) +' on result with index '+str(index_id), 1)
              raise
        
        ###############
        ## filtering prediction
        filtering_predictions_outfile=filtered_list_folder+family+'.tab'
        write('\n'); write('FILTER:', how=terminal_colors['routine']); write(' removing redundancy of results and filtering --> '+filtering_predictions_outfile+' ')
        filtered_out_ids_hash={}    #this keeps the ids of filtered out hits
        refiltered_out_ids_hash={}  #this keeps the ids of refiltered out hits
        duplicated_chosen_hits_ids_hash={} #so these are skipped when looping through results
        if opt['filter'] or not is_file(filtering_predictions_outfile):
          write('RUN', 1)
          duplicated_chosen_hits=[] #this will contain the list of results identical or overlapping to some other 
          remove_overlapping_genes(chosen_predictions.values(), cmp_fn=choose_among_overlapping_p2gs_intrafamily, phase=False, out_removed_genes=duplicated_chosen_hits, remember_overlaps=True) #dropping output.. we're interested in the removed ones
          for d in duplicated_chosen_hits:     
            duplicated_chosen_hits_ids_hash[ int(d.id) ] = d  #annotating duplicated hits ids in this hash
            chosen_predictions[d.id] = d  #replacing the hit with an identical one, except for the .overlapping attribute which tells us which gene was overlapping with this gene

          filter_decision='';         filtering_predictions_outfile_text=''
          for index_id in considered_indexes:
            chosen_hit=chosen_predictions[str(index_id)]
            if  duplicated_chosen_hits_ids_hash.has_key(index_id):
              filtering_predictions_outfile_text+=str(index_id)+'\tredundant\t'+str(duplicated_chosen_hits_ids_hash[index_id].overlapping.id)+'\n'
              chosen_hit.filtered='redundant'
            else:
              passes_filtering=False
              try:
                passes_filtering=profile_ali.p2g_filtering_eval()(chosen_hit) # execute refiltering !! 
                if not passes_filtering:
                  filtered_out_ids_hash[index_id]=chosen_hit
                  filtering_predictions_outfile_text+=str(index_id)+'\tfiltered\n'
                  chosen_hit.filtered='filtered'
                else:
                  ## passing filtering. Let's check for refiltering.
                  passes_refiltering=False
                  passes_refiltering=profile_ali.p2g_refiltering_eval()(chosen_hit)
                  if not passes_refiltering:
                    filtering_predictions_outfile_text+=str(index_id)+'\trefiltered\n'          
                    refiltered_out_ids_hash[index_id]=chosen_hit
                    chosen_hit.filtered='refiltered'
                  else:
                    filtering_predictions_outfile_text+=str(index_id)+'\tkept\n'          
                    chosen_hit.filtered='kept'
              except Exception, e: 
                printerr('filtering/refiltering ERROR on '+chosen_hit.output_id()+' : '+str(sys.exc_info()[1]), 1)
                raise  

          # writing file with results from this step          
          if filtering_predictions_outfile_text: filtering_predictions_outfile_text=filtering_predictions_outfile_text[:-1] #removing last \n
          write_to_file(filtering_predictions_outfile_text, filtering_predictions_outfile)

        else:
          write('(just loading file)', 1)
          considered_indexes_check_hash={}
          for index_id in considered_indexes: considered_indexes_check_hash[index_id]=False
          for line in open(filtering_predictions_outfile, 'r'):
            if not line.rstrip(): break 
            index_id=int(line[:-1].split("\t")[0])
            filter_decision=line[:-1].split('\t')[1]

            chosen_predictions[str(index_id)].filtered=filter_decision # updating "filtered" attribute
            if filter_decision=='redundant':  # if the results was called redundant, let's set the .overlapping attribute pointing to the gene (overlapping to this) that was kept 
              chosen_predictions[str(index_id)].overlapping=chosen_predictions[line[:-1].split('\t')[2]]
              duplicated_chosen_hits_ids_hash[index_id]=chosen_predictions[str(index_id)]
            elif filter_decision=='filtered':     filtered_out_ids_hash[index_id]=chosen_predictions[str(index_id)]
            elif filter_decision=='refiltered':   refiltered_out_ids_hash[index_id]=chosen_predictions[str(index_id)]
            elif filter_decision=='kept':         'nothing to do'
            else: raise notracebackException, 'ERROR loading '+chosen_predictions_outfile+' : program "'+program_chosen+'" not recognized! '
            considered_indexes_check_hash[index_id]=True

          if not all( considered_indexes_check_hash.values() ):  
            indexes_not_found_in_file=[] #strings
            for index_id in considered_indexes_check_hash: 
              if not considered_indexes_check_hash[index_id]: indexes_not_found_in_file.append(str(index_id))
            raise notracebackException,  "ERROR loading "+filtering_predictions_outfile+' some hit indexes were not found: '+join(indexes_not_found_in_file, ' ')

        if not considered_indexes: write(' -- no blast hit passed filtering --', 1)        
        if not opt['no_db']:        results_db.clear_results(family)
    
        #############
        ### post_filtering ACTIONS: including SECISEARCH 
        for action_id in sorted(actions['post_filtering'].keys()):
          for index_id in considered_indexes:
            x=chosen_predictions[str(index_id)]
            #write('checking action '+str(action_id)+' : '+actions['post_filtering'][action_id], 1)
            try: exec(actions['post_filtering'][action_id])  ### running action
            except:   
              printerr('selenoprofiles ERROR trying to run post_filtering action: '+str(actions['post_filtering'][action_id]) +' on result with index '+str(index_id), 1)
              raise
        
        ########################
        ####  DATABASE write down everything
        # printing information to screen and saving to database
        for index_id in considered_indexes:
          chosen_hit=chosen_predictions[str(index_id)]
          chromosomes_with_results[chosen_hit.chromosome]=1
          if duplicated_chosen_hits_ids_hash.has_key(index_id):  write( chosen_hit.output_id().ljust(max_chars_per_column['id']+3)+'DROPPED  since overlaps with '+chosen_hit.overlapping.output_id(), 1)        
          elif filtered_out_ids_hash.has_key(index_id):          write( chosen_hit.output_id().ljust(max_chars_per_column['id']+3)+'DROPPED  by filtering', 1)   
          else:                                                                                           
            write( (family+'.'+str(index_id)+'.'+chosen_hit.label).ljust(max_chars_per_column['id']+3)+'KEPT     by filtering    | ')
            if refiltered_out_ids_hash.has_key(index_id):                        write('DROPPED     by refiltering', 1)
            else:                                                                write('KEPT        by refiltering    |--> OK', 1)
          if not opt['no_db']:
            ### adding result to database
            if opt['full_db'] or chosen_hit.filtered=='kept':            
              chosen_hit.compute_misc_sequence_data()         ## computing miscellaneous sequence data to be included in the database as feature
              results_db.add_result(chosen_hit)            
        if not opt['no_db']:  
          results_db.set_has_results_for_profile(family, 'UNCHECKED'); 
          results_db.save()
          write('\n'); write('DATABASE', how=terminal_colors['database']); write(' wrote results for profile: '+family, 1)
          
        remove_items_from_list(considered_indexes, duplicated_chosen_hits_ids_hash.keys()+filtered_out_ids_hash.keys()+refiltered_out_ids_hash.keys(), inplace=True)  #removing duplicates and filtered out from considered_indexes
      elif db_has_results_for_this_profile=='YES' or opt['merge']:
        write('\n'); write('DATABASE', how=terminal_colors['database']); write(' loading results! skipping all pipeline...', 1)    
      elif db_has_results_for_this_profile in ['UNCHECKED-2', 'WAITING'] and not opt['merge']:
        write('\n'); write('DATABASE', how=terminal_colors['database']); write(' WARNING loading results that were not checked for inter-family redundancy! consider rerunning with option -merge', 1)

    except skipprofileException:
      skipped_profiles[family]=True
      if not opt['no_db']:
        results_db.set_has_results_for_profile(family, 'NO')   ;         results_db.save()    
      printerr( sys.exc_info()[1], 1)
    except: 
      # something went wrong. here we make sure that the db doesn't say that some profile is still computing
      if not opt['no_db']:
        results_db.set_has_results_for_profile(family, 'NO')   ;         results_db.save()          
        for completed_profile_name in profiles_names[:profile_index]:
          if results_db.has_results_for_profile(completed_profile_name)=='UNCHECKED': results_db.set_has_results_for_profile(completed_profile_name, 'UNCHECKED-2')
        results_db.save()
      raise

    ## cleaning memory
    del blast_hits_hash; del exonerate_hits_hash; del genewise_hits_hash; del blast_hits_of_empty_exonerates_hash; del considered_indexes; del chosen_predictions

  if skipped_profiles:
    for profile_index in range(  len(profiles_names)-1, -1, -1  ):
      if  profiles_names[profile_index] in skipped_profiles:   profiles_names.pop(profile_index)

  ### managing workflow depending on db options
  if opt['no_db']:
    write('\nOption -no_db:   the files are now ready, but results are not stored in the database and were not checked for inter-family redundancy. Please run with option -database. Now quitting...', 1)            
    sys.exit()  
  if opt['stop']:  ### if this option was specified, we don't go to the output phase since we expect more profiles to be computed, and the redundancy must be checked before outputing
    for family in profiles_names:        results_db.set_has_results_for_profile(family, 'UNCHECKED-2');
    results_db.save()
    write('\nOption -stop:   results are now stored in the database, but they were not checked for inter-family redundancy. Now quitting...', 1)            
    sys.exit()
  
  #### merging inter-profile results, directly on the database. but we have to check for other instances of selenoprofiles
  if not all_results_are_loaded_from_db or opt['merge']:
    if not opt['merge'] and  not opt['no_merge'] and results_db.computation_in_progress(apart_from=profiles_names):
      # there's a process computing. this (self) program will just wait until that process will finish and check redundancy. then, it will pass to the output phase
      try:
        for family in profiles_names:        results_db.set_has_results_for_profile(family, 'WAITING')
        results_db.save()  
        write('DATABASE', how=terminal_colors['database']);  write(" detected other(s) selenoprofiles instance(s) computing on this target. Waiting for the last one to end ...", 1)
        index_waited=0
        while results_db.computation_in_progress():
          index_waited+=1
          time.sleep(sleep_time)
          service('waited: '+str(sleep_time*index_waited)+' seconds ')
      except: 
        for family in profiles_names:        results_db.set_has_results_for_profile(family, 'UNCHECKED-2')
        results_db.save()        
        raise
      for family in profiles_names:        results_db.set_has_results_for_profile(family, 'YES')
      results_db.save()      
      write('DATABASE', how=terminal_colors['database']);  write(" Proceeding to output phase", 1)
    elif not opt['no_merge']: 
      ## if opt['merge']
      # this is the last active selenoprofiles. let's check overlaps and unlock so every possible process waiting will proceed
      write('DATABASE', how=terminal_colors['database']); write(' checking if results from different profiles are overlapping...', 1)
      try:
        if opt['merge']:    chromosomes=None # will trigger processing all chromosomes
        else:               chromosomes=chromosomes_with_results.keys()
        #################
        ###checking if any result is overlapping with anyone else in the database... for those, it prints a message to screen, and the state is set to overlapping
        results_db.remove_redundancy( chromosome_list=chromosomes )   
        for family in profiles_names:       results_db.set_has_results_for_profile(family, 'YES')
        results_db.set_stopped_profiles_as_ok()  ## if we run instances of selenoprofiles with option -stop, we can now set them to YES, as they were checked for redundancy
        results_db.save()  
      except: 
        for family in profiles_names:        results_db.set_has_results_for_profile(family, 'UNCHECKED-2')
        results_db.save()        
        raise        
    else: 
      #opt['no_merge']
      write('Option -no_merge: proceeding to output phase without checking if from different profiles are overlapping...', 1)      
      
  if opt['clean']:
    # if this option is active, we delete all intermediate files now that results are stored in the database.
    write('CLEAN:', how=terminal_colors['routine']); write(' removing intermediate files...', 1)
    for family in profiles_names:
      profile_ali=profiles_hash[family]
      write('Cleaning files for profile: '+family, 1)
      bash('rm -r '+blast_folder+family+' '+exonerate_folder+family+' '+genewise_folder+family+' '+prediction_choice_folder+family+'.tab '+filtered_list_folder+family+'.tab ')#+blast_nr_folder+family)
  
  write('\n'); write("="*56); write(' OUTPUT ', how=terminal_colors['routine']); write("="*56, 1)
 ###############
 ## OUTPUT predictions # I cycle through profiles again, cause after filtering we need to check hits overlapping with some other family
  output_something=False
  for family in profiles_names:
    profile_ali=profiles_hash[family]
    profile_ali_with_predictions=profile_ali.copy()

    if not 'blast_nr_folder_profile_subfolder' in globals(): blast_nr_folder_profile_subfolder=Folder(blast_nr_folder+family) #setting this to global variable so it can be read from the blast_against_tag_db method in the p2ghit class.
         
    #write( '>>> profile: '+family, 1)
    output_indexes=[]                       #index collection of the results which have active state
    all_predictions={}                      #original_blast_it id:str ->  p2g_hit object chosen (blast, exonerate or genewise)
    kept_predictions={}                        #original_blast_it id:str -> hits kepts by filtering/refiltering/redundancy filter
    duplicated_chosen_hits_ids_hash={}         #original_blast_it id:str -> hits removed after filtering because redundant
    filtered_out_ids_hash={}                   #original_blast_it id:str -> hits removed with the filtering 
    refiltered_out_ids_hash={}                 #original_blast_it id:str -> hits removed with the refiltering
    duplicated_in_other_families_ids_hash={}   #original_blast_it id:str -> hits removed with remove_redundancy method of the db

    #determining active states. These arethe filtreing decisions for which we decide to give output. So for example if the user provide -state filtered, only those predictions filtered away will be processed.   Multiple states can be provided by separating with commas, e.g. :   -state kept,refiltered 
    states_hash={'kept':kept_predictions, 'redundant':duplicated_chosen_hits_ids_hash, 'filtered': filtered_out_ids_hash, 'refiltered':   refiltered_out_ids_hash, 'overlapping':duplicated_in_other_families_ids_hash}
    active_states=opt['state'].split(',')
    if opt['state']=='all': active_states=states_hash.keys()
    if opt['output_filter']:             
      profile_ali_with_predictions_file= output_folder+family+'.custom_filter.ali'
      active_states=[]
    elif active_states == ['kept']:      profile_ali_with_predictions_file= output_folder+family+'.ali'
    else:                                profile_ali_with_predictions_file= output_folder+family+'.'+join(active_states, ',')+'.ali'
   
    for db_result in results_db.get_results(family, active_states):
      chosen_hit = load_p2g_from_db_result(db_result, profile=profile_ali)
      state= chosen_hit.filtered
      if opt['output_filter']:               
        x=chosen_hit
        try: 
          if eval(opt['output_filter']):         #running filter
            output_indexes.append(int(chosen_hit.id))           
            all_predictions[chosen_hit.id]=chosen_hit
        except:    printerr('selenoprofiles ERROR trying to run output_filter: '+str(opt['output_filter']) +' on result with index '+str(chosen_hit.id), 1);  raise
      else:    all_predictions[chosen_hit.id]=chosen_hit

    ###### check if they are already sorted
    output_indexes=  [str(j) for j in     sorted( [int(i) for i in all_predictions.keys()])     ]   #ordering indexes numerically
 
    for index_id in output_indexes:
      chosen_hit=all_predictions[str(index_id)] 
      prefix_outfile=output_folder+family+'.'+str(index_id)+'.'+chosen_hit.label # without final dot
      
      #############
      ### pre_output ACTIONS
      for action_id in sorted(actions['pre_output'].keys()):
        x=all_predictions[str(index_id)]
        try: exec(actions['pre_output'][action_id])  ### running action
        except:   
          printerr('selenoprofiles ERROR trying to run pre_output action: '+str(actions['pre_output'][action_id]) +' on result with index '+str(index_id), 1)
          raise
     
      if opt['output_ali']:
        #adding sequence to global alignment
        profile_ali_with_predictions=   chosen_hit.alignment_with_profile(profile_ali=profile_ali_with_predictions, dont_shrink=True)
      for keyword in allowed_output_formats:
        this_outfile=prefix_outfile+'.'+keyword # e.g. selenoprofiles_results/Homo_sapiens/output/sps_temp.25.threonine.fa
        if (opt['output_'+keyword] and (not is_file(this_outfile) or opt['O'])) or opt['output_'+keyword+'_file']:
          textout=getattr(chosen_hit, 'output_'+keyword)() # like it was:             chosen_hit.fasta() 
          if not textout is None:        
            write( chosen_hit.output_id().ljust(max_chars_per_column['id']+3)+(' -'+upper(keyword)+'> ').ljust(4+max([len(i) for i in allowed_output_formats])) ) # e.g. >sps.25.threonine -FASTA> 
            output_something=True
            if opt['output_'+keyword+'_file']:
              write( opt['output_'+keyword+'_file']+'  ') 
              print >> output_file_handlers[keyword], textout #like it was           print >> output_fasta_file_h, text
            if  opt['output_'+keyword] and (not is_file(this_outfile) or opt['O']) :
              write(this_outfile)
              write_to_file( textout, this_outfile)
            write(' ',1)

    if opt['output_ali'] and output_indexes and (not is_file(profile_ali_with_predictions_file) or opt['O'] or opt['output_filter']): #deciding if to write .ali ; if output_filter is active, we write a .custom_filter.ali file and this should be done anyway, even if the user produced another before
      output_something=True
      prediction_titles= profile_ali_with_predictions.titles()[  profile_ali.nseq(): ]
      profile_ali_with_predictions.shrink(only_titles=prediction_titles  )
      write( (family+'.*').ljust(max_chars_per_column['id']+3)+(' -ALI> ').ljust(4+max([len(i) for i in allowed_output_formats]))+profile_ali_with_predictions_file, 1)
      profile_ali_with_predictions.display(profile_ali_with_predictions_file)
      
  if not output_something: write(' -- nothing to output -- ', 1)
  write('\nPipeline workflow completed.   Date: '+bbash('date'), 1)    
  
  #compute global alignments for indexes in     indexes_passing_p2g_filter     and also only for those     indexes_passing_p2g_refilter

     
  #remove every file created for this family in the temp folder
     
#######################################################################################################################################

### main routines functions! They are used later on.
def psitblastn(profile, target_file, outfile='', blast_options={}):
  """ This function runs psitblastn on the target and store results on outfile (if indicated) or in a temporary file, and returns a blast_parser to the results.
      Blast options are read first from blast_options and then from the profile alignment (latter overriding the former).
  """  
  blast_options_used=blast_options.copy()
  profile_blast_options=profile.blast_options_dict()
  for option_name in profile_blast_options:     blast_options_used[option_name]=profile_blast_options[option_name]
  if not outfile:                  outfile=fileid_for_temp_folder(profile.filename)+'_BLAST_'+fileid_for_temp_folder(target_file)

  b=bash('which blastall')
  if b[0]: raise notracebackException, "ERROR blastall not found! Please install it"
  blastall_bin= dereference( b[1] )
  bbash(blastall_bin+' -p psitblastn  -d '+target_file+' -R '+profile.pssm()+ ' -i '+profile.blast_query_file()+' -I  '+join(['-'+k+' '+str(blast_options_used[k]) for k in blast_options_used], ' ')+' > '+temp_folder+'psitblastn_out')
  bbash('mv '+temp_folder+'psitblastn_out '+outfile)
  if not is_valid_blast_output(outfile): raise notracebackException, "psitblastn ERROR the blast output "+outfile+" doesn't seem complete, although blast exit status was not an error! check the options used"
  return parse_blast(outfile)

def tblastn(ss_profile, target_file, outfile='', blast_options={}):
  """ This function runs tblastn with a single sequence profile on the target and store results on outfile (if indicated) or in a temporary file, and returns a blast_parser to the results.
      Blast options are read first from blast_options and then from the profile alignment (latter overriding the former).
  """  
  blast_options_used=blast_options.copy()
  profile_blast_options=ss_profile.blast_options_dict()
  for option_name in profile_blast_options:     blast_options_used[option_name]=profile_blast_options[option_name]
  if not outfile:                  outfile=fileid_for_temp_folder(ss_profile.filename)+'_BLAST_'+fileid_for_temp_folder(target_file)
  b=bash('which blastall')
  if b[0]: raise notracebackException, "ERROR blastall not found! Please install it"
  blastall_bin= dereference( b[1] )
  bbash(blastall_bin+' -p tblastn  -d '+target_file+' -i '+ss_profile.filename+' -I  '+join(['-'+k+' '+str(blast_options_used[k]) for k in blast_options_used], ' ')+' > '+temp_folder+'tblastn_out')
  bbash('mv '+temp_folder+'tblastn_out '+outfile)
  if not is_valid_blast_output(outfile): raise notracebackException, "tblastn ERROR the blast output "+outfile+" doesn't seem complete, although blast exit status was not an error! check the options used"
  return parse_blast(outfile)


def multi_tblastn(ms_profile, target_file, outfile='', blast_options={}):
  """ This function runs a tblastn with a multiple sequence profile compressed to a single consensus query sequence, on the target and store results on outfile (if indicated) or in a temporary file, and returns a blast_parser to the results.
      Blast options are read first from blast_options and then from the profile alignment (latter overriding the former).
  """  
  blast_options_used=blast_options.copy()
  profile_blast_options=ms_profile.blast_options_dict()
  for option_name in profile_blast_options:     blast_options_used[option_name]=profile_blast_options[option_name]
  if not outfile:                  outfile=fileid_for_temp_folder(ms_profile.filename)+'_BLAST_'+fileid_for_temp_folder(target_file)
  b=bash('which blastall')
  if b[0]: raise notracebackException, "ERROR blastall not found! Please install it"
  blastall_bin= dereference( b[1] )
  bbash(blastall_bin+' -p tblastn  -d '+target_file+' -i '+ms_profile.blast_query_file()+' -I  '+join(['-'+k+' '+str(blast_options_used[k]) for k in blast_options_used], ' ')+' > '+temp_folder+'tblastn_out')
  bbash('mv '+temp_folder+'tblastn_out '+outfile)
  if not is_valid_blast_output(outfile): raise notracebackException, "tblastn ERROR the blast output "+outfile+" doesn't seem complete, although blast exit status was not an error! check the options used"
  return parse_blast(outfile)

def exonerate(query_file, target_file, outfile='', exhaustive=False, exonerate_options={}, dont_parse=False, mode='p2g'): 
  """ see below in cyclic exonerate; note: in selenoprofiles, dont_parse=True """
  if not outfile:  outfile=temp_folder+'tempout.exonerate'
  exonerate_options_str = ''
  for k in exonerate_options: exonerate_options_str+=' -'+k+' "'+exonerate_options[k]+'"'
  b=bash('exonerate -m '+mode+' --showtargetgff 1 '+exonerate_options_str+' '+int(exhaustive)*'--exhaustive '+' "'+query_file+'" "'+target_file+'" > '+outfile)
  if b[0]:
    #an error occured
    if "Expected protein query (not DNA) for model" in b[1]: 
      printerr('WARNING exonerate crashed because the protein sequence looks too much like a nucleotide sequence. Skipping this result!', 1)
    else: raise Exception, "ERROR running exonerate: " +b[1]
  if not dont_parse:  return parse_exonerate(outfile)

def cyclic_exonerate(profile_ali, target_file, outfile='', seed='', extension=30000, exonerate_options={}, max_iterations=10, chromosome_lengths_hash={}, mode='p2g', merge_multiple=True):
  """ This function runs cyclic_exonerate (see Mariotti and Guigo, 2010) on the target using the specified profile and using as seed (determining the positions of search) a gene object, whose boundaries are used. 
  If outfile is not indicated, it is put on a temporary file.   A exonerate_parser object is returned.
  seed is a gene object which contains all the information to initiate the cyclic exonerate routine. the chromosome name, strand and boundaries positions are used. If a .query attribute is present, its .chromosome attribute is used to obtain the full query sequence and so use the profile to choose the best query before running exonerate. If it hasn't, the first query chosen is simply the first sequence of the alignment.
  If the exonerate object is empty or it has load error, an empty exonerate object with a descriptive .error_message attribute is returned
  """
  global chromosome_lengths
  if not chromosome_lengths_hash: chromosome_lengths_hash=chromosome_lengths
 
  #replacing the input profile with a copy which change in all sequences all aminoacids in sec columns with *
  profile=profile_ali.copy()
  
  #determining exonerate options: the ones specified for the profile override the ones specified as arguments of cyclic_exonerate
  exonerate_options_used=exonerate_options.copy()
  if issubclass(profile.__class__, profile_alignment):
    profile_exonerate_options=profile.exonerate_options_dict()
    for option_name in profile_exonerate_options:     
      if profile_exonerate_options[option_name]=='DEFAULT':      del exonerate_options_used[option_name]
      else:       exonerate_options_used[option_name]=profile_exonerate_options[option_name]
  else: # allowing use of simple alignment class as profile_ali in input
    p=profile_alignment()
    for k in profile.__dict__:   p.__dict__[k]=profile.__dict__[k]
    profile=p
    try: profile.name=str(profile.name)
    except AttributeError: profile.name='NoName'
    if not profile.queries: profile.queries=range(profile.nseq())

  for title in profile.titles():
    seq= profile.fix_multiple_secs(  profile.seq_of(title, sec_columns='*')  )    
    profile.set_sequence(title, seq)

  #initializing
  done_left = 1 ; done_right = 1 ; done_change_query = 0 ; iterations = 1 ; cyclic_not_worth_doing=0; empty_exonerate=0
  chromosome_file = fastafetch(split_folder, seed.chromosome, target_file)
  chromosome_length= chromosome_lengths_hash[seed.chromosome]
  if issubclass(seed.__class__, blasthit):
    current_alignment = alignment()
    complete_query_title=profile.fill_title(seed.query.chromosome)
    current_alignment.add(complete_query_title, replace_chars(seed.alignment.seq_of('q'), 'U', '*'))
    current_alignment.add('target', seed.alignment.seq_of('t'))
    current_alignment.fill_sides(profile, inplace=True, wild_chars='UX*')
    profile_plus_prediction= profile.transfer_alignment(current_alignment)
    query_name = choose_query_from_profile_plus_prediction(profile_plus_prediction, profile, target_name='target')    
  else:
    query_name= profile.queries_titles()[0]
  
  # computing filenames, all in temporary folder
  query_filename = temp_folder+'cyclic_exonerate_query.fa'
  target_filename =temp_folder+'cyclic_exonerate_target.fa' #unless no cyclic procedure is used
  tempout_filename=temp_folder+'cyclic_exonerate_output'
  
  #current_range is the gene object which is used to generate the target fasta file given as input to exonerate
  current_range=seed.boundaries_gene().extend(left=extension, right=extension)
  current_range.check_boundaries( chromosome_length )
  if float(chromosome_length)/ current_range.span() <2:   ## in this case (chromosome not much bigger than current_range) no cyclic extensions are performed: the whole chromosome is used
    done_left = 1 ; done_right = 1
    target_filename=chromosome_file
    cyclic_not_worth_doing=1
  if current_range.boundaries()[0]==1:                     done_left=1
  if current_range.boundaries()[1]==chromosome_length:     done_right=1
  e=None

  ## main while loop: it will exit from here only once it's done
  while not ( done_left and done_right and done_change_query ) and iterations<max_iterations and not empty_exonerate:
    #preparing query : it is query_name
    query_full_sequence=nogap(profile.seq_of(query_name))
    write_to_file(">"+query_name+'\n'+query_full_sequence, query_filename)
    #preparing target
    if not ( cyclic_not_worth_doing ):      current_range.fasta_sequence(to_file=target_filename, chromosome_file=chromosome_file, title='fasta_title')

    ### run exonerate!!! 
    exonerate(query_filename, target_filename, outfile=tempout_filename, exhaustive=False, exonerate_options=exonerate_options_used, dont_parse=True, mode=mode)

    e=superexoneratehit(); error_message=e.load(tempout_filename, seed=seed, query_full_sequence=query_full_sequence, merge_multiple=merge_multiple)
    if not e:      empty_exonerate=True; e.error_message=error_message
    else:
      #computing current_alignment and profile_plus_prediction -> this will be used to compute the most similar query in the profile to the target. This will be then given as query to exonerate
      current_alignment=alignment()
      complete_query_title=profile.fill_title(e.query.chromosome)
      current_alignment.add(complete_query_title, replace_chars(e.alignment.seq_of('q'), 'U', '*'))
      current_alignment.add('target', e.alignment.seq_of('t'))
      current_alignment.fill_sides(profile, inplace=True, wild_chars='UX*')
      profile_plus_prediction=  profile.transfer_alignment(current_alignment) 

      new_query_name = choose_query_from_profile_plus_prediction(profile_plus_prediction, profile, target_name='target') #choosing query
      
      if new_query_name!= query_name and not iterations==max_iterations-1:        query_name=new_query_name #the iterations condition is to be sure to have annotated the corrected query name even when we stop because of the iteration number
      else:                                  done_change_query=True
     
      ## checking if we filled the boundaries on left/right
      if current_alignment.seq_of('target')[0]!='-':      done_left=1
      if current_alignment.seq_of('target')[-1]!='-':     done_right=1
      if current_range.boundaries()[0]==1:                     done_left=1
      if current_range.boundaries()[1]==chromosome_length:     done_right=1

    iterations+=1
  ##
  #  end cycling
  
  ### if outfile is provided, now we move our output in the temporary folder to that desired destination. If it is not, anyway we have all the information loaded in the returned object 
  header='#blast_seed:'+seed.header(no_species=True)+'\n#profile_alignment: '+profile.name+' ; query_name: '+query_name+' ; query_full_seq: '+query_full_sequence
  if outfile:
    write_to_file(header, outfile)
    bbash('cat '+tempout_filename+' >> '+outfile)

  return e

def genewise(profile_ali, target_file, outfile='', seed='', extension=15000, genewise_options={}):
  """ This function runs genewise on the target using the specified profile and using as seed (determining the positions of search) a gene object, whose boundaries are used, after being extended by the extension parameter.
  If outfile is not indicated, it is put on a temporary file.   A genewise object is returned."""

  global chromosome_lengths
  
  #replacing the input profile with a copy which change in all sequences all aminoacids in sec columns with U (in exonerate, it is *)
  profile=profile_ali.copy()
  for title in profile.titles():
    seq= profile.fix_multiple_secs( profile.seq_of(title, sec_columns='U') )
    profile.set_sequence(title, seq)

  #determining genewise options: the ones specified for the profile override the ones specified as arguments of the genewise function
  genewise_options_used=genewise_options.copy()
  profile_genewise_options=profile.genewise_options_dict()
  for option_name in profile_genewise_options:     
    if profile_genewise_options[option_name]=='DEFAULT':      del genewise_options_used[option_name]
    else:      
      genewise_options_used[option_name]=profile_genewise_options[option_name]

  chromosome_file = fastafetch(split_folder, seed.chromosome, target_file)
  chromosome_length= chromosome_lengths[seed.chromosome]
  
  ### choosing the query. procedure depends on the type of the seed provided : if it is exonerate, just pick same query. If it is blast, map the alignment to the profile and choose most similar query
  seed_type=''
  if issubclass(seed.__class__, blasthit):
    current_alignment = alignment()
    complete_query_title=profile.fill_title(seed.query.chromosome)
    current_alignment.add(complete_query_title, seed.alignment.seq_of('q'))
    current_alignment.add('target', seed.alignment.seq_of('t'))
    current_alignment.fill_sides(profile, inplace=True, wild_chars='UX*') # correcting for modified U (the blast query may not have U but it is put in all sec columns, so than the sequence read may be different than this)  
    profile_plus_prediction= profile.transfer_alignment(current_alignment)
    query_name = choose_query_from_profile_plus_prediction(profile_plus_prediction, profile, target_name='target')    
    seed_type='blast_'
  elif issubclass(seed.__class__, exoneratehit):
    query_name= profile.fill_title(seed.query.chromosome)
    seed_type='exonerate_'
  else:
    query_name= profile.queries()[0]
  
  query_filename = temp_folder+'genewise_query.fa'
  target_filename =temp_folder+'genewise_target.fa' 
  tempout_filename=temp_folder+'genewise_output'
  error_file=temp_folder+'genewise_errors_unclean'
  error_file_cleaned=temp_folder+'genewise_errors'
  
  # preparing query  
  query_full_sequence=nogap(profile.seq_of(query_name))
  write_to_file(">"+query_name+'\n'+query_full_sequence, query_filename)
  # preparing target
  current_range=seed.boundaries_gene().extend(left=extension, right=extension)
  current_range.check_boundaries( chromosome_length )
  current_range.fasta_sequence(to_file=target_filename, chromosome_file=chromosome_file, title='gw_target') #negative strand is already fastarevcomp'ed

  header='#'+seed_type+'seeded; Target_range:'+seed.boundaries_gene().header(no_species=True)+'\n#profile_alignment: '+profile.name+' ; query_name: '+query_name+' ; query_full_seq: '+query_full_sequence+' ; target_for_this_genewise_run: '+current_range.fasta_title()
  write_to_file(header, tempout_filename)
  report=''
  
  b=bash('which genewise')
  if b[0]: raise notracebackException, "ERROR genewise not found! Please install it or skip its execution with -dont_genewise"
  df_genewise= dereference( b[1] )
  add_wisecfg=''
  if df_genewise.endswith('src/bin/genewise'): 
    cfg_dir= df_genewise[:-16]+'wisecfg'
    add_wisecfg= 'export WISECONFIGDIR="'+cfg_dir+'"; '        
  # b holds the exit status of genewise (in pos 0)
  b=bash(add_wisecfg+' '+df_genewise +' -pretty -sum -gff '+join(['-'+str(k)+' '+str(genewise_options_used[k]) for k in genewise_options_used], ' ')+' '+query_filename +' "'+target_filename+'"  >> '+tempout_filename+' 2> '+   error_file)
  bash("sed 's/Cells done \[[ ]*-*[0-9][0-9]*%\]//g' "+error_file+" |  sed 's/Cells [0-9]*%\\x08\[ [0-9]*\]//g' | sed  's/\[[0-9][0-9]*% done\]Before mid-j [0-9][0-9]* Cells done [0-9][0-9]*%\\x08//g' | sed 's/\[[0-9][0-9]*% done\]After[ ]* mid-j[ ]* [0-9][0-9]* Cells done [0-9][0-9]*%\\x08//g' |   sed 's/Explicit read off.*\[[0-9][0-9]*,[0-9][0-9]*\]\[[0-9][0-9]*,[0-9][0-9]*\]//g' | sed 's/\\x[0-9]*//g' | sed 's/ [ ]*/ /g'   > "+error_file_cleaned)
  if b[0]: # genewise exited with exit status != 0. Some error occured
    if int(bash('grep -c "Segmentation fault" '+error_file_cleaned)[1])>0:      report="Segmentation fault" # known problem 
    else:                                                                       report= join(open(error_file_cleaned, 'r').readlines(), '\n')
    if 'Could not read a GeneFrequency file in human.gf' in report:
      raise notracebackException, "ERROR genewise is not installed properly. Please visit  http://big.crg.cat/news/20110616/installing_programs_and_modules_needed_by_selenoprofiles to fix it. This is the error message: "+report+'\n'
    g=genewisehit() 
    g.error_message=report
  else:
    ## preparing genewisehit object, loading it from the file and returning it. Some information are read from the header that we just put in it.
    g=genewisehit() 
    g.load(tempout_filename)
  if outfile: bbash('mv '+tempout_filename+' '+outfile)
  return g

def choose_prediction(candidates, check_filtered=False):  #also called choose_prediction_selenoprofiles so you can link that even if you define a new choose_prediction 
  """ Decision routine. Returns the best object among the provided and a text with the reason why it returned that. 
  The challenge is generally between exonerate and genewise, but blast is treated as equal, apart from the fact that it is not considered when frameshifts are predicted by any program, since blast would be chosen by accident.
  If check_filtered flag is set to True, then the .filtered attribute of p2gs is also checked. The ones with "filtered" and "refiltered" are considered lower level than "kept". "overlapping", "redundant" are not considered since they are not run with this function for construnction.
  """      
  ## putting all non-empty prediction in candidates (list). Then, when we don't like any of these for some reason (checkpoints below) we remove it from the candidate list. When we have a single candidate left, there we have chosen.
  if len(candidates)==0: return empty_p2g(), 'None available'
  if len(candidates)==1: return candidates[0], 'only one available'

  #this is used with database remove redundancy. the candidates with .filtered attribute which is "filtered" or  "refiltered" are not chosen, unless all candidates are like this
  if check_filtered:
    if not  all( [ c.filtered in ["filtered", "refiltered"]  for c in candidates]) :
      for candidate_i in range(len(candidates)-1, -1, -1):
        if candidates[candidate_i].filtered in ["filtered", "refiltered"]: candidates.pop(candidate_i)    

  # frameshifts: if any is predicted, we don't want blast prediction. but if we have frameshifts only in a program (!=blast), then we want the other one
  non_blast_candidates=[]
  for c in candidates: 
    if c.prediction_program()!='blast': non_blast_candidates.append(c)
  any_frameshifts_predicted=any([ c.frameshifts() for c in non_blast_candidates ]) #are there any frameshifts predicted in any candidate? -> False or True
  if any_frameshifts_predicted :
    candidates=non_blast_candidates
    if len(candidates)==1: return candidates[0],  "frameshifts exclude blast"
    min_frameshifts_predicted=min([ c.frameshifts() for c in candidates])
    if min_frameshifts_predicted==0:
      for candidate_i in range(len(candidates)-1, -1, -1):
        candidate=candidates[candidate_i]
        if candidate.frameshifts() > min_frameshifts_predicted:            candidates.pop(candidate_i)
      if len(candidates)==1: return candidates[0],  'other has frameshift'

  # aligned SecTGAs
  max_aligned_secTGAs=max([ c.aligned_secTGAs() for c in candidates])
  for candidate_i in range(len(candidates)-1, -1, -1):
    candidate=candidates[candidate_i]
    if candidate.aligned_secTGAs() < max_aligned_secTGAs:            candidates.pop(candidate_i)
  if len(candidates)==1: return candidates[0], 'SecTGA aligned'

  # stop codons: is there a prediction without any?
  min_stop_codons=min( [c.n_stop_codons() for c in candidates] )
  if min_stop_codons==0: 
    for candidate_i in range(len(candidates)-1, -1, -1):
      candidate=candidates[candidate_i]
      if candidate.n_stop_codons() > min_stop_codons:            candidates.pop(candidate_i)
  if len(candidates)==1: return candidates[0], 'other has in-frame stops'
  
  # protein length
  max_prot_length= max( [len(replace_chars(c.protein(), 'x', '')) for c in candidates] )
  for candidate_i in range(len(candidates)-1, -1, -1):
    candidate=candidates[candidate_i]
    if len(replace_chars(candidate.protein(), 'x', '')) < max_prot_length:     candidates.pop(candidate_i)
  if len(candidates)==1: return candidates[0], 'longest CDS predicted'
  
  if all( [candidates[-1].prediction_program() == c.prediction_program() for c in candidates]): #all candidate predicted by same program .. this is basically for blast when merging results from different clusters of the same profile
    candidates.sort(   key=lambda x:x.id  )
 
  return candidates[-1], 'highest priority program' #if we arrived this far, predictions left in candidates have similar features, and same protein length

choose_prediction_selenoprofiles=choose_prediction

def assign_label_selenoprotein_families(p2g_hit):
  """ This function is used to assign a label for hits predicted for selenoprotein families (containing at least a U in the profile)"""
  if not p2g_hit: raise Exception, "ERROR can't assign label to empty p2ghit!"
  if not p2g_hit.frameshifts() and p2g_hit.stop_codons() and all([stop_codon=='UGA' for stop_pos, stop_codon in p2g_hit.stop_codons()]): return 'uga_containing'  # stop_pos in 1 based
  if p2g_hit.n_stop_codons() or p2g_hit.frameshifts(): return 'pseudo'
  if 'U' in p2g_hit.protein(): return "selenocysteine"
  else:
    u_pos_in_ali=p2g_hit.alignment.all_positions_of('U')
    if not u_pos_in_ali :  return 'unaligned'
    aligned_aa_single_lett= p2g_hit.alignment.seq_of('t')[u_pos_in_ali[0]]  # checking the first of sec positions to determine the label of this homologue
    if aligned_aa_single_lett in '-x': return 'gapped'
    else:
      if aligned_aa_single_lett=='X': return 'unknown'
#      elif aligned_aa_single_lett=='*': return 'whattha'
      elif not one_letter_to_three_letter_aa_diz.has_key(aligned_aa_single_lett): 
        raise notracebackException, "assign_label_selenoprotein_families ERROR unknown amino acid aligned to Sec position: "+aligned_aa_single_lett
      return one_letter_to_three_letter_aa_diz[aligned_aa_single_lett] # e.g. alanine


def assign_label_non_selenoprotein_families(p2g_hit):
  """ This function is used to assign a label for hits predicted for non selenoprotein families (not containing any Us in the profile)"""
  if not p2g_hit: raise Exception, "ERROR can't assign label to empty p2ghit!"
  if p2g_hit.stop_codons() or p2g_hit.frameshifts(): return 'pseudo'
  return 'homologue'

############################################# MERGE FUNCTIONS ##################################################################
def merge_p2g_hits_by_colinearity(p2g_list, inplace=True, already_sorted=False, max_overlap_in_target=0, max_overlap_in_query=6, max_distance_in_target=10000, max_distance_in_query=None, post_function=None,  sequence_collection=None, get_sequence=None, no_seq=False): 
  """ Input can be any iterable of gene objects having a query attribute which is a gene itself. examples: [super]blasthit, exonerate, genewise classes.
  
  The post_function, if defined,  is used to add attributes from the colinear objects to the new one. The function must take 3 arguments (gene objects or subclasses): the upstream gene object, the downstream one, the resulting one, and return a single gene object which will overwrite the resulting object.
  By default the .alignment is modified to include the sequence of the joined p2g hits.
    When portions of the query are missing, one or more "x" are inserted to fill the gap.  
    If you provide a sequence_collection, the method seq_of will be used to obtain the missing sequence in the query. The target is filled just with gaps.
    For the same purpose, you may use instead any get_sequence function which given a fasta title, returns the sequence.
  To ignore the .alignment attribute (when loading from tab blast files), use no_seq=True    
  """  
  if inplace:    worklist=p2g_list
  else:          worklist=list(p2g_list)
  if not already_sorted:        worklist.sort(cmp=order_genes_for_chr_strand_pos)
  index=1 ; 
  full_seq_query=None
  while index<len(worklist):
    this=worklist[index-1]  ;     next=worklist[index]  
#    print [this], [next], index
    if this.chromosome == next.chromosome and this.strand == next.strand:
      if this.strand=='-':         this, next= next, this
      distance_in_target= this.is_upstream_of(next, max_overlap=max_overlap_in_target)
      distance_in_query = this.query.is_upstream_of(next.query, max_overlap=max_overlap_in_query)
      if (distance_in_target and (max_distance_in_target==None or distance_in_target <= max_distance_in_target)) and (distance_in_query and (max_distance_in_query==None or distance_in_query <= max_distance_in_query)):
        
        new_one=this.union_with(next);                    new_one.id=this.id
        new_one.query=this.query.union_with(next.query);  new_one.query.id=this.query.id
        if not no_seq:
          if distance_in_query>0: 
            seq_target= this.alignment.seq_of('t')+(distance_in_query-1)*'-'+next.alignment.seq_of('t')

            if  sequence_collection is None and get_sequence is None:               seq_query=  this.alignment.seq_of('q')+(distance_in_query-1)*'x'+next.alignment.seq_of('q')
            else:
              if full_seq_query is None:
                try: 
                  if not sequence_collection is None:                get_sequence= lambda x:nogap(sequence_collection.seq_of(x))
                  full_seq_query= get_sequence( this.query.chromosome )
                except: printerr('merge_p2g_hits_by_colinearity ERROR can\'t obtain the sequence of: '+this.query.chromosome ); raise
              seq_query=  this.alignment.seq_of('q')+   full_seq_query[this.query.boundaries()[1]:next.query.boundaries()[0]-1]  +next.alignment.seq_of('q')
          else:

            distance_in_query_corrected_with_gaps_in_q =  distance_in_query; p=0; gaps_encountered_query=0; gaps_encountered_target=0
            ali_positions_to_remove=0; 

            for position in range( next.alignment.length() ):
              if next.alignment.seq_of('q')[position] == '-':               gaps_encountered_query+=1
              if position-gaps_encountered_query==-distance_in_query: # we find the right position
                first_position_of_ali_to_keep=position #0 based 
                break
              if next.alignment.seq_of('t')[position] == '-':               gaps_encountered_target+=1

            seq_query=  this.alignment.seq_of('q')+next.alignment.seq_of('q')[first_position_of_ali_to_keep:]
            seq_target= this.alignment.seq_of('t')+next.alignment.seq_of('t')[first_position_of_ali_to_keep:]

            if distance_in_query<0: # we corrected the alignment sequence. We have to correct the gff (the numbers inside .exons ) as well
              to_remove= next.subseq( 1, (first_position_of_ali_to_keep-gaps_encountered_target)*3)
              new_one=new_one.subtracted_of(to_remove)

          new_one.alignment=alignment(); new_one.alignment.add('q', seq_query); new_one.alignment.add('t', seq_target)

        #replacing current ones in list with new one.
        if post_function:           new_one=post_function(this, next, new_one)
        worklist[index-1:index+1]=[new_one]
        index-=1
    index+=1
  return worklist    

def merge_p2g_hits_by_colinearity_post_function_blast_hits(upstream_g, downstream_g, result_g):
  """ Post function for function merge_p2g_hits_by_colinearity when used on blast hits. Returns a superblasthit.  Beware: this means that the merge_p2g_hits_by_colinearity function returns a mixture of blasthits and superblasthits. """
  #transforming blast hit in a superblasthit class. 
  out_g=superblasthit()
  for i in result_g.__dict__:      out_g.__dict__[i]=result_g.__dict__[i]
  out_g.query=result_g.query.boundaries_gene()
  out_g.evalue= min (upstream_g.evalue, downstream_g.evalue)
  out_g.merged=[]
  for i in upstream_g, downstream_g:
    if i.__class__.__name__ == 'superblasthit':
      for blast_h in i.merged:        out_g.merged.append(blast_h)
    else:                           out_g.merged.append(i)
  out_g.reset_derived_data() #the cds of this must be computed again. That will happen automatically next time .cds() is called
  out_g.join_adjacent_exons()
  return out_g

def merge_p2g_hits_by_colinearity_post_function_exonerate_hits(upstream_g, downstream_g, result_g):
  """ Post function for function merge_p2g_hits_by_colinearity when used on exonerate hits. Returns a superexoneratehit.  Beware: this means that the merge_p2g_hits_by_colinearity function returns a mixture of exoneratehit and superexoneratehit. """
  out_g=superexoneratehit()
  for i in result_g.__dict__:      out_g.__dict__[i]=result_g.__dict__[i]
  out_g.query=result_g.query.boundaries_gene()
  out_g.score= upstream_g.score+ downstream_g.score
  out_g.merged=[]

  for i in upstream_g, downstream_g:
    if i.__class__.__name__ == 'superexoneratehit':
      for exonerate_h in i.merged:        
        #write('appending hits in s.e. '+str(i.id)+' : '+str(exonerate_h.id)+ '  to out gene: '+str(out_g.id ), 1)
        out_g.merged.append(exonerate_h)
    else:                           
      #write('appending hit e.   '+str(i.id)+'   to out gene: '+str(out_g.id ), 1)
      out_g.merged.append(i)  
  out_g.reset_derived_data() #the cds of this must be computed again. That will happen automatically next time .cds() is called

  return out_g

############################################# OTHER FUNCTIONS ##################################################################

def choose_query_from_profile_plus_prediction(profile_plus_prediction, profile_ali, target_name='target'):
  """ used with the profile and the profile with the prediction, returns the most similar query tagged protein"""
  possible_tagged_queries=profile_ali.queries_titles()
  if not possible_tagged_queries: raise notracebackException, "cyclic_exonerate ERROR no sequences are tagged as queries in the profile " +profile_ali.name
  ordered_queries=profile_plus_prediction.order_by_similarity_with(target_name)  
  for q_title in ordered_queries:
    if q_title in possible_tagged_queries:        return q_title

def choose_among_overlapping_p2gs_intrafamily(p2g_hit_A, p2g_hit_B):
  """ this function is used when removing intra-profile redundancy to decide which prediction is kept when two overlap. This is basically a link to the choose_prediction function which is used again here """
  if p2g_hit_A.positions_summary() == p2g_hit_B.positions_summary(): return cmp(p2g_hit_A.id, p2g_hit_B.id) ## taking care of ties
  g=choose_prediction_selenoprofiles([p2g_hit_A, p2g_hit_B])[0]
  if g==p2g_hit_A: return -1
  else:            return +1

def choose_among_overlapping_p2gs_interfamily(p2g_hit_A, p2g_hit_B):
  """ function called to choose which prediction to choose during the inter-profile redundancy check step. it also check the .filtered attribute (counts only if -full_db is active) """
  #filtered attribute:                                                                                                                                                                                                           
  if     p2g_hit_A.filtered in ['filtered', 'refiltered'] and not p2g_hit_B.filtered in ['filtered', 'refiltered']: return +1
  elif   not p2g_hit_A.filtered in ['filtered', 'refiltered'] and p2g_hit_B.filtered in ['filtered', 'refiltered']: return -1
  #seq id                                                                                                                                                                                                                        
  for i in [ p2g_hit_A, p2g_hit_B ]: 
    if not i.awsi_data.has_key(True): raise Exception, "ERROR the database does not have stored the sequence identity value for prediction: "+i.output_id()+" Please rerun all profiles with option -D"
  awsi_A=p2g_hit_A.awsi_data[True]
  awsi_B=p2g_hit_B.awsi_data[True]
  if    awsi_A > awsi_B :     return -1
  elif  awsi_A < awsi_B :     return +1
  return cmp(p2g_hit_A.id, p2g_hit_B.id) #taking care of ties

def description_format(line):
  """ This function take as input a description line for an exonerate/genewise hit, obtained using the header function of the gene class (option no species and no id) and turns it into a nicer looking string so that all those passing through this function will be graphically aligned.
  example: sps_temp.25.genewise_tbs -> chromosome:10 strand:- positions:13426764-13426956,13420711-13420814,13418249-13418356,13415823-13415977,...,13401151-13401362
  """
  global max_chars_per_column
  out=''
  out+=line.split()[0].ljust(max_chars_per_column['id'])+' ' #sps_temp.25.genewise
  out+=line.split()[1]+' ' #->
  if 'chromosome:' in line.split()[2]:
    out+=line.split()[2][11:].ljust(max_chars_per_column['chromosome'])+' ' #chromosome:10
    out+=line.split()[3][7:].ljust(max_chars_per_column['strand'])+'   '     #strand:-
    out+=line.split()[4][10:]                                            #positions:13426764-13426956,13420711-13420814,13418249-13418356,13415823-13415977,...,13401151-13401362
  else: out+= join(line.split()[2:], ' ')
  return out

def load_p2g_from_db_result( db_result, profile=None ):
    """ Note: uses global variable results_db
    function to transform the object coming out from the db to a p2ghit of the right class. 
    example of db_result: 
    db_id,  profile, program, target_header, query_header, chromosome, alignment_query, alignment_target, label, state, awsi
    
    the profile name is ignored cause we need something to go from the profile name to the profile instance.  If provided as argument of this function, it's added.
    If profile is not provided as argument, a .profile_name attribute is added to the object. This is useful cause sometimes you deal with results for which you don't want to load the profile. In this case, some methods (such .output_id() ) can still be called, since they look for this .profile_name attribute instead of profile.name
    
    When the results.state attribute is stored in the database with an index ( e.g.   redundant_89 , overlapping_SelK.7 ), the global variable results_db is quiered to retrieve the target gene and this is loaded (through this same function) and put in the .overlapping ; the recursion will only go the a single level, since "redundant" and "overlapping" results are not checked for redundancy so cannot be pointed by .overlapping attributes
    """
    db_id= db_result[0]
    profile_name= db_result[1]
    program= db_result[2]
    target_header= db_result[3]
    query_header= db_result[4]
    chromosome=   db_result[5]
    alignment_query=  db_result[6]
    alignment_target= db_result[7]
    label= db_result[8]
    state= db_result[9] #see below
    awsi= db_result[10]

    if program == 'blast':                p=blasthit()
    elif program == 'exonerate':          p=exoneratehit()
    elif program == 'genewise':           p=genewisehit()
    else: raise Exception, "load_p2g_from_db_result ERROR program not recognized: "+str(program)
    try:     p.load_from_header(target_header)
    except:  raise Exception, "ERROR malformed header stored in database! profile: "+profile_name+' ; db_id: '+str(db_id)+' ; target_header: '+target_header
    p.chromosome= chromosome
    p.query.load_from_header(query_header)
    p.query.id= p.id+'_query'
    p.query.strand=  '+'

    p.alignment.add('q', alignment_query )
    p.alignment.add('t', alignment_target )
    p.label=label
    if profile: p.profile=profile
    else: p.profile_name=profile_name
    p.species=target_species
    p.target= target_file
    if '_' in state:  #either redundant or overlapping. let's add the .overlapping attribute to the object before returning it
      if state.startswith('redundant'):
        other_profile_name=profile_name
        index_id = state.split('_')[1]
        state=state.split('_')[0]

      elif 'overlapping' in  state:
        other_profile_name=join(state.split('overlapping_')[1:], '_').split('.')[0]
        index_id = state.split('.')[1]
      if 1: #target: #otherwise the overlapping result cannot be loaded
        res=   results_db.get_result(other_profile_name, index_id)
        if res:        
          p.overlapping = load_p2g_from_db_result( res, profile=profile ) 
          state='overlapping'
        else:
          ###   the result overlapping to the one we're loading was not found in the database. it means some profile was run more than once changing parameters.
          if not state.split('overlapping_')[0]:  #normally the state attribute in the database contains the filtered state followed by "overlapping" and the output id of the result to which is overlapping. e.g. :  refiltered_overlapping_SelK_34  ; anyway this control is to ensure compatibility with versions of selenoprofiles previous to 2.1 , in which the state attribute was not containing the filtered attribute, e.g. overlapping_SelK_34 ; in this case we don't know what was the original filtered attribute. 
            state='unknown'
          else:
            state=state.split('overlapping_')[0]
    p.filtered=state
    p.features.extend( results_db.features_of(    db_id, add_parent=p   ) )
    p.awsi_data[True]=awsi
    return p  

def parse_blast_subject_and_evalues(blast_file):
  """ Simple function to parse a blast output and return the list of [subject name, evalue].
  This function was built to be used with tag_blast, since the normal parse_blast(er) methods do not keep the full subject names.
  """
  if not is_valid_blast_output(blast_file): raise Exception, "parse_blast_subject_and_evalues ERROR not a valid blast output file:"+str(blast_file)
  outlist=[]
  blast_file_h=open(blast_file, 'r')  
  bline=blast_file_h.readline()
  while bline:
    while bline and not bline.startswith('>'):     bline=blast_file_h.readline()
    if bline.startswith('>'):       
      title=bline[1:-1]
      bline=blast_file_h.readline()  
      while bline and not  'Length =' in bline:
        title+=del_white(bline[:-1])
        bline=blast_file_h.readline()  
      while bline and  not 'Expect =' in bline:  bline=blast_file_h.readline()   
      if not bline: raise Exception, "parse_blast_subject_and_evalues ERROR not a valid blast output file: "+str(blast_file)+' ; no evalue found for title:'+title
      evalue= e_v(del_white(bline.split('Expect =')[1].split(',')[0]))
      outlist.append([title, evalue])    
    bline=blast_file_h.readline()
  blast_file_h.close()
  return outlist  

##############################################################################################################################################
############################################ SELENOPROFILES CLASSES ##########################################################################
##############################################################################################################################################
#PROFILE

class profilechangedException(Exception):          
  """ just another exception class"""
class profiledatafilemissingException(Exception):     
  """ just another exception class"""
class clusteringseqidchangedException(Exception):  
  """ just another exception class"""
      
class profile_alignment(alignment):
  """ Subclass of alignment.
  Additional attributes:
  .filename             -> link to the main file loaded.  
  .queries              -> keywords, or list of indices (Starting with 0) of sequences being used as queries of exonerate/genewise -- see make_profile doc for possible syntaxes
  .name                 -> family name (filename without extension, if not specified)
  .blast_options        -> options string concatenated with every command line running psitblastn
  .exonerate_options    -> options string concatenated with every command line running exonerate
  .genewise_options     -> options string concatenated with every command line running genewise
  .tag_blast_options    -> options string concatenated with every command line running blastp for tag_blast (tag_score and go_score methods)
  .blast_filtering      -> python procedure to filter blast hits
  .p2g_filtering        -> python procedure to filter p2g candidates
  .p2g_refiltering      -> python procedure to refilter rp2g candidates
  .tags                 -> list of regexp tags (as strings) used for tag_score. Matches with any of these give a positive score
  .neutral_tags         -> list of regexp tags (as strings) used for tag_score. Matches with any of these avoid putting a negative score (the entry needs to be positive or is skipped)
  .tag_db               -> protein database used for tag_blast (tag_score and go_score methods)
  .clustering_seqid     -> sequence identity threshold for pre-psitblastn clustering
  .max_columns_gaps     -> maximum percent of gaps allowed in a column to be included in the consensus query computed for psitblastn searches
  .max_blast_hits       -> maximum number of blast hits allowed for this profile
  """
  parameters=['blast_filtering', 'p2g_filtering', 'p2g_refiltering',  'blast_options', 'exonerate_options', 'genewise_options', 'tag_db', 'tag_blast_options', 'neutral_tags', 'uniref2go_db', 'clustering_seqid', 'max_column_gaps', 'max_blast_hits', 'tags']
  psitblastn_relative_options=['blast_filtering','blast_options','clustering_seqid', 'max_column_gaps', 'max_blast_hits']

  def __init__(self, filename='', **build_options):
    """ This class can be initialized with an alignment file. If a .config file if found, the profile parameters are read from there, otherwise (or in case any build_option is specified) it is build and the file is created. """
    self.filename='';      
    self.queries=None;
    for category in profile_alignment.parameters:
      exec('self.'+category+'="DEFAULT"')
    self.tags=[]; self.go_terms=[]
    self.sec_pos_data=None; self.name=''
    self.clusters_data=None
    self.blast_queries_data=None
    self.conservation_data=None
    alignment.__init__(self)  #this is initalising also : self.conservation_map_data = None
    if filename:      self.load(filename, **build_options)

  def __setitem__(self, key, value):    self.__dict__[key]=value
  def __getitem__(self, key):           return self.__dict__[key]
  
  def load(self, filename, **build_options):
    """ Loads a alignment into the profile object. If the config file is not present (or if any build_options are provided), the profile is built using defaults when necessary."""
    alignment.__init__(self, filename)
    if not self.check_length():            raise notracebackException, "profile_alignment -> load ERROR can't load "+filename+" : is not an aligned fasta file"
    elif not self:                         raise notracebackException, "profile_alignment -> load ERROR can't load "+filename+" : empty"
    self.filename=filename
    if not is_file(self.filename+'.config') or build_options:   self.make_profile(**build_options)
    else:                                                       self.load_profile()
    if not self.first_words_are_uniq():    raise skipprofileException, "ERROR the first word of titles in the profiles must be unique! This is not true for profile: "+self.name

  def first_words_are_uniq(self):
    """ Returns False if profile names do no respect this rule"""
    first_words_hash={}
    for title in self.titles():
      f_word=title.split()[0]
      if first_words_hash.has_key(f_word): return False
      first_words_hash[f_word]=1
    return True

  def load_profile(self):
    """ This function loads all options from the config filename. If something is not specified in the config, it gets the default value.
    The "queries" attribute may have several forms: 
    - it can be a list of 0-based indexes   e.g.   [0, 1, 3, 4, 5, 7]
    - it can be the word "all"  , which means that all sequences are used as queries (but see below)
    - it can be an expression like best:FLOAT , which means that the first X sequences are used, where X is FLOAT*number of sequences in the alignment. Remember that the sequences are ordered by "completeness".  e.g.   best:0.5  -> the first half of the sequences are used
    - it can be an expression like last_query:TITLE, which means that all the first sequences up to the one named TITLE are used (TITLE can be even partial: it is searched in each title, the first one matching will be the last query).    e.g.  last_query:gi|73946000|ref|XP_854030.1
    NOTE: in all forms apart from the first one (indexes listing), if the profile is a selenoprotein family, the sequences which lack one or more selenocysteine positions (meaning that they contain a gap in any of the sec_positions) are ignored. To avoid this, you can precede the expression with not_only_u:  ; e.g.   not_only_u:all    e.g.#2  not_only_u:last_query:XP_854030.1
    """
    config_diz={} #keeping track of those changed in this run of load_profile, so then I can set defaults only for the remaining ones.
    for line in open(self.filename+'.config', 'r'):
      if line[0]!='#' and "=" in line:
        key=    del_white( line[:-1].split('=')[0] )
        value=  del_white( join(line[:-1].split('=')[1:], '=') )
        if not ( value[0] in '["\''  or all([c in digits+'.' for c in value])  ) : value='"'+value+'"'        
        exec('config_diz[ key ] = '+ value)
        exec('self.__dict__[ key ] = '+value)
    for this_option in  profile_alignment.parameters:
      if not this_option in config_diz:        self.__dict__[this_option] = 'DEFAULT'
    #if type(self.queries)==str and lower(self.queries)=='all':   self.qui
    if not self.queries_titles():      raise notracebackException, "profile_alignment ERROR the queries list cannot be empty! profile: "+self.name   
    try:   self.load_profile_data()             
    except (profilechangedException, profiledatafilemissingException): 
      printerr('Profile '+self.name+' : profile_data not found or detected change in alignment file! Recomputing profile_data...', 1)
      self.compute_all_profile_data()
      self.save_profile_data()

  def compute_all_profile_data(self, keep_awsi=False):
    """ compute all profile data, for saving it """
    self.build_conservation(keep_awsi=keep_awsi)
    self.conservation_map() #ignoring returned value, just to compute it
    self.clusters()         #ignoring returned value, just to compute it

  def conservation_map(self, **kargs):
    if 'exclude' in kargs: del kargs['exclude']
    return alignment.conservation_map(self, exclude={'BLAST_QUERY_MASTER':1}, **kargs)
      
  def build_conservation(self, silent=False, keep_awsi=False):
    """ compute average and std deviation for the two methods (with coverage or w/o)"""
    self.conservation_data={}    
    if len([ t for t  in self.titles() if t != 'BLAST_QUERY_MASTER'])==1:
      printerr('Single sequence profile!', 1)
      self.conservation_data['with_coverage'] =  [ 1.0, 0.0 ]
      self.conservation_data['without_coverage'] =  [ 1.0, 0.0 ]
      return 
    title2score={};     scores=[]
    for title in self.titles():      
      if title!='BLAST_QUERY_MASTER':      title2score[title]=  self.weighted_seq_identity_of(title, True, exclude={'BLAST_QUERY_MASTER':True})   
    for title in title2score.keys():      scores.append(title2score[title])
    average_score= average(scores)
    std_dev= std_deviation(scores)
    self.conservation_data['with_coverage'] =  [ average_score, std_dev ]
    if keep_awsi:      self.conservation_data['awsi_scores_with_coverage']= title2score
    title2score={};     scores=[]
    for title in self.titles():
      if title!='BLAST_QUERY_MASTER':          title2score[title]=  self.weighted_seq_identity_of(title, False, exclude={'BLAST_QUERY_MASTER':True})
    for title in title2score.keys():      scores.append(title2score[title])
    average_score= average(scores)
    std_dev= std_deviation(scores)
    self.conservation_data['without_coverage'] =  [ average_score, std_dev ]
    if keep_awsi:      self.conservation_data['awsi_scores_without_coverage']= title2score

  def conservation_properties(self, with_coverage=True): 
    """ Returns average awsi score, deviation """
    if self.conservation_data is None: self.build_conservation()
    if with_coverage:  return self.conservation_data['with_coverage']
    else:              return self.conservation_data['without_coverage']

  def has_selenocysteine(self):    return bool(self.sec_pos())

  def make_profile(self, **options):                                                  #replacing this later below
    """As options, you may provide anything that can appear in a profile config file: PARAMETERS 
If any of those are not set, they are not written on the .config file, which is equivalent to set them on DEFAULT
There are 3 ways to define the queries in the profile, ordered by hierarchy:
all             -> All sequences are queries (default)    
[...]           -> You provide directly the list of queries titles, or indices (starting from 0) in the ordered titles list
last_query      -> You provide the last title in the ordered titles list that will be labelled as query. This may be a incomplete title: it just have to be found in a title
query_threshold -> you provide the proportion of the full list of queries in the ordered list that will be labelled

NOTE that in the first, third and fourth cases, those titles for which at least one U position in the alignment contains a gap in the sequence are ignored, unless the option "not_only_U" is defined.  """
    try:
      self.convert_sequences(upper)
      self.remove_empty_columns()      
      self.sort_by_completeness()
      self.display(self.filename) #writing the alignment with the new order of titles, by completeness
      config_file_text=''
      #family name
      if 'name' in options:      self.name=options['name']
      else:                      self.name=join(base_filename(self.filename).split('.')[:-1], '.')
     
      allowed_chars=lowercase+uppercase+'1234567890'+'_-'
      unvalid_chars=set()
      for char in self.name: 
        if char not in allowed_chars:           unvalid_chars.add(char)
      if unvalid_chars:
        self.name=replace_chars(self.name, list(unvalid_chars), '_')
        printerr('WARNING character'+'s'*int(bool(len(unvalid_chars)>1))+' '+"'"+join(list(unvalid_chars), '')+"'"+' are not allowed in a profile name and are replaced by "_"  --> '+self.name, 1)
      config_file_text+='name = '+(self.name)+'\n'      

      if   'last_query'      in options:       last_query=str(options['last_query'])
      elif 'query_threshold' in options:       query_threshold=float(options['query_threshold'])
      else:                               query_threshold=1.0
      #queries tagging
      if 'queries' in options.keys():
        self.queries=options['queries']
        # the list of queries was provided as indices 
      else:
        query_expression=''
        if 'not_only_U' in options and options['not_only_U']:   query_expression+='not_only_u:'
      
        if 'last_query' in options:    
          if all([  not str(options['last_query']) in title    for title in self.titles()]): raise notracebackException, "make_profile ERROR cannot find last query specified: "+str(options['last_query'])
          query_expression+='last_query:'+str(options['last_query'])                
        elif 'query_threshold' in options:
          query_expression+='best:'+str(options['query_threshold'])
        else: 
          query_expression+='all'
        self.queries=query_expression
        
      config_file_text+='queries = '+str(self.queries)+'\n'
      if not self.queries_titles():      raise notracebackException, "profile_alignment ERROR the queries list cannot be empty! profile: "+self.name
      #other options... those not specified are not written in the config file, and will be set to defaults. NB these options can also be set to keywords which are set in the selenoprofiles config file.
      for this_option in profile_alignment.parameters:
        if not this_option.endswith('_options') and this_option in options:
          self[this_option]=options[this_option]
          config_file_text+= this_option+' = '+str(options[this_option])+'\n'

      # options for running blast, exonerate, genewise. using the SELENO keyword for selenoprotein families, DEFAULT values for the others
      for this_option in ['blast_options', 'exonerate_options', 'genewise_options']:
        value=None
        if this_option in options:          value=options[this_option]          
        elif self.has_selenocysteine():              value='SELENO'        
        if value:
          self[this_option]=value
          config_file_text+= this_option+' = '+str(value)+'\n'
      write_to_file(config_file_text, self.filename+'.config')
      self.add( 'BLAST_QUERY_MASTER', self.blast_queries().seq_of('BLAST_QUERY_MASTER')  )
      keep_awsi=False
      if 'keep_awsi' in options and options['keep_awsi']: keep_awsi=True
      self.compute_all_profile_data(keep_awsi)
      self.save_profile_data()
      if self.nseq()<3 and not 'silent' in options: printerr('WARNING building profile: '+self.name+' has very few sequences! If default filtering is active (awsi_filter), the distribution of AWSI scores will not be taken into account for this profile.', 1)      
      
    except Exception, e:
      printerr('profile_alignment->make_profile ERROR building profile: '+self.name )
      raise
  make_profile.__doc__ =   replace( make_profile.__doc__ , 'PARAMETERS', join(parameters, ', ') )

  def save_profile_data(self):
    """ dump a .profile_data file to avoid repeating the computation. """
    fileout=self.filename+'.profile_data'
    try:
      fh=open(fileout, 'w')
      print >> fh, "# md5sum id for alignment: "+self.md5sum_id()
      print >> fh, "# Clustering seq id value: "+str(round( self.clustering_seqid_value(), 3 ))
      print >> fh, "# Number of clusters: "+str(self.n_clusters())
      print >> fh, "# Blast queries alignment:\n"+str(self.blast_queries())
      for cluster_index in range(self.n_clusters()):
        print >> fh, "# Cluster alignment n."+str(cluster_index+1)+":\n"+str(self.clusters_data[cluster_index].display(return_it=True))
      print >> fh, "# Average weighted sequence identity with coverage: "+str(self.conservation_properties(with_coverage=True)[0])+' +- '+str(self.conservation_properties(with_coverage=True)[1])
      print >> fh, "# Average weighted sequence identity without coverage: "+str(self.conservation_properties(with_coverage=False)[0])+' +- '+str(self.conservation_properties(with_coverage=False)[1])
      print >> fh, "# Conservation map: ["
      for index, dictionary in enumerate(self.conservation_map()):
        print >> fh, str(dictionary)+',  #'+str(index+1) #each line has a commented position number, one-based, to be able to read it
      print >> fh, "]"
      fh.close()
    except IOError:    printerr('WARNING trying to save profile data, couldn\'t write '+fileout+'! The profile data will be recomputed when necessary! ', 1)

  def load_profile_data(self):
    """ It loads all saved data for this profile. It fills the attributes: .clusters_data, .blast_queries_data, .conservation_data, self.conservation_map_data
    It will recompute data if the profile appears to have changed (or the clustering_seqid value changed -> recompute clusters).
    """
    profile_data_file= self.filename+'.profile_data'
    try:
      if not is_file(profile_data_file): raise profiledatafilemissingException
      fh=open(profile_data_file, 'r' )
      md5sum_line= fh.readline()
      if md5sum_line.split()[-1]!=self.md5sum_id():                   raise profilechangedException
      seqid_value= round(float(fh.readline().split()[-1]), 3)

      if seqid_value != round(self.clustering_seqid_value(), 3):      
        printerr('WARNING The value of clustering_seqid changed, recomputing the clusters for profile: '+self.name, 1)
        self.clusters()
        line=fh.readline() 
        while not line.startswith('# Average weighted sequence identity with coverage:'):               line=fh.readline() #setting filehandler in the right place
      else:
        useless_line=fh.readline(); useless_line=fh.readline();
        #### loading clusters data 
        #now parsing blast queries alignment
        line=fh.readline(); title=''; self.blast_queries_data = alignment()
        while not line.startswith('# Cluster alignment'):
          if line[0]=='>': 
            if title:  self.blast_queries_data.add(title, seq)
            title=line.rstrip()[1:]
          else: seq=line.rstrip()
          line=fh.readline()
        self.blast_queries_data.add(title, seq) # adding last seq (blast query master)
        #now parsing clusters alignment
        cluster_index=0; self.clusters_data=[profile_alignment()]        
        while not line.startswith('# Average weighted sequence identity with coverage:'):
          #now line is: # Cluster alignment n.X:
          line=fh.readline() #now in title of first sequence
          title=''; seq=''; 
          while line and not line.startswith('#') :
            if line[0]=='>': 
              if title:  self.clusters_data[cluster_index].add(title, seq)
              title=line.rstrip()[1:]
            else: seq=line.rstrip()
            line=fh.readline()
          self.clusters_data[cluster_index].add(title, seq) #adding last seq of this cluster
          if line.startswith('# Cluster alignment'):
            cluster_index+=1
            self.clusters_data.append(profile_alignment())
            title=''; seq=''          
        for ali in self.clusters_data:
          for parameter in profile_alignment.psitblastn_relative_options:          exec('ali.'+parameter+'=self.'+parameter)
          ali.name=self.name+'.'+str(cluster_index+1)

      ##### reading conservation data
      self.conservation_data={}
      self.conservation_data['with_coverage']=   [      float(line.split(': ')[1].split()[0])   ,      float(line.strip().split()[-1])   ]   # average, std_dev
      line=fh.readline()    
      self.conservation_data['without_coverage']=[      float(line.split(': ')[1].split()[0])   ,      float(line.strip().split()[-1])   ]   # average, std_dev
      line=fh.readline() ;       line=fh.readline()      
      #      now in conservation map data.. 
      cons_map=[]
      while not line.startswith(']'):
        dict_this_pos_string=line.split('#')[0]
        exec('dict_this_pos = '+dict_this_pos_string)
        cons_map.append(dict_this_pos)
        line=fh.readline()     
      self.conservation_map_data={'BLAST_QUERY_MASTER':cons_map}
      fh.close()
    except  (profiledatafilemissingException, profilechangedException) as e:  
      if isinstance(e, profilechangedException):
        write( "WARNING the alignment file for the profile "+self.name+" has a different checksum than before! The derived information will be recomputed.", 1)
      self.clusters_data=None ## so it will be recomputed next time it's needed
      self.conservation_map_data=None
      self.conservation_data=None
      self.blast_queries_data=None
      raise 
#    except:
#      raise skipprofileException, "ERROR loading profile_data for profile: "+self.name+' ; try removing the .profile_data file and rerun.'

  def md5sum_id(self):
    """ Returns the md5sum output on the alignmetn file. useful to check if the alignment has changed from a previous run"""
    return bbash('md5sum '+self.filename).split()[0]
        
  def queries_titles(self):
    """ Returns the titles of the sequences tagged as queries. see make_profile help for the accepted expressions """
    if type(self.queries)==str: 
      q_titles=[]
      expression= lower(self.queries).split('not_only_u:')[-1]
      if expression=='all':
        for title in self.titles():
          if  not self.sec_pos() or 'not_only_u:' in  lower(self.queries)  or all([self.seq_of(title)[p]!='-' for p in self.sec_pos()]):     q_titles.append(title)
      elif expression.startswith('best:'):
        threshold= float(  self.queries.split(':')[1].strip()  )
        index_last_seq_to_tag=int(  threshold  *self.nseq())
        for index, title in enumerate( self.titles() ):
          if index <= index_last_seq_to_tag and (not self.sec_pos() or 'not_only_u:' in  lower(self.queries)  or all([self.seq_of(title)[p]!='-' for p in self.sec_pos()]) ):
            q_titles.append(title)
      elif expression.startswith('last_query:'):
        last_query_expression= expression.split('last_query:')[1].strip()
        for title in self.titles():
          if not self.sec_pos() or 'not_only_u:' in  lower(self.queries)  or all([self.seq_of(title)[p]!='-' for p in self.sec_pos()]):     
            q_titles.append(title)        
          if last_query_expression in title: break
      else: raise Exception, "ERROR I can't recognize the expression for the queries attribute of the profile "+self.name+' : '+expression
      return q_titles
    else: 
      try:                    return [self.titles()[i] for i in self.queries]
      except IndexError:      raise skipprofileException, "selenoprofiles ERROR the configuration file for profile "+self.name+" appears incorrect: one or more indexes for the queries are larger than the number of sequences in the profile. Please rebuild your profile (see script selenoprofiles_build_profile.py )."
    
  def sec_pos(self):
    """ Returns the positions in which you have at least a U (0 based)"""
    if self.sec_pos_data is None :      self.sec_pos_data=self.all_positions_of('U')
    return self.sec_pos_data

  def seq_of(self, title, sec_columns=None):
    """ Overriding the standard seq_of method to give the possibility to automatically substitute all aminoacids in the requested seq that align to U columns with a provided character (e.g. * or X).  """
    seq=alignment.seq_of(self, title)
    if sec_columns:
      for sec_pos in self.sec_pos():
        if seq[sec_pos]!='-':
          seq=seq[:sec_pos] + sec_columns + seq[sec_pos+1:]
    return seq

  def fasta(self, only_queries=0, tag_queries=0, to_file=''):
    """ This function return a string with the fasta representation of the profile alignment. Query are indicated with # QUERY only if tag_queries is True."""
    out=''
    for title_index, title in enumerate( self.titles() ): 
      is_query=title_index in self.queries
      if not only_queries or is_query:
        if tag_queries and is_query:          title+=' # QUERY'
        out+='>'+title+'\n'+self.seq_of(title)+'\n'
    if to_file:      write_to_file(out[:-1], to_file) 
    else:            return out

  def clusters(self):
    """ Returns the list of alignments after clustering, already processed adding the blast query """
    if not self.clusters_data:
      self.clusters_data=self.clustering(self.clustering_seqid_value(), reclustering=True, dont_remove_empty_columns=True, outclass=profile_alignment)

      self.blast_queries_data=alignment()
      for cluster_index, ali in enumerate(self.clusters_data):
        ##debug        
        #ali.display(temp_folder+'cluster'+str(cluster_index+1)+'.fa')
        #raw_input(temp_folder+'cluster'+str(cluster_index+1)+'.fa')

        q_title = 'BLAST_QUERY_'+str(cluster_index+1)
        q_seq   = ali.blast_query_seq(sec_char='U') #like .consensus_sequence( self.max_column_gaps_for_blast_query_value() ) 
        #nonetheless  since we'd miss the ones that are not Us in this subcluster but that are U in the main cluster we do also:
        for sec_pos in self.sec_pos():          
          if not self.is_gap_column(sec_pos):          q_seq= q_seq[:sec_pos] + 'U' + q_seq[sec_pos+1:]
        q_seq= self.fix_multiple_secs(q_seq)

        ali.add(q_title, q_seq)
        for parameter in profile_alignment.psitblastn_relative_options:          exec('ali.'+parameter+'=self.'+parameter)
        ali.name          = self.name+'.'+str(cluster_index+1)
        ali.remove_useless_gaps()
        self.blast_queries_data.add(q_title, q_seq)
      self.blast_queries_data.add('BLAST_QUERY_MASTER', self.blast_master_sequence(self.blast_queries_data) )
      #self.save_clusters_data()
    return self.clusters_data

  def clusters_relative_positions(self): 
    """it's a hash of hashes. for each cluster you have a k: title (blast_query something) -> value:  hash like {k: 1-based position -> value: 1-based positions in ali}     -- 1-based! """
    out={}
    for title in self.blast_queries().titles()[:-1]:
      out[title]={}
      for pos in range(len( nogap(self.blast_queries().seq_of(title)) )): #pos is 0 based
        out[title][pos+1]=self.blast_queries().position_in_seq('BLAST_QUERY_MASTER',  self.blast_queries().position_in_ali(title, pos+1) ) 
    return out

  def blast_master_sequence(self, blast_queries_alignment):
    """ Compute the blast master sequence from the various blast queries of the different clusters. This function is used both when the clusters are computed or loaded from a .cluster file. """
    blast_master_seq=''
    for pos in range( blast_queries_alignment.length()  ):
      if blast_queries_alignment.is_gap_column(pos): blast_master_seq+='-'
      else: 
        cons_pos = self.conservation_map()[pos]
        most_cons_aas= sorted( cons_pos.keys(), key=cons_pos.get, reverse=True )
        if 'U' in cons_pos:            blast_master_seq+= 'U'
        elif most_cons_aas[0]!= '-':   blast_master_seq+= most_cons_aas[0]
        else:                          blast_master_seq+= most_cons_aas[1]
    return blast_master_seq

  def n_clusters(self):
    """ Return the number of clusters computed for this profile"""
    return len(self.clusters())

  def pssm(self, fileout=''):
    """ This function builds a pssm from the profile to be used with psiblast. If the file is not specified, it is created in the temp file with the base name of the profile + .pssm.
    In both cases, a file called as the pssm, but with the added extension .ascii is created, containing the pssm in human readable format.
    """
    if fileout:      pssm_filename = fileout
    else:            pssm_filename = temp_folder+self.name+'.pssm.chk'
    pssm_ascii_filename  =  pssm_filename[:-4]+'.ascii'
    ## debug
    #self.display(temp_folder+'temp_doing_pssm.fa')
    #raw_input(temp_folder+'temp_doing_pssm.fa')
                
    blast_ali_file       =  self.blast_format_alignment()
    query_filename_withX =  self.blast_query_file(sec_char='X')
    write_to_file('>noseq\nXXXX', temp_folder+'tiny_db'); bbash('formatdb -i '+temp_folder+'tiny_db -p T -o T')
    b=bash('which blastpgp')
    if b[0]: raise notracebackException, "ERROR blastpgp not found! Please install it"
    blastpgp_bin= dereference( b[1] )
    cmnd=blastpgp_bin+'  -i '+query_filename_withX+' -B '+blast_ali_file+' -j 1 -d  '+temp_folder+'tiny_db '+' -C '+pssm_filename+'  -Q '+pssm_ascii_filename
    b=bash(cmnd)
    if b[0] and "[blastpgp] WARNING: SetUpBlastSearch failed." in b[1]:
      raise notracebackException, "ERROR blastall is not properly installed. Please visit http://big.crg.cat/news/20110616/installing_programs_and_modules_needed_by_selenoprofiles to fix it. This is the error message: "+b[1]+'\n'
    elif b[0]: raise Exception, "unknown ERROR with blastpgp, cmnd:  "+cmnd+' ERROR: '+b[1]
    return pssm_filename

  def pssm_matrix(self, modify_sec=False):
    """ This function returns a python representation of the pssm, as a list of hashes with integers as values and aminacids as keys. An additional key is present, 'self', containing the aminoacid in the query. This can be used as a control. The length of the list is the length of the blast query sequence. """
    if is_file(temp_folder+self.name+'.pssm.ascii'):      ascii_pssm=temp_folder+self.name+'.pssm.ascii'
    else:     ascii_pssm= self.pssm()[:-4]+'.ascii' 

    ## parse pssm
    pssm_hash_list=[] # returned object
    ascii_pssm_fh=open(ascii_pssm, 'r')
    line=ascii_pssm_fh.readline()
    while line =='\n' or line.split()[0]!='1':     
      if not line =='\n' and line.split()[0]=='A':   
        aa_keys=line.split()
        aa_keys=aa_keys[:len(aa_keys)/2  ]  #now in aa_keys we have a list of the aminoacids in the same order as the values will appear in the nextt lines, and not replicated
      line=ascii_pssm_fh.readline()
    #now we're in first value line
    while line!='\n':
      pssm_hash_list.append({})
      splt=line.split()
      pssm_hash_list[-1]['self']=splt[1]
      for i in range( len( aa_keys  )  ):
        aa   = aa_keys[i]
        value= int(splt[i+2])
        pssm_hash_list[-1][aa]=value
      line=ascii_pssm_fh.readline()
    ascii_pssm_fh.close()

    if modify_sec:
      sec_cols=self.sec_pos()
      blast_query_title=self.titles()[-1]
      for pos in range(len(pssm_hash_list)):        
        if pssm_hash_list[pos]['self']=='X' and self.position_in_ali(blast_query_title, pos+1)-1 in sec_cols:
          #sec position
          pssm_hash_list[pos]['UGA']= 8
          pssm_hash_list[pos]['C']=   4
          pssm_hash_list[pos]['R']=   1
          pssm_hash_list[pos]['*']=  -4
          pssm_hash_list[pos]['X']=  -1
        else:
          pssm_hash_list[pos]['*']=  -4
          pssm_hash_list[pos]['UGA']=-4
          pssm_hash_list[pos]['X']=  -1
    else:
      for pos in range(len(pssm_hash_list)):
        pssm_hash_list[pos]['*']=  -4
        pssm_hash_list[pos]['X']=  -1

    return pssm_hash_list

  def blast_queries(self):
    """ Returns a list of the blast queries for this profile, in the form [ [title1, seq1], [title2, seq2] ... ] """
    self.clusters() #just to ensure blast_queries_data has been created and filled.
    return self.blast_queries_data

  def blast_query_title(self ): #, with_suffix=False):
    """ Returns the title of the query used for psiblast"""
    last_title= self.titles()[-1]
    if last_title.startswith('BLAST_QUERY'): return last_title
    return "BLAST_QUERY"
  
  def blast_query_file(self, fileout='', sec_char='U'): #3.1a
    """ This function builds a file that can be used as query of psiblast. Us are replaced with * (or with the argument of sec_char). 
    If fileout is not specified, it is built in the temporary folder and the filename is returned. """
    if not fileout:      fileout=temp_folder+self.name+'.blast_query'
    seq=self.blast_query_seq(sec_char=sec_char)
    seq=nogap(seq)
    write_to_file(">"+self.blast_query_title()+'\n' +   seq  , fileout)
    return fileout

  def blast_query_seq(self, sec_char='U'):   #3.1a
    """ This is used by the previous function to compute the sequence of the blast query. careful: it contains gaps!
    The sequence returned is a consensus position by position of the profile. If a U is found at any column, a U is always put in that position.
    The function utilize a parameter to determine which columns contains to many gaps to be considered in the blast query:  the option defined in the profile (or the default in the main) configuration file max_column_gaps.
    This function is actually never called for the profile alignment: it is called with the clusters subalignments. 
    
    """
    last_title= self.titles()[-1]
    if last_title.startswith('BLAST_QUERY'): 
      seq=self.seq_of(last_title)
      for sec_pos in self.sec_pos():
        seq=seq[:sec_pos] + sec_char + seq[sec_pos+1:]
    else:      seq=self.consensus_sequence( threshold=self.max_column_gaps_for_blast_query_value(), sec_char=sec_char )                
    seq= self.fix_multiple_secs( seq )
    return seq

  def fix_multiple_secs(self, seq):
    """ Bugfix to avoid consecutive U or *, which makes blastall explode"""    
    seq_no_gap=replace(seq, '-', '')
    index=0; index_no_gap=-1
    fixed=False
    while index< len(seq):
        char=seq[index]
        if char != '-': index_no_gap+=1
        if char in '*U' and  index_no_gap >2 and  seq_no_gap[index_no_gap] == seq_no_gap[index_no_gap-1]  and seq_no_gap[index_no_gap-1] == seq_no_gap[index_no_gap-2]:        
          # three consecutive sec chars!
          seq=seq[:index] + 'X' + seq[index+1:]
          seq_no_gap= seq_no_gap[:index_no_gap] + 'X' + seq_no_gap[index_no_gap+1:]
          fixed=True
        index+=1
    #if fixed: print "FIXED: new seq_no_gap = "+seq_no_gap
    return seq     
    
  def blast_format_alignment(self, fileout=''):
    """ This function builds a alignment in the format wanted by blastpgp to build a pssm. It is basically a Tab format, with sequences names that cannot contain any space """
    max_chars=120
    if not fileout:      fileout=temp_folder+self.name+'.blast_format_alignment'
    fileout_h=open(fileout, 'w')
    for title in self.titles():
      title_out=replace_chars(title, ' ', '_')
      title_out=title_out[:min(len(title_out), max_chars)]
      seq=self.seq_of(title)
      print >> fileout_h,  title_out+'\t'+seq
    fileout_h.close()

    return fileout

  def blast_filtering_eval(self):
    """ This returns a function of evaluation of a blast hit. this returns True or False if the blasthit provided passes or not the blast filtering specified for this profile. """
    return self.filtering_eval(category='blast_filtering')
  def p2g_filtering_eval(self):
    """ This returns a function of evaluation of a p2g. this returns True or False if the p2g provided passes or not the filtering specified for this profile. """
    return self.filtering_eval(category='p2g_filtering')
  def p2g_refiltering_eval(self):
    """ This returns a function of evaluation of a p2g. this returns True or False if the p2g provided passes or not the refiltering specified for this profile. """
    return self.filtering_eval(category='p2g_refiltering')

  def filtering_eval(self, category  ):
    """ Utility to save code for the previous 3 functions. category is p2g_refiltering or blast_options etc. ; it checks if the value defined in the profile is a keyword, if it is not, it tries to evaluate it as the second part of a lambda statement, with the object being filtered called x. If not value is defined in the profile, the default of this category is used."""
    if self.__dict__[category]:
      if keywords[category].has_key(self.__dict__[category]):           return keywords[category][self.__dict__[category]]
      else:
        function_text=self.__dict__[category]
        function_text='lambda x:'+function_text
        try:
          exec('f= '+function_text)
          return f
        except:
          printerr("profile_alignment.filtering_eval ERROR can't return "+category+"  function "+function_text, 1)
          raise       
    else: #set default
      if not keywords[category].has_key('DEFAULT'): raise notracebackException, "profile_alignment.filtering_eval ERROR the default value is not defined for " +category
      return keywords[category]['DEFAULT']

  def tag_db_filename(self):
    """ Returns the tag database specified for this profile. This includes translating a keyword, if the database is specified as keyword and not as a complete path. It also checks if the file exists and raises and exception if it doesn't """
    possible_keyword=self.tag_db
    if not self.tag_db: possible_keyword='DEFAULT'
    if keywords['tag_db'].has_key(possible_keyword): possible_keyword= keywords['tag_db'][possible_keyword]
    if not is_file(possible_keyword): raise notracebackException, "tag_db_filename ERROR "+possible_keyword+' not found'
    return possible_keyword

  def uniref2go_db_filename(self):
    """ Returns the uniref2go database for this profile. This includes translating a keyword, if the database is specified as keyword and not as a complete path. It also checks if the file exists and raises and exception if it doesn't """
    possible_keyword=self.uniref2go_db
    if not self.uniref2go_db: possible_keyword='DEFAULT'
    if keywords['uniref2go_db'].has_key(possible_keyword): possible_keyword= keywords['uniref2go_db'][possible_keyword]
    if not is_file(possible_keyword): raise notracebackException, "uniref2go_db_filename ERROR "+possible_keyword+' not found'
    return possible_keyword
  
  def blast_options_dict(self):
    """ Returns a dictionary with the blast options like option: value. """
    return self.options_dict(category='blast')
  
  def exonerate_options_dict(self):    return self.options_dict(category='exonerate')
  def genewise_options_dict(self):     return self.options_dict(category='genewise')
  def tag_blast_options_dict(self):    return self.options_dict(category='tag_blast')

  def options_dict(self, category):
    """ Utility to run blast_options_dict, exonerate_options_dict or any other with the same logic"""
    if keywords[category+'_options'].has_key(self.__dict__[category+'_options']):       options_string=keywords[category+'_options'][self.__dict__[category+'_options']]
    else:                                                                       options_string=self[category+'_options']
    out={}
    while options_string:
      try:
        ## fixing this code after years. Boy this code looks primitive. Putting a bad looking patch.
        s=options_string.find('-')+1
        assert s
        e=s+options_string[s:].find(' ')
        option_name    = options_string[s:e] #.split('-')[1].split()[0]
        options_string = options_string.split('-'+option_name)[1]
        if options_string.split('-')[0][-1]==' ':       value = del_white( options_string.split('-')[0]   ) 
        else:                                           value=options_string.split()[0]          
        options_string = join(options_string.split(value)[1:], value).strip()
        out[option_name]=value
      except: raise notracebackException, 'selenoprofiles ERROR parsing '+category+' options: '+str(options_string)
    return out

  def clustering_seqid_value(self):
    """ This takes the attribute self.clustering_seqid and translates it into something that can be used right away (it changes the potential label into a value)"""
    return self.options_non_str(category='clustering_seqid', istype=float)
  def max_column_gaps_for_blast_query_value(self):
    """ This takes the attribute self.max_column_gaps and translates it into something that can be used right away (it changes the potential label into a value)"""
    return self.options_non_str(category='max_column_gaps', istype=float)
  def max_blast_hits_number_value(self):
    """ This takes the attribute self.max_blast_hits and translates it into something that can be used right away (it changes the potential label into a value)"""
    return self.options_non_str(category='max_blast_hits', istype=int)

  def options_non_str(self, category, istype):
    """ generalized function to return a number"""
    if type(self.__dict__[category]) ==str:
      if keywords[category].has_key(self.__dict__[category]):        return keywords[category][self.__dict__[category]]
      else:  raise notracebackException, "selenoprofiles ERROR can't find keyword "+self.__dict__[category]+' for parameter: '+category
    elif type(self.__dict__[category]) ==istype:                        return self.__dict__[category]
    else: raise notracebackException, "selenoprofiles ERROR bad type for profile parameter "+category+': '+str(self.__dict__[category])

  def nseq(self):    return len ([t for t in self.titles() if t!='BLAST_QUERY_MASTER'])
  
  def summary(self, descriptor=''):
    """ This returns a summary of the attributes of the object. If descriptor is set as None, the first line of output (like >>> profile: profilename) is skipped """
    if descriptor!=  None and not descriptor: descriptor=self.name
    tab_width=120
    o='' #' .'+'.'*118+'.\n'
    if not descriptor == None: o+=(' .'+'-'*49+' Profile: '+descriptor+' ').ljust(tab_width-1, '-')+". \n"
    line=' | n_seqs: '+str(self.nseq())+' ; n_queries: '+str(len(self.queries_titles()))
    if self.sec_pos(): line+=' ; sec_position'+'s'*int(bool(len(self.sec_pos())>1))+': '+str(join([str(i) for i in self.sec_pos()], ', '))
    line+=' ; n_clusters: '+str(self.n_clusters())    
    o+=line.ljust(tab_width-1, ' ')+'|\n'    
    cp=self.conservation_properties(with_coverage=True); ncp=  self.conservation_properties(with_coverage=False)
    line=' | average AWSIc: '+str(round(cp[0], 3))+' +- '+str(round(cp[1], 3))+'  ;  average AWSIw: '+str(round(ncp[0], 3))+' +- '+str(round(ncp[1], 3))+' '
    o+=line.ljust(tab_width-1, ' ')+'|\n '
    options_lines='| '
    for a in profile_alignment.parameters:
      value=getattr(self, a); add=''
      if type(value)!=list:  values=[value]
      else: values= value
      for value in values:
        if value and value!='DEFAULT' and value!=[]: add=a+': '+str(value)+' | '
        if len(options_lines)+len(add)>tab_width-3: 
          o+=options_lines[:-3].ljust(tab_width-3, ' ')+' |\n '
          options_lines='| '+add
        else: options_lines+=add
    if options_lines!='| ':    o+=options_lines[:-3].ljust(tab_width-3, ' ')+' |\n'
    o+=' |'+(tab_width-3)*'-'+'|'
    return o

  __repr__=summary
  
  """
  ## debug
  def add(self, title, seq):
    if title.startswith('BLAST_'): 
      write('adding '+title+' '+seq, 1)
      #raise Exception, '...'
    alignment.add(self, title, seq)
  """    
#    """ Extend the normal add function, it can control if the sequence introduced is a query"""
#    alignment.add(self, title, seq)
#    self.queries.append(len(self.titles())-1) 
    
  def neutral_tags_list(self):
    """ Utility to return the list of neutral tags contained in self.neutral_tags but checking whether a keyword is provided: in this case, this is translated to the corresponding list defined in the config file """
    if type(self.neutral_tags)==list: return self.neutral_tags
    if not keywords['neutral_tags'].has_key(self.neutral_tags): raise notracebackException, "profile_alignment->neutral_tags_list ERROR can't find keyword  "+str(self.neutral_tags)+' for neutral_tags in the config file. A line like this must be present: neutral_tags.'+self.neutral_tags+" = ['tag1', 'tag2', 'tag3']"
    return keywords['neutral_tags'][self.neutral_tags]
  

##############################################################################################################################################
# DATABASE

from sqlite3 import dbapi2 as sqlite, OperationalError, Connection, Cursor
import time
class overriding_cursor(Cursor):
  """ just to replace the execute command with something lock-proof"""
  def execute(self, *args):
    """ ovverides the normal execute command allowing for waiting if the db is locked"""
    attempts=0;     success=False
    while not success and attempts <= max_attempts_database:    
      attempts+=1
      try:      
        out=Cursor.execute(self, *args)
        success=True
      except OperationalError:
        printerr('Cannot execute the command, the database is locked. I\'m waiting to see if it unlocks... max attempts remaining: '+str( max_attempts_database-attempts+1), 1)
        time.sleep(sleep_time)
    if not success:      raise notracebackException, "ERROR locked database"
    else:                return out

class selenoprofiles_db(sqlite.Connection):
  """ Class to manage the selenoprofiles database. """
  results_fields= [ 'id', 'profile', 'program', 'target_header', 'query_header', 'chromosome', 'alignment_query', 'alignment_target', 'label', 'state', 'awsi'] 

  def cursor(self):    return overriding_cursor(self)
  
  def original_cursor(self): return sqlite.Connection.cursor(self)

  def update_to_last_version(self):
    """if the database does not contain the newly introduced fields / tables, an ERROR is printed and the program exits"""
    #checking number of columns for results
    db_cursor=self.cursor()
    existing_fields= [i[1] for i in db_cursor.execute('pragma table_info(results);').fetchall()]
    if existing_fields != self.results_fields:       raise notracebackException, 'ERROR database was created with an older version of selenoprofiles. It is recommended to delete the results.sqlite file and rerun to populate the database with option -D'
    
  def has_table(self, name):
    """ Utility to check if the database as a table with such name"""
    try:
      db_cursor=self.original_cursor()
      db_cursor.execute('SELECT * FROM '+name)
      return True
    except sqlite.OperationalError:      return False

  def table_entry(self, table_name, **keyargs):
    """ Utility to return an entry from a table with certain characteristics, definable by key arguments.  """
    db_cursor=self.cursor()
    line='SELECT * FROM '+table_name+ ' WHERE '+join([k+'=="'+str(keyargs[k])+'"' for k in keyargs], " AND ")
    db_cursor.execute(line )
    all_results= db_cursor.fetchall()
    if len(all_results)>1: raise Exception, "selenoprofiles_database-> table_entry ERROR more than one entries found with conditions: "+join([k+'=="'+str(keyargs[k])+'"' for k in keyargs], " AND ")
    elif len(all_results)==0: return False
    return all_results[0]

  def initialise_db(self):
    """ Add necessary tables """
    db_cursor=self.cursor()
    db_cursor.execute('CREATE TABLE has_results_for_profile (id INTEGER PRIMARY KEY, profile VARCHAR(50), has_output VARCHAR(11))')  
    ## values for  has_output:  'NO' or entry not present -> doesn't have ouput  ;   'ONGOING'   -> computation is in progress right now  ;    'UNCHECKED' -> computation is over, but other profiles are running within the same instance of selenoprofile that completed thsi profile, before remove_redundancy method can be run;   'WAITING' -> actual computation is over, but a selenoprofiles instance for this profile is waiting another isntance before remove_redundancy can be run  ;          'YES': results are ready 
    db_cursor.execute('CREATE TABLE results (id INTEGER PRIMARY KEY, profile VARCHAR(50), program VARCHAR(12), target_header VARCHAR(1000), query_header VARCHAR(1000), chromosome VARCHAR(100), alignment_query VARCHAR(2000), alignment_target VARCHAR(2000), label VARCHAR(30) , state VARCHAR(30) , awsi FLOAT )')
    db_cursor.execute('CREATE TABLE features (id INTEGER PRIMARY KEY, resultid INTEGER,  type VARCHAR(15), text VARCHAR(2000) )')       
    
  def add_entry(self, table_name, entry_from_db, make_id_null=True):
    """ This add a result from another db in the table with the same name"""
    db_cursor=self.cursor()
    if make_id_null:  entry_to_put=entry_from_db[1:]
    else:             entry_to_put=entry_from_db
    parentesis_part= '('
    if make_id_null:       parentesis_part+='null, '
    for t in entry_to_put: parentesis_part+='?, '
    parentesis_part=parentesis_part[:-2]
    parentesis_part+=')'
    write('INSERT INTO '+table_name+' VALUES '+parentesis_part +' --- '+str(entry_to_put), 1)
    db_cursor.execute('INSERT INTO '+table_name+' VALUES '+parentesis_part , (entry_to_put))        

  def has_results_for_profile(self, profile):
    """ Return a value among NO YES ONGOING UNCHECKED WAITING UNCHECKED-2.    NO   is returned even if there's no entry in the db"""
    db_cursor=self.cursor()
    db_cursor.execute( 'SELECT * FROM has_results_for_profile WHERE profile=="'+profile+'"' ) 
    a= db_cursor.fetchone() 
    if not a: return 'NO'
    return a[-1]

  def set_has_results_for_profile(self, profile_name, value):
    """ Utility to set this property: has_results_for_profile. It finds automatically if there was an entry already and replace it.
    value should be a string, among NO YES ONGOING UNCHECKED WAITING UNCHECKED-2  """
    db_cursor=self.cursor()
    old_entry= self.table_entry('has_results_for_profile',profile=profile_name)  
    if old_entry:  #...if we already had an entry for this target and profile, we must replace it
      entry_id=old_entry[0]    
      db_cursor.execute('UPDATE has_results_for_profile SET  profile = ?, has_output=? WHERE id == '+str(entry_id), (profile_name, value))
    else:           db_cursor.execute('INSERT INTO has_results_for_profile VALUES ( null , ?, ?)', (profile_name, value))

  def set_stopped_profiles_as_ok(self):
    """ This is used in combination with option -stop: it will set the values of the profiles that were with option -stop as checked"""
    db_cursor=self.cursor()
    db_cursor.execute('UPDATE has_results_for_profile SET  has_output="YES" WHERE has_output == "UNCHECKED-2"')          

  def computation_in_progress(self, apart_from=[]):
    """ Returns True if the database reports of any   has_results_for_profile   property which is either UNCHECKED or ONGOING, which means that computation is in progress in this target."""
    db_cursor=self.cursor()
    line='SELECT * FROM has_results_for_profile WHERE has_output == "UNCHECKED" OR has_output == "ONGOING"'
    db_cursor.execute(line);     all_entries= db_cursor.fetchall()
    for i in all_entries:
      if not i[-2] in apart_from:  return True
    return False       

  def save(self):    
    """ saves the changes (commit). it tries more than once if it finds locked  """
    attempts=0;   success=False
    while not success and attempts <= max_attempts_database:    
      attempts+=1
      try:      
        out=self.commit()
        success=True
      except OperationalError:
        printerr('Cannot commit changes, the database is locked. I\'m waiting to see if it unlocks... max attempts remaining: '+str( max_attempts_database-attempts), 1)
        time.sleep(sleep_time)
    if not success:      raise notracebackException, "ERROR the database was found locked in all "+str(max_attempts_database)+ " attempts!"
    else:                return out

#  def unlock(self):
#    """If you have errors like: database table is locked but you know it isn't, this should fix it """  ##it doens't!
#    self.commit()
#    self.close()

  def add_result(self, p2g_hit):
    """ Add a result to the database, annotating its data (prediction program etc), and also  the chromosome is in. The majority of information is in the header"""
    state=p2g_hit.filtered
    if state=='redundant': state+='_'+str(p2g_hit.overlapping.id)
    db_cursor=self.cursor()
    target_header= p2g_hit.header(no_species=True, no_target=True, no_chromosome=True)
    query_header=  p2g_hit.query.header(no_id=True, no_species=True, no_target=True, no_strand=True)
    db_cursor.execute('INSERT INTO results VALUES (null  , ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', (p2g_hit.profile.name, p2g_hit.prediction_program(),  target_header, query_header,  p2g_hit.chromosome,  p2g_hit.alignment.seq_of('q'), p2g_hit.alignment.seq_of('t'), p2g_hit.label,  state, p2g_hit.weighted_seq_identity_with_profile() ))    
    for obj in p2g_hit.features:
      self.add_feature_to_result(obj, p2g_hit)
      
  def add_feature_to_result(self, obj, p2g_hit):
    """ Add a feature (obj) to this result in the database """
    db_cursor=self.cursor()
    result_id=self.result_id(  profile_name=p2g_hit.profile.name, output_index=p2g_hit.id)    
    feature_name= obj.__class__.__name__
    text=obj.dump_text()
    db_cursor.execute('INSERT INTO features VALUES (null, ?,  ?, ?)', (result_id, feature_name, text))

  def result_id(self, profile_name, output_index):
    """ Utility to parse the db and return the database id of a certain prediction. the argument output_index is:  5 in the example prediction sps.5.selenocysteine """
    db_cursor=self.cursor()
    db_cursor.execute( 'SELECT * FROM results WHERE profile=="'+profile_name+'" AND target_header LIKE "'+str(output_index)+' %"' ) 
    all_results= db_cursor.fetchall()
    if len(all_results)>1: 
      raise Exception, "selenoprofiles_database-> result_id ERROR redundant entries in results table; check:  "+profile_name+'.'+str(output_index)
    elif len(all_results)==0: return None
    return all_results[0][0]
    
  def clear_results(self, profile_name): 
    """ Clean the table from the entries relative to a certain profile. """
    db_cursor=self.cursor()
    db_cursor.execute('SELECT id FROM results WHERE profile == "'+profile_name+'" ')
    ids_removed=db_cursor.fetchall()
    db_cursor.execute('DELETE FROM results WHERE profile == "'+profile_name+'" ')
    for object_id in ids_removed:       db_cursor.execute('DELETE FROM features WHERE resultid == "'+str(object_id[0])+'" ')      
    
  def get_results(self, profile_name, states=[]):
    """ Utility to get all results in the database relative to this profile """
    db_cursor=self.cursor()
    sql_cmnd='SELECT * FROM results WHERE profile == "'+profile_name+'" '
    if states: 
      sql_cmnd+=' AND ( '
      for state in states:
        if   state =='redundant':            sql_cmnd+=' state LIKE "redundant_%" OR'
        elif state =='overlapping':          sql_cmnd+=' state LIKE "%_overlapping_%" OR'
        else:                                sql_cmnd+=' state == "'+state+'" OR'
      sql_cmnd = sql_cmnd[:-2]+')' #removing last OR
    db_cursor.execute(sql_cmnd)
    return db_cursor.fetchall()

  def get_result(self, profile_name, index_id):
    """ Utility to get a single result from the database. The index_id must be specified """
    db_cursor=self.cursor()
    db_cursor.execute('SELECT * FROM results WHERE profile == "'+profile_name+'" AND  target_header LIKE "'+str(index_id)+' %"')
    all_results= db_cursor.fetchall()
    if len(all_results)>1:    raise Exception, "selenoprofiles_database-> result_id ERROR redundant entries in results table; check:  "+profile_name+'.'+str(index_id)
    elif len(all_results)==0: return None
    return  all_results[0]

  def get_results_on_chromosome(self, chromosome, strand=None):
    """ Utility to get all results in the database on a certain chromosome (or scaffold). If strand is defined to "+" or "-", only those in the desired strand are reported """
    db_cursor=self.cursor()
    strand_add=''
    if strand in ['+', '-']: strand_add=' AND target_header LIKE "% strand:'+strand+' %"'
    db_cursor.execute('SELECT * FROM results WHERE chromosome == "'+chromosome+'" '+strand_add)
    return db_cursor.fetchall()

  def all_chromosomes_with_predictions(self):
    """ Returns a list of the chromosomes with something predicted"""
    db_cursor=self.cursor()
    db_cursor.execute('SELECT DISTINCT chromosome FROM results ')
    db_results=db_cursor.fetchall()
    if not db_results: return []
    return [   r[0] for r in db_results]

  def features_of(self, result_db_id, add_parent=False):
    """ Returns the list of all features as python objects associated in the db with the desired result object, given its db id. 
    With add_parent,  you can provide a p2ghit object that will be linked in the .parent attribute of feature object. """
    list_out=[]
    db_cursor=self.cursor()
    db_cursor.execute('SELECT * FROM features WHERE resultid == "'+str(result_db_id)+'"')  
    for db_feat in db_cursor.fetchall():
      feat_id, result_id, feat_type, text = db_feat
      try:        feat_class= eval(feat_type)
      except NameError:      raise Exception, "ERROR feature '"+feat_type+"' found in database is not recognized!"
      f= feat_class();   
      if add_parent:        f.parent=add_parent
      try:        f.load_dumped_text(text)
      except:     raise Exception, "ERROR while loading dumped text in database for feature of result with db_id:  "+str(result_db_id)+' ; text failed to load:\n'+str(text)
      list_out.append(f)
    return list_out    

  def remove_overlapping_states(self, silent=False):
    """ Clean the database from all interfamily overlap calls """
    if not silent: 
      write('DATABASE', how=terminal_colors['database']); write(' removing overlap calls. ');    write('fetching overlaps... ')
    db_cursor=self.cursor()
    db_cursor.execute('SELECT id, state FROM results WHERE state LIKE "%_overlapping_%"') 
    if not silent:     write('now removing overlaps... ', 1)
    n_removed=0
    for obj in db_cursor.fetchall():
      db_entry_id = obj[0]  
      state =       obj[1]  
      state_no_overlap=  state.split('_overlapping_')[0]
      db_cursor.execute('UPDATE results SET state="'+state_no_overlap+'" WHERE id=="'+str(db_entry_id)+'"')
      n_removed+=1  
      if not silent: service(str(n_removed)+' overlaps removed')
    self.save()        
    if not silent:  write('DATABASE', how=terminal_colors['database']); write(' removed '+str(n_removed)+' overlap calls.', 1)  

  def remove_redundancy(self, chromosome_list=None, silent=False):
    """ Parse all results in the database for this target and removes the overlapping ones (updating their state attribute to "redundant", deciding which to keep with the usual choose_prediction function"""
    db_cursor=self.cursor()
    if chromosome_list is None: chromosome_list=self.all_chromosomes_with_predictions()
    for chromosome in chromosome_list:
      for strand in ['+', '-']:
        results_on_this_chromosome=[] #and strand, actually
        for db_result in self.get_results_on_chromosome(chromosome, strand=strand):
          p2g_result=load_p2g_from_db_result(db_result) 
          p2g_result.db_id=db_result[0]
          if not p2g_result.filtered in  ['redundant', 'overlapping']:    results_on_this_chromosome.append( p2g_result ) #loading all results on this chromosome
        removed_overlapping_hits=[]
        non_redundant_genes= remove_overlapping_genes(results_on_this_chromosome,  cmp_fn=choose_among_overlapping_p2gs_interfamily,  phase=False, out_removed_genes=removed_overlapping_hits, remember_overlaps=True) ## phase==False for gene with frameshifts
        for p2g_result in removed_overlapping_hits:
          if not silent:        write('DATABASE', how=terminal_colors['database']); write(' remove redundancy: '+p2g_result.output_id()+' removed since found overlapping with '+p2g_result.overlapping.output_id(), 1)
          db_cursor.execute('UPDATE results SET state=? WHERE id == ? ', (p2g_result.filtered+'_overlapping_'+str(p2g_result.overlapping.profile_name)+'.'+str(p2g_result.overlapping.id), str(p2g_result.db_id)) )

##############################################################################################################################################
####GENE SUBCLASSES
##########
#P2GHIT

class p2ghit(gene):
  """ This class includes blasthit, exoneratehit and genewisehit class. It contains the code for all operations that we want to perform on any selenoprofiles result, so it is the mother class for output. It is a gene class which contains a .query attribute which is a second gene instance. """
  def __init__(self):
    gene.__init__(self)
    self.query=gene()
    self.alignment=alignment()
    self.cds_sequence=''
    self.profile=''
    self.label=''
    self.tag_score_data=None
    self.go_score_data =None
    self.features=[]
    self.sequence_identity_with_profile_data={}
    self.awsi_data={}
  
  def reset_derived_data(self):
    """ To be called everytime a prediction set of coordinates is modified"""
    self.cds_sequence=''
    self.awsi_data={};     self.sequence_identity_with_profile_data={}
    for f in self.features:
      if hasattr(f, 'reset'):       f.reset()    
  
  def summary(self, description='', **keyargs):
    if not description:   description=self.__class__.__name__.upper()
    return gene.summary(self, description, **keyargs)+'\n'+gene.summary(self.query, description+'_query', **keyargs)
  
  def cds(self):
    """ Returns the coding sequence for the prediction. If frameshifts are predicted, these inserted nucleotides are not returned: the sequence returned always translates to the predicted protein sequence (considering SecTGAs) """
    ##### Although the coding sequence is completely obtainable in lazy computing from the misc_sequence_data, the codign sequence is available even before populating the misc_sequence_data, since it is parsed from exonerate and genewise output. So I'm keeping the .cds_sequence attribute just to take advantage of this, and not perform any subsequence fetch if user wants only the coding sequence
    if not self.__dict__.has_key('cds_sequence') or not self.cds_sequence: 
      self.cds_sequence= self.subsequence(1, self.length() )
    return self.cds_sequence

  def dna(self):
    """ Returns the full DNA sequence for the prediction, including introns and eventual frameshifts."""
    return self.subsequence(1, self.span(), include_introns=True)

  def aligned_cds(self):
    """ Returns the cds like the previous function, but with gaps in positions as the .alignment attribute. This is useful if you want to index the cds with indexes coming from the aligned sequences."""
    seq=self.cds()
    for pos in range(len(self.alignment.seq_of('t'))):
      if self.alignment.seq_of('t')[pos]=='-':        seq=seq[:pos*3]+'---'+seq[pos*3:]
    return seq
  
  def frameshifts(self):
    """ Return the number of frameshifts found in this prediction """
    introns_g= self.introns(minimal=True);     count=0
    for st, end in introns_g.exons: 
      if end-st+1 <=5: count+=1
    return count

  def protein(self):
    """ Returns the protein sequence predicted in the target"""
    return nogap( self.alignment.seq_of('t') )

  def is_complete_at_three_prime(self):
    """ returns True if the next codon is a stop codon """
    seq_until_stop, translation, length_of_three_prime  = self.three_prime_coding_tail()    
    if translation and translation[0]=='*': return True
    return False
  
  def output_header(self):
    """ used in all fasta files """
    add=''
    if opt['fasta_add']:
      x = self
      try: add= ' '+str(eval(opt['fasta_add']))
      except: 
        printerr("ERROR with fasta_add option: can't evaluate: "+opt['fasta_add'], 1)
        raise
    return self.header(no_id=True, function_prediction_program=1)+add

  def output_fasta(self):
    """ Returns the text for its fasta file"""
    return ">"+self.output_id()+' '+self.output_header()+'\n'+fasta(self.protein())
  fasta=output_fasta
  
  def output_cds(self):
    """ Returns the text for a cds fasta file"""
    seq= self.cds()
    seq_until_stop, translation, length_of_three_prime  = self.three_prime_coding_tail()      #    if self.is_complete_at_three_prime():
    if translation and translation[0]=='*':     seq+=seq_until_stop
    return ">"+self.output_id()+' '+self.output_header()+'\n'+seq

  def output_dna(self):
    """ Returns the text for a dna fasta file"""
    seq=self.dna()
    seq_until_stop, translation, length_of_three_prime  = self.three_prime_coding_tail()      #    if self.is_complete_at_three_prime():
    if translation and translation[0]=='*':     seq+=seq_until_stop
    return ">"+self.output_id()+' '+self.output_header()+'\n'+seq

  def output_aligned_cds(self):
    """ Returns a pairwise alignment in fasta between the coding sequence and a fake back translation of the blast_query_master where gaps in the original alignment have been converted to X, while new gaps are as -.  """
    t='UnMaTcHaBlE'
    prot_profile_ali=self.alignment_with_profile(dont_shrink=True, title=t)
    c_pos_target=-1 ; seq_query='' ; seq_target='' 
    for pos in range(  prot_profile_ali.length()  ):
      if prot_profile_ali.seq_of(t)[pos]!='-': 
        c_pos_target+=1
        seq_target+= self.cds()[c_pos_target*3:c_pos_target*3+3]      
      else:
        seq_target+='---'
      all_gaps_apart_in_target = prot_profile_ali.seq_of(t)[pos]!='-' and all (  [ prot_profile_ali.seq_of(title)[pos]=='-'  for title in prot_profile_ali.titles() if title !=t ]   )
      if all_gaps_apart_in_target: seq_query+='---'
      else:
        seq_query+= retrotransl( prot_profile_ali.seq_of('BLAST_QUERY_MASTER')[pos], gaps_to='XXX')
    return ">"+self.output_id()+' '+self.output_header()+'\n'+seq_target+'\n>BLAST_QUERY_MASTER backtranslated with Xs instead of original gaps\n'+seq_query

  def output_introns(self):
    """Returns the text for the introns output file """
    if len(self.exons)<2: return None
    out=''; fs_found=0
    introns=self.introns()
    for intron_index in range(len(self.exons)-1):
      intron=introns.get_exon(intron_index)
      if intron.length()>5:        
        out+=">intron_"+str(intron_index+1-fs_found)+'_of:'+self.output_id()+' '+self.output_header()+'\n'+intron.fasta_sequence()[1]+'\n'
      else: fs_found+=1
    return out

  def output_id(self):
    """ Returns the base for output: e.g. sps.1.selenocysteine"""
    try:      return self.profile.name+'.'+str(self.id)+'.'+self.label
    except:   return self.profile_name+'.'+str(self.id)+'.'+self.label      
    
  def output_feature(self, feature_name):
    """ Generic function for outputing a feature, like secises """    
    out=''
    for obj in self.features:
      if obj.__class__.__name__ == feature_name:         out+= obj.output()+'\n'
    return out
    
  def output_secis(self):    return self.output_feature('secis')
    
  def stop_codons(self):
    """ Returns a list of elements like[position_in_prot_seq_1_based, codon_in_upper_RNA_letters] """
    out=[]
    for pos in range(len( self.protein() )):
      if self.protein()[pos] in '*':
        codon=replace_chars(upper(self.cds()[pos*3:pos*3+3]), 'T', 'U')
        if codon in STOP_CODONS:   out.append([pos+1, codon])
        else: 
          #print self
          raise Exception, "ERROR the codon for \"*\" is not a stop codon: "+codon+' ; result_id:'+str(self.id)+' ; prediction_program: '+self.prediction_program()+' GFF: \n'+self.gff()
    return out

  def n_stop_codons(self):
    """ Returns the number of stop codons. This was implemented to speed up procedures in which you just care about how many stop codons you have, and not how many do you have."""
    return self.protein().count('*')
  
  def aligned_secTGAs(self):
    """ Returns the number of aligned secTGAs in the query """
    secTGAs=0
    for pos in  self.alignment.all_positions_of('U'):
      if not self.alignment.seq_of('t')[pos] in '-x': secTGAs+=1
    return secTGAs
  
  def prediction_program(self):
    """ Recognize the type of p2ghit and returns the program used for the prediction """
    if not self: return 'None'
    if issubclass(self.__class__, blasthit):        return 'blast'
    elif issubclass(self.__class__, exoneratehit):  return 'exonerate'
    elif issubclass(self.__class__, genewisehit):   return 'genewise'
    else: raise Exception, "ERROR can't determine prediction program for p2ghit object: "+str(self)

  def query_full_name(self):
    """ Return the full name of the query, as it appears in the profile. This differs from self.query.chromosome by the fact that this is sometimes only the first word of the full title """
    return self.profile.fill_title(self.query.chromosome)
  
  def coverage(self):
    """ Returns the coverage respect to the profile, i.e. how much profile is spanned"""
    if not self.__dict__.has_key('profile') or (not self.profile):  raise Exception, "coverage ERROR no profile alignment is defined for this p2g hit. Please add this as its property "
    prediction_start_in_ali, prediction_end_in_ali =       self.boundaries_in_profile()
    coverage_value= float(prediction_end_in_ali-prediction_start_in_ali+1)/self.profile.length()
    return coverage_value    

  def ali_fasta_title(self):
    """Returns the fasta title normally displayed for the prediction in the .ali file """
    return self.output_id()+' '+self.output_header()   
  
  def filled_alignment(self, target_title=''):
    """ Returns a pairwise alignment like the one in .alignment, but filled with the full query sequence (self.alignment contains only the query sequence which is found in the prediction file. This method returns an alignment which contains it all) 
    The titles are: 1.the full title for the query , 2.the output title for the target (=   self.output_id()+' '+self.header(no_id=True, function_prediction_program=1)   ), unless a different target title is specified as argument
    """
    ali=alignment()
    query_full_title=self.profile.fill_title(self.query.chromosome)
    query_seq=self.alignment.seq_of('q') 
    if not target_title: target_title=self.ali_fasta_title()
    ali.add(query_full_title, query_seq)
    ali.add(target_title, self.alignment.seq_of('t'))
    ali=ali.fill_sides( self.profile, wild_chars='UX*')
    return ali    

  def alignment_with_profile(self, profile_ali='', dont_shrink=False, title=''):
    """ This returns the alignment for the p2g hit along with the profile alignment, produced by the transfer alignment method of alignment class. A profile_ali must be defined either as argument or be present as .profile attribute of the p2g_hit. The alignment returned contains the same titles as the profile_ali, plus the standard fasta title containing the protein sequence of the target.
    """
    if self.__dict__.has_key('profile') and not profile_ali:        profile_ali=self.profile
    elif not self.__dict__.has_key('profile') and not profile_ali:  raise Exception, "p2ghit->alignment_with_profile ERROR no profile alignment is defined for this p2g hit. Please add this as its property "
    ali=self.filled_alignment(target_title=title)
    return ali.transfer_alignment( profile_ali, dont_shrink=dont_shrink )

  def sequence_identity_with_profile(self, dont_count_gaps=0):
    """This maps the prediction back to the profile and returns the average sequence identity with the sequences of the profile. dont_count_gaps is flag for the following behavior:
    0 -> terminal gaps are not counted
    1 -> gaps are not counted
    2 -> everything is counted (resulting score will be much lower than the other methods for uncomplete hits)
  """
    if not self.sequence_identity_with_profile_data.has_key(dont_count_gaps):
      dont_count_gaps=int(dont_count_gaps)
      dcg= bool(dont_count_gaps==1)
      dctg=bool(dont_count_gaps==0)
      self.sequence_identity_with_profile_data[dont_count_gaps]=self.alignment_with_profile(dont_shrink=True, title='UnMaTcHaBlE').average_sequence_identity_of(title='UnMaTcHaBlE', dont_count_gaps=dcg, dont_count_terminal_gaps=dctg)
    return self.sequence_identity_with_profile_data[dont_count_gaps]

  seq_id=sequence_identity_with_profile

  def weighted_seq_identity_with_profile(self, with_coverage=True):
    """ Scores this prediction against the rest of the alignment, returning a float from 0.0 to 1.0; also called AWSI - awsi. See function with same name in MMlib.
    The score is based on pairwise comparisons, in which a WSI score is computed as a sequence identity with weights on the different columns of the profiles. During a pairwise comparison between a profile sequence and a candidate sequence, the weight is given by the representation of the aminoacid in this profile sequence across all the profile. More conserved positions are given more weight. Unless with_coverage=False, this weight is also multiplied by the column coverage, that is to say, the number of total characters which are not gaps divided by the total number of profile sequences.
    """
    if not self.awsi_data.has_key(with_coverage):
      try:
        if self.profile.nseq()==1:       
          a=self.alignment_with_profile(dont_shrink=True, title='UnMaTcHaBlE')          
          out=a.sequence_identity_of(0, 1)  
        else:      
          cons_map_profile=self.profile.conservation_map()
          coverage_per_position=[]  #meaning: the proportion of seqs with something not a gap / tot n seqs 
          for element in cons_map_profile:
            coverage=  sum([  element[k]   for k in element if k !='-'  ])
            coverage_per_position.append(coverage)
          a=self.alignment_with_profile(dont_shrink=True, title='UnMaTcHaBlE')  # getting alignment of this prediction with the profile
          columns_all_gaps={}
          for pos in range(a.length()):
            is_gap=True
            for title in a.titles():
              if is_gap and title!='UnMaTcHaBlE':   ## is_gap is to make the loops as efficient as possible
                is_gap= a.seq_of(title)[pos]=='-' 
            if is_gap and a.seq_of('UnMaTcHaBlE')[pos]!='-':         columns_all_gaps[pos]=1  # 0 based !   self check to avoid problems but if alignemnt doens't contain unaligned columns it's not a problem
          scores=[];     
          for t in self.profile.titles():
            if t=='BLAST_QUERY_MASTER': continue
            total_weight=0; score_this_seq=0
            pos_in_ali=-1    ### position in original profile alignment (0-based)
            for pos in range(a.length()):
              if not columns_all_gaps.has_key(pos):          pos_in_ali+=1
              else:                                          continue
              #now pos_in_ali reflects position we're in, 1 based, in the original alignment with no seq added.        
              candidate_seq_here=a.seq_of('UnMaTcHaBlE')[pos]
              profile_seq_here=  a.seq_of(t)[pos]
              if  profile_seq_here =='-': continue 
              if not with_coverage and candidate_seq_here=='-': continue        
              #print self.output_id(), pos, pos_in_ali, candidate_seq_here, profile_seq_here
              weight=       cons_map_profile [pos_in_ali][profile_seq_here] 
              if with_coverage: weight*= coverage_per_position[pos_in_ali]
              is_identity=  candidate_seq_here == profile_seq_here
              score_this_seq+=  int( is_identity  )*weight
              total_weight+=weight
            if total_weight:   scores.append(       score_this_seq/float(total_weight)     )

          if not len(scores): raise Exception, "p2ghit-> weighted_seq_identity_with_profile ERROR the sequence "+self.output_id()+' has not any aligned position with any of the other sequences!'
          out=sum(scores)/float(len(scores))
        self.awsi_data[with_coverage]= out
      except:        raise   #just for debugging, add msg here if you want
    return     self.awsi_data[with_coverage]
  awsi=weighted_seq_identity_with_profile

  def awsi_z_score(self, with_coverage=True):
    """ min: -999  max 999 """
    awsi =   self.weighted_seq_identity_with_profile( with_coverage )
    average, std_dev=   self.profile.conservation_properties( with_coverage )
    if std_dev==0.0: 
      if awsi< average:    z_scorez=-999.0
      elif awsi> average:  z_scorez=999.0
      else:                z_scorez=0.0
    else:
      z_scorez= (awsi-average) / std_dev
      if z_scorez<-999:  z_scorez=-999.0
      elif z_scorez>999: z_scorez=999.0
    return z_scorez
  z_score= awsi_z_score
  
  def awsi_filter(self, awsi=0.9, z_score=-3,  few_sequences_awsi=0.3, with_coverage=True):
    """ default awsi-based filter. Predictions with values of awsi or z_score greater than arguments pass the filter. by default AWSIc is computed. with_coverage=False activates AWSIw (see manual).  If the profile contains only one or two sequences, the filter checks that AWSI is greater than argument few_sequences_awsi """
    if self.profile.nseq()<3: return self.weighted_seq_identity_with_profile( with_coverage )  >= few_sequences_awsi
    if self.weighted_seq_identity_with_profile( with_coverage ) >= awsi: return True
    return self.awsi_z_score(with_coverage=with_coverage)>z_score

  def most_similar_sequence(self, dont_count_gaps=0):
    """ Returns the title of the most similar sequence in the profile, computed with dont_count_gaps as in sequence_identity_with_profile""" 
    dont_count_gaps=int(dont_count_gaps)
    dcg= bool(dont_count_gaps==1)
    dctg=bool(dont_count_gaps==0)    
    ali=self.alignment_with_profile(dont_shrink=True, title='UnMaTcHaBlE')
    max_title=''; max_seqid=0.0
    for title in ali.titles():
      if title!='UnMaTcHaBlE':
        seqid= ali.sequence_identity_of('UnMaTcHaBlE', title,  dont_count_gaps=dcg, dont_count_terminal_gaps=dctg)
        if seqid>max_seqid: 
          max_title=title; max_seqid=seqid
    if max_title: return max_title    

  def sequence_identity_in_range(self, start, end, dont_count_gaps=0, min_aligned=1):
    """Returns the average sequence identity computed only along the profile range start,end both 1-based and included.
 dont_count_gaps is flag for the following behavior:
    0 -> terminal gaps are not counted
    1 -> gaps are not counted
    2 -> everything is counted (resulting score will be much lower than the other methods for uncomplete hits)
   min_aligned is the minimum number fo aligned residues to return something != 0.0 ; this is a trick to avoid returning a high seq id when just a single or few aminoacids are aligned.                                                                                                        """
    dont_count_gaps=int(dont_count_gaps)
    dcg= bool(dont_count_gaps==1)
    dctg=bool(dont_count_gaps==0)
    a=self.alignment_with_profile(dont_shrink=True, title='UnMaTcHaBlE')
    start_cut=None; end_cut=None # 0 based
    pos_in_ali=0
    for pos in range(a.length()):
      all_gaps_in_this_column = all( [a.seq_of(title)[pos]=='-'   for title in a.titles() if title!='UnMaTcHaBlE'])
      if not all_gaps_in_this_column: pos_in_ali+=1
      #now pos_in_ali reflects position we're in, 1 based, in the original alignment with no seq added.
      if start == pos_in_ali and start_cut is None:    start_cut=pos
      if end   == pos_in_ali : 
        end_cut= pos
        break
    if start_cut is None or  end_cut is None: raise Exception, "sequence_identity_in_range ERROR start / end called are not valid: "+str(start)+','+str(end)+' ; alignment is long: '+str(self.profile.length())
    cut= a.columns(start_cut, end_cut-start_cut+1)
    if len(nogap(cut.seq_of('UnMaTcHaBlE')))>=min_aligned: return cut.average_sequence_identity_of(title='UnMaTcHaBlE', dont_count_gaps=dcg, dont_count_terminal_gaps=dctg)
    else: return 0.0

  def conservation_score(self):
    """ Returns a score computed using the conservation_score method of the alignment class in MMlib: """+alignment.conservation_score.__doc__
    return self.alignment_with_profile(title='UnMaTcHaBlE').conservation_score(title='UnMaTcHaBlE')
  
  def gff(self, is_gtf=False):
    """ Returns a gff for this prediction. it includes the positions of selenocysteine in the line of the exon where it belongs"""
    g=self.copy()
    g.id=self.output_id()
    g.program='selenoprofiles_'+self.prediction_program()
    comment=''
    sec_pos_features=[]
    u_pos_in_ali_list=self.alignment.all_positions_of('U')
    for sec_index, u_pos_in_ali in  enumerate(u_pos_in_ali_list):
      homologue_residue_pos= self.alignment.position_in_seq( 't', u_pos_in_ali+1  )  #0-based to 1-based
      if homologue_residue_pos:  #exluding super-rare case in which sec is in the ali but as a Nterm tail
        residue_type=      self.protein()[homologue_residue_pos-1] 
        first_pos_codon=    self.subseq(   (homologue_residue_pos-1)*3 + 1 , 1 , minimal=True).exons[0][0]  #first pos of codon for hom residue
        sec_pos_features.append(  [first_pos_codon, "sec_position."+str(sec_index+1)+":"+str(first_pos_codon)+'-'+residue_type  ]   )
      
    if self.is_complete_at_three_prime():
      if   self.strand=='+':  g.extend(right=3, inplace=True)
      elif self.strand=='-':  g.extend(left=3,  inplace=True)
        
    out= gene.gff(g, tag='cds', is_gtf=is_gtf, comment=comment, position_features=sec_pos_features)
    ### adding features to gff
    for obj in self.features:
      if obj.included_in_gff:
        if 'program' in dir(obj): program= obj.program
        else: program='generic_program'
        if 'gff_add' in dir(obj) and obj.gff_add: 
          comment=''
          for k in obj.gff_add:  ## each k element is an attribute name, or a method
            if k.startswith('M:'): #it's a method
              k=k[2:]
              value= eval("obj."+k+'()')
            else: 
              value= eval("obj."+k)
            comment+=k+':'+str(value)+' '
          comment=comment[:-1]       
        out+='\n'+obj.gff( program=program, tag = obj.__class__.__name__,   is_gtf=is_gtf, comment=comment )
    return out
  output_gff=gff
  
  def gtf(self):
    """ Returns a gtf3 text for this prediction"""
    return self.gff(is_gtf=True)
  output_gtf=gtf
  
  def seq_in_profile_pos(self, pos):
    """ Return the aa in the target aligned to the position pos of the profile alignment, 1 based!"""
    if not self.__dict__.has_key('profile') or (not self.profile):  raise Exception, "profile_alignment->seq_in_profile_pos ERROR no profile alignment is defined for this p2g hit. Please add this as its property "
    query_full_name  = self.query_full_name()
    pos_in_query_seq = self.profile.position_in_seq(query_full_name , pos)   # 1 based
    pos_in_self_ali  = self.filled_alignment().position_in_ali(query_full_name, pos_in_query_seq) # 1 based 
    return self.filled_alignment().seq_of(1)[pos_in_self_ali-1]
  
  def sec_is_aligned(self, to_chars='UC', query_aa='U'):
    """ Returns True if at least one U in the query (or another aa provided with query_aa) is aligned to either a U or C (or any other list of characters provided as to_chars argument). If a * is detected in the target in a candidate position, runs internally place_selenocysteine to be able to evaluate."""
    for pos in range(self.alignment.length()):
      if self.alignment.seq_of('q')[pos] in query_aa+"*":
        if 'U' in to_chars and  self.alignment.seq_of('t')[pos] == '*':            self.place_selenocysteine()
        if self.alignment.seq_of('q')[pos] in query_aa  and  self.alignment.seq_of('t')[pos] in to_chars: return True
    return False


  def introns_summary(self):
    """ Returns a string with a summary of the positions of the introns, if any are present. The positions are cds-based, 1-based and denote the last position of each exon. example:   '123,456' for a 3 exons gene, '' for a single exon gene """
    out=''
    if len(self.exons)>1:
      length_cds_so_far=0
      for start, end in self.exons[:-1]:
        length_cds_so_far+= ( end - start + 1 )
        out+=str(length_cds_so_far)+','
      out=out[:-1]
    return out

  def sec_positions_list(self):
    """ This function returns what was the sec_pos attribute, that is to say, a list of positions of all the Sec aminoacids relative to the coordinates along the target chromosome. The position is precisely the one of the first position of the UGA codon.  The positions are guessed from the positions of U letter in the alignment inside the p2ghit object      """
    sec_pos=[];     aa_seq=nogap(self.alignment.seq_of('t'))
    pos_u_in_aa_seq=aa_seq.find('U')
    while pos_u_in_aa_seq!=-1:
      g_subseq=self.subseq(1+ pos_u_in_aa_seq*3, 3 )
      if self.strand=='+':   sec_pos.append(   min(g_subseq.boundaries()) )
      elif self.strand=='-':    sec_pos.append( max(g_subseq.boundaries()) )
      pos_u_in_aa_seq=aa_seq.find('U',    pos_u_in_aa_seq+1)
    return sec_pos

  def tag_score(self, silent=False, max_n_titles=0, neutral=False, verbose=False):
    """ This function computes the tag score specific to this prediction, which is computed by running blastp with the predicted protein sequence against a reference database, typically uniref, and parsing the results using profile-defined tags.
    Blast is run only if necessary (output file not present) using function blast_against_tag_db (see below).
    For each blast hit, the full title of the subject and its evalue are considered.
    Initially the title is run a set of neutral tags. If any of these are matching, this blast hit is not considered for the final score. Then, the set of positive tags (profile defined tags) are tried to match. If any of these match, a positive score is assigned to this blast hit. If not, a negative score is assigned. The absolute value of the score attributed is the negative logarithm of the evalue.
    The final score is the sum of all scores attributed to all blast hits.    """    
    if self.tag_score_data is None:
      if self.profile.tags:
        positive= self.profile.tags
        if not neutral: neutral=  self.profile.neutral_tags_list()
        tag_blast_outfile= self.blast_against_tag_db(silent=silent)
        names_and_evalues= parse_blast_subject_and_evalues(tag_blast_outfile)
        if max_n_titles==0:      max_n_titles=len(names_and_evalues)
        score=0
        #print positive
        #for p in positive: print "***", p, "***"
        for title, evalue in names_and_evalues[:max_n_titles]:
          score_of_this_tag=shortcut_log(evalue)          
          if   match_any_word(title, positive):      pass
          elif match_any_word(title, neutral):       score_of_this_tag=0
          else:                                      score_of_this_tag=-score_of_this_tag
          score+=score_of_this_tag          
          if verbose: write(title+' --->'+str(score_of_this_tag), 1)          
        self.tag_score_data=score
    return self.tag_score_data
  
  def go_score(self, silent=False, max_n_titles=0, verbose=False, terms=[]):
    """ This function computes the GO (gene ontology) score of this prediction, which si computed running blast with the predicted protein sequence against a reference database, typically uniref, and parsing the results using a set of user-defined GO terms associated to this family.
    Blast is run only if necessary (output file not present) using function blast_against_tag_db (see below).
    For each blast hit, the full title of the subject and its evalue are considered; if the protein has no functional GO term, it is treated as neutral, i.e., it is assigned a score of 0.
    If it has some, then the profile-GO are matched with those of this protein (and with the parents of those annotated for this protein). If any is found, the title is assigned a positive a score, while if it hasn't any, a negative score is assigned.
    The absolute value of the score attributed is the negative logarithm of the evalue.
    The final score is the sum of all scores attributed to all blast hits.    """    
    if self.go_score_data is None:
      if not terms: terms= self.profile.go_terms
      if terms:       
        if not is_file(self.profile.uniref2go_db_filename()):  raise notracebackException, 'selenoprofiles-> go_score ERROR uniref2go database file not found: '+self.profile.uniref2go_db_filename()+'  This file is necessary to call the method go_score. See installation script.'
        tag_blast_outfile= self.blast_against_tag_db(silent=silent)
        names_and_evalues= parse_blast_subject_and_evalues(tag_blast_outfile)
        if not names_and_evalues:   total_score=0
        else:
          if max_n_titles==0:      max_n_titles=len(names_and_evalues)
          #writing ids of the blast titles found in nr in a temporary file
          id_list=[]      
          for name, evalue in names_and_evalues[:max_n_titles]:        id_list.append(  name.split()[0]  )
          id_list_file=temp_folder+'id_list'
          write_to_file( join(id_list, '\n'), id_list_file )
          #searching for the go associated to these ids
          id_go_associations_hash={} #key: id; value: list of GO code strings e.g, "GO:0055114"
          id_go_associations_string= bbash("gawk -v id_file="+id_list_file+""" -F"\\t" 'BEGIN{ while ((getline idline < id_file)>0){ GI_INPUT[idline]=1 } } { split($1, GI, "; "); split("", gi_match); for (i=1; i<=length(GI); i++){ if (GI[i] in GI_INPUT) { gi_match[GI[i]]=1}   }; o=""; for (g in gi_match) o=o"; " g; if (o)  print substr(o, 3) "\\t" $2  }' """+self.profile.uniref2go_db_filename(), dont_die=1)
          if id_go_associations_string:
            for line in id_go_associations_string.split('\n'):
              for id_code in line.split('\t')[0].split('; '): #more than one id can be present in a line
                try:
                  id_go_associations_hash[id_code]=line.split('\t')[1].split('; ')
                except:
                  printerr('WARNING go_score error parsing line: '+line, 1)
          total_score=0
          for name, evalue in names_and_evalues[:max_n_titles]:        
            id_code=name.split()[0]
            go_list_for_this_title=[]
            if not id_go_associations_hash.has_key(id_code):           score_this_title=0
            else:
              go_list_for_this_title=id_go_associations_hash[id_code]
              #checking if any of the GO terms for this family is found in the GOs for this blast title
              found_go=False;               any_molecular_function_annotated=False
              for go_for_this_title in go_list_for_this_title:
                try:
                  all_parents_of_go_for_this_title=  gene_ontology().get_all_parents_ids(go_for_this_title)
                  #checking if the GO terms are referred to molecular function
                  if 'GO:0003674' in all_parents_of_go_for_this_title:
                    any_molecular_function_annotated=True
                    if any( [t in all_parents_of_go_for_this_title+[go_for_this_title]  for t in terms]): found_go=True
                except NoSuchTermError:
                  printerr('WARNING GO found in annotated protein was not found among known GO terms: '+go_for_this_title+' . GO databases should be both updated to ensure consistency!', 1)
              if not any_molecular_function_annotated:        score_this_title=0
              elif found_go:                                  score_this_title=  shortcut_log(evalue)
              else:                                           score_this_title= -shortcut_log(evalue)
            if verbose: write(name+' GO_TERMS: '+str(go_list_for_this_title)+' --> '+str(score_this_title), 1)
            total_score+=score_this_title
        self.go_score_data=total_score
    return self.go_score_data
  
  def blast_against_tag_db(self, outfile='', silent=False, force=False):
    """ Utility to perform a 'tag blast' of a predicted sequence (p2ghit) against a reference database (typically nr); the output can then be used to determine the tag_score of the prediction. If outfile is not defined, it is set to: blast_nr_folder_profile_subfolder+  self.profile.name+'.'+self.id+'.blast_nr'  (e.g. selenoprofiles_results/A.anophagefferens/tab_blast/SelR/SelR.15.blast_nr).
    If silent is not True, a service message is displayed while blast is running
    Normally, the output file is tested, and if it is a valid output, blast is not run. 
    If force is True, blast is run anyway.
    Returns the path to the outfile produced
    """
    if not outfile:  outfile= Folder(blast_nr_folder+self.profile.name)+ self.profile.name+'.'+str(self.id)+'.blast_nr'
    if force or not is_valid_blast_output(outfile): #checking if the output file is already present
      query_file=temp_folder+self.output_id()+'.fasta';     write_to_file(self.fasta(), query_file)
      blast_db=self.profile.tag_db_filename()
      if not is_file(blast_db): raise notracebackException, "selenoprofiles ERROR nr database file not found: "+blast_db+'  -- This is necessary to use the tag_score and go_score methods. Please change the .config file of this profile and do not use these methods, or download the nr database and make the tag_db.DEFAULT variable in the main configuration file point to it.'
      tag_blast_options= self.profile.tag_blast_options_dict();     tag_blast_options_string=join(['-'+k+' '+str(tag_blast_options[k]) for k in tag_blast_options], ' ')
      if not silent: service('Running tag blast for '+self.output_id()+' against db: '+blast_db+'  -> '+outfile)
      b=bash('which blastall')
      if b[0]: raise notracebackException, "ERROR blastall not found! Please install it"
      blastall_bin= dereference( b[1] )
      bbash(blastall_bin+' -I -p blastp -d '+blast_db+' -i '+query_file+' '+tag_blast_options_string+' > '+temp_folder+'tagblast_out')
      bbash('mv '+temp_folder+'tagblast_out '+outfile)      
      if not silent:     service('')
      bbash('rm '+query_file)
    return outfile
  tag_blast=blast_against_tag_db

  def three_prime(self, length=-1, chromosome_length=False):
    """ Returns the fasta text for the sequence downstream this p2ghit """
    if length==-1:
      if 'three_prime_length' in globals():        length=three_prime_length
      else: 
        try:     length=get_MMlib_var('three_prime_length')
        except:  raise Exception, "p2ghit -> three_prime ERROR length was not specified nor the variable three_prime_length was present in MMlib"
    ### getting a gene object with right boundaries, also checking chromosome length
    downstream_g=self.downstream(distance=0, region_length=length)
    if not chromosome_length: 
      if 'chromosome_lengths' in globals():        chromosome_length=chromosome_lengths[self.chromosome]
      else: 
        try:     chromosome_length= get_MMlib_var('chromosome_lengths')[self.chromosome]
        except:   raise Exception, "p2ghit -> three_prime ERROR chromosome_length was not specified nor the variable chromosome_lengths was present in MMlib"
    downstream_g.check_boundaries(chromosome_length)
    if not downstream_g: return None #out of boundaries
    ### getting sequence, possibly getting the one in memory already
    if length <= get_MMlib_var('three_prime_length'):
      title= downstream_g.header()
      seq=   self.sequence_of_internal_gene(downstream_g)
    else: 
      title, seq  = downstream_g.fasta_sequence()
    ## out title
    new_title   = 'region_of_'+ str(downstream_g.length())+'_bp_downstream_of:'+self.output_id()+' '+join(title.split()[1:-1], ' ')
    return ">"+new_title+'\n'+fasta(seq) 
  output_three_prime = three_prime 

  def five_prime(self, length=-1, chromosome_length=False):
    """ Returns the fasta text for the sequence upstream this p2ghit. option/config parameter five_prime_length must be defined """
    if length==-1:
      if 'five_prime_length' in globals():        length=five_prime_length
      else:
        try:     length=get_MMlib_var('five_prime_length')
        except:  raise Exception, "p2ghit -> five_prime ERROR length was not specified nor the variable five_prime_length was present in MMlib"
    ### getting a gene object with right boundaries, also checking chromosome length
    upstream_g=self.upstream(distance=0, region_length=length)
    if not chromosome_length:
      if 'chromosome_lengths' in globals():        chromosome_length=chromosome_lengths[self.chromosome]
      else:
        try:     chromosome_length= get_MMlib_var('chromosome_lengths')[self.chromosome]
        except:   raise Exception, "p2ghit -> five_prime ERROR chromosome_length was not specified nor the variable chromosome_lengths was present in MMlib"
    upstream_g.check_boundaries(chromosome_length)
    if not upstream_g: return None #out of boundaries
    if length <= get_MMlib_var('five_prime_length'):    
      title= upstream_g.header()
      seq=   self.sequence_of_internal_gene(upstream_g)
    else: 
      title, seq  = upstream_g.fasta_sequence()
    new_title   = 'region_of_'+ str(upstream_g.length())+'_bp_upstream_of:'+self.output_id()+' '+join(title.split()[1:-1], ' ')
    return ">"+new_title+'\n'+fasta(seq)
  output_five_prime = five_prime
  
  def boundaries_in_profile(self):
    """ Returns the boundaries that the target spans in the profile [start, end], both 1 based and included"""
    if not self.__dict__.has_key('profile') or (not self.profile):  raise Exception, "profile_alignment->boundaries_in_profile ERROR no profile alignment is defined for this p2g hit. Please add this as its property "
    query_start, query_end =  self.query.boundaries() #all positions in this function are 1 based
    query_full_name=self.query_full_name()
    prediction_start_in_ali = self.profile.position_in_ali(query_full_name, query_start)
    prediction_end_in_ali =   self.profile.position_in_ali(query_full_name, query_end)    
    return [prediction_start_in_ali, prediction_end_in_ali]

  def is_contained_in_profile_range(self, pos_start, pos_end):
    """ This function tells if the prediction spans only a certain portion of the profile. Positions are 1 based. 
    This is thought to filter out prediction which spans uniquely a regino of the profile with low sequence information/repetitive sequence    """
    prediction_start, prediction_end= self.boundaries_in_profile() #positions are 1 based and referred to the original profile
    return prediction_start>=pos_start and prediction_end<=pos_end

  def spans_profile_range(self, pos_start, pos_end):
    """ This function tells if the prediction spans a certain portion of the profile. Positions are 1 based.    """
    prediction_start, prediction_end= self.boundaries_in_profile() #positions are 1 based and referred to the original profile
    return are_overlapping_ranges( [prediction_start, prediction_end], [pos_start, pos_end])

  def target_sequence_is_repetitive(self, minimal_aa_set_size=5):
    """ This function is useful to filter out spurious hits. It checks if the target is just the repetition of a few aminoacids:
    if the set of distinct aminoacids contained in the sequence predicted for tha target is < than minimal_aa_set_size (def:5), it returns True)    """
    return len( all_chars_in(self.alignment.seq_of('t')) ) < minimal_aa_set_size

  def alignment_is_repetitive(self, minimal_aa_set_size=5, only_identities=False):
    """ This function (similar to target_sequence_is_repetitive) considers the target sequence in the positions in which it is conserved, defined as the positions in which query and target carries similar aminoacids (defined be similar_aas function in MMlib): the aminoacids in these positions are concatenated to form an "alignment core". If this is formed only by few aminoacids (less than minimal_aa_set_size, def:5), True is returned, otherwise False is returned.  The function is thought for filtering out spurious hits. """
    if only_identities:    
      all_conserved_positions=[]; id_matrix=self.alignment.identity_matrix()
      for pos in range(len( id_matrix )):     
        if id_matrix[pos]: all_conserved_positions.append(pos+1)
    else:                  all_conserved_positions=self.alignment.positions_of_similar_aas() #1 based!
    concatenated_seq_in_conserved_positions= join([ self.alignment.seq_of('t')[p-1] for p in all_conserved_positions   ], '')
    return len( all_chars_in(concatenated_seq_in_conserved_positions) ) < minimal_aa_set_size
  
  def show_conservation_in_profile_range(self, profile_start, profile_end, min_aas=4):
    """ This function is useful to filter out spurious hits. It returns true only if the prediction has at least min_aas (def:4) which are conserved in the query in the indicated profile range. Conserved here means they must be similar amino acids according to similar_aas function. It is a more strict version of  is_contained_in_profile_range  (but reversed in boolean value). positions in input are 1 based """
    prediction_start, prediction_end= self.boundaries_in_profile()
    if  prediction_start>=profile_end or prediction_end<=profile_start : return False
    query_full_name=self.profile.fill_title( self.query.chromosome )
    start_pos_in_query = self.profile.position_in_seq(query_full_name, profile_start)
    end_pos_in_query = self.profile.position_in_seq(query_full_name, profile_end)
    # the function position_in_seq returns 0 when the sequences contains only gaps in the region upstream and including profile_start. In this case we want to set the initial position in the query to 1, unless the query spans a region not overlapping with profile_start --- profile_end
    if    start_pos_in_query==0 and end_pos_in_query==0:      return False    
    elif  start_pos_in_query==0: start_pos_in_query=1
    prediction_pairwise_ali= self.filled_alignment(target_title='t')
    start_pos_in_prediction_pairwise_ali =   prediction_pairwise_ali.position_in_ali(query_full_name, start_pos_in_query)
    end_pos_in_prediction_pairwise_ali   =   prediction_pairwise_ali.position_in_ali(query_full_name, end_pos_in_query)
    number_of_similar_aas=0
    for p in range(start_pos_in_prediction_pairwise_ali, end_pos_in_prediction_pairwise_ali+1):
      aa_in_query= prediction_pairwise_ali.seq_of(query_full_name)[p-1]
      aa_in_target=prediction_pairwise_ali.seq_of('t')[p-1]
      if aa_in_query==aa_in_target or similar_aas( aa_in_query ,   aa_in_target  ):      number_of_similar_aas+=1
    return number_of_similar_aas>= min_aas

  def has_redox_box(self, sec_index=1):
    """ This returns truth if the prediction show a redox box including the sec_position in the profile. If more than one sec is present, sec_index can be used to specify which position. 1 stands for first"""
    sec_pos= self.profile.sec_pos()[sec_index-1]+1   #sec pos is 1 based
    p_start, p_end = self.boundaries_in_profile()
    if sec_pos >= p_start and sec_pos <= p_end:
      return self.seq_in_profile_pos(sec_pos) in 'UC' and (self.seq_in_profile_pos(sec_pos-3) in 'UC' or self.seq_in_profile_pos(sec_pos+3) in 'UC' ) 
    return False
    

  def pretty_alignment(self, chars_per_line=100):
    """ Shows the pairwise alignemnt between query and target"""
    o=''
    target_cds=self.aligned_cds()
    a=''; b=''; c=''; d=''; e=''; f=''; U_positions_in_out_strings=[]; stops_positions_in_out_strings=[]
    exon_index=0;     bases_walked_in_this_exon=0;  nt_walked_in_target=0; frameshift_index=0; intron_index=0
    nt_in_other_exon=-1 # -1 means: we're not in a intron boundary
    splice_sites_seq= self.splice_site_sequences()
    frameshifts_seq = self.frameshift_sequences()
    for pos in range(self.alignment.length()): #0 based
      if self.strand=='+':    codon_start=self.exons[exon_index][0]+bases_walked_in_this_exon; codon_end=codon_start+2
      elif self.strand=='-':  codon_start=self.exons[exon_index][1]-bases_walked_in_this_exon; codon_end=codon_start-2
      if   self.strand=='+':  nt_in_other_exon = codon_end-self.exons[exon_index][1]     # can be 0, 1, 2
      elif self.strand=='-':  nt_in_other_exon = self.exons[exon_index][0]-codon_end     # can be 0, 1, 2
      aa_query= self.alignment.seq_of(0)[pos]
      aa_target=self.alignment.seq_of(1)[pos]
      if aa_query==  'U' or aa_target=='U': U_positions_in_out_strings.append( len(a) )
      if aa_target in 'u*': stops_positions_in_out_strings.append(len(a))
      if aa_query== aa_target:                    b_char='|'
      elif similar_aas(aa_query, aa_target):      b_char='/'
      else:                                       b_char=' '
      a+= aa_query
      b+= b_char
      c+= aa_target 
      d+= lower(target_cds[pos*3])
      if nt_in_other_exon<2 :   e+= lower(target_cds[pos*3+1])
      else:                     e+=' '
      if nt_in_other_exon<1 :   f+= lower(target_cds[pos*3+2])
      else:                     f+=' '

      if aa_target!='-':
        nt_walked_in_target+=3
        bases_walked_in_this_exon += 3
        if self.strand=='+':      
          next_codon_start=self.exons[exon_index][0]+bases_walked_in_this_exon;
          if next_codon_start>  self.exons[exon_index][1]: #we are beyond this exon... let's go to the next one. Drawing intron
            if nt_in_other_exon in [1,2]:              bases_walked_in_this_exon=nt_in_other_exon
            else:                                      bases_walked_in_this_exon=0
            exon_index+=1
            if exon_index<len(self.exons):
              intron_length=self.exons[exon_index][0]-self.exons[exon_index-1][1]-1
              if intron_length>5:
                a+=' <---Intron---> '
                b+=' <'+(str(intron_length)+'nt').center(12)+'> '
                c+='                '
                d+='                '
                try:                   eseq='   '+lower(  splice_sites_seq[intron_index][:2]  ) +'     '+lower(splice_sites_seq[intron_index][-2:])+'   '
                except:                eseq='               '
                if nt_in_other_exon==2: e+=eseq + lower(target_cds[pos*3+1])
                else:                   e+=eseq+' '
                if nt_in_other_exon>=1: f+='               '+ lower(target_cds[pos*3+2])
                else:                   f+='                '                         
                intron_index+=1
              else: #short intron is actually a frameshift!
                a+=' ! FRAME ! '
                b+=' ! SHIFT ! '
                c+=' '+(str(intron_length)+'nt').center(9) +' '
                d+='           '
                try:  eseq=lower( frameshifts_seq[frameshift_index]  ).center(11)
                except:   eseq='           '
                e+=eseq
                f+='           '
                frameshift_index+=1
            nt_in_other_exon=-1           
        elif self.strand=='-':    
          next_codon_start=self.exons[exon_index][1]-bases_walked_in_this_exon;        
          if next_codon_start<self.exons[exon_index][0]: #we are beyond this exon... let's go to the next one. Drawing intron
            if nt_in_other_exon in [1,2]:              bases_walked_in_this_exon=nt_in_other_exon
            else:                                      bases_walked_in_this_exon=0
            exon_index+=1
            if exon_index<len(self.exons):
              intron_length=self.exons[exon_index-1][0]-self.exons[exon_index][1]-1
              if intron_length>5:
                a+=' <---Intron---> '
                b+=' <'+(str(intron_length)+'nt').center(12)+'> '
                c+='                '
                d+='                '
                try:                   eseq='   '+lower(  splice_sites_seq[intron_index][:2]  ) +'     '+lower(splice_sites_seq[intron_index][-2:])+'   '
                except:                eseq='               '
                if nt_in_other_exon==2: e+=eseq + lower(target_cds[pos*3+1])
                else:                   e+=eseq+' '
                if nt_in_other_exon>=1: f+='               '+ lower(target_cds[pos*3+2])
                else:                   f+='                '                         
                intron_index+=1
              else: #short intron is actually a frameshift!
                a+=' ! FRAME ! '
                b+=' ! SHIFT ! '
                c+=' '+(str(intron_length)+'nt').center(9) +' '
                d+='           '
                try:  eseq=lower( frameshifts_seq[frameshift_index]  ).center(11)
                except:   eseq='           '
                e+=eseq
                f+='           '
                frameshift_index+=1
            nt_in_other_exon=-1           
    g=' '*len(a)
    for p in U_positions_in_out_strings:          g=g[:p]+'*'+g[p+1:]
    for p in stops_positions_in_out_strings:      g=g[:p]+'X'+g[p+1:]    
    #splitting to tot chars_per_line       
    for index_returns in range( 1+(len(a)-1)/chars_per_line ):
      o+= '\nQuery   '+a[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'
      o+= '        '+b[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'        
      o+= 'Target  '+c[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'
      o+= '        '+d[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'
      o+= '        '+e[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'        
      o+= '        '+f[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'
      o+= '        '+g[  chars_per_line*(index_returns) :  chars_per_line*(index_returns+1) ]+'\n'
    return o


  def output_p2g(self, chars_per_line=100):
    """ Selenoprofiles custom output. It contains the protein alignment and the underlining sequence of the target, in a format similar to genewise. All introns shorter than 6 bp are shown as frameshifts. """
    o=''
    try:    
      assert issubclass(self.profile.__class__, alignment)
      assert self.profile
    except: raise Exception, "ERROR output_p2g: the profile attribute must be set to a profile alignment"
    try: o+="Output_id:  "+self.output_id()+'\n'
    except: pass
    try: o+="----------  "+"-"*len(self.output_id())+"\n"
    except: pass
    try: o+=("-Species        "+str(self.species)).ljust(60)+' -Taxid '+str(self.species.taxid)+'\n'
    except: pass
    try: o+="-Target         "+str(self.target) +'\n'
    except: pass
    o+="-Chromosome ("+self.strand+") "+str(self.chromosome)+'\n'
    o+="-Program        "+self.prediction_program() +'\n'  
    try: o+="-Query name     "+self.query_full_name() +'\n'  
    except: pass
    q_bounds=self.query.boundaries(); q_length=len(nogap(self.profile.seq_of(self.query_full_name())));    q_bonds_in_profile=self.boundaries_in_profile()
    sec_pos_in_profile_1_based= [i+1 for i in self.profile.sec_pos() ]
    try:
      o+="-Query range    "+join([str(i) for i in q_bounds], '-').ljust(10)           +' length:'+str(q_length).ljust(5)+' coverage: '+str(round(float((q_bounds[1]-q_bounds[0]+1))/q_length, 2)) +'\n'
      o+="-Profile range  "+join([str(i) for i in q_bonds_in_profile], '-').ljust(10) +' length:'+str(self.profile.length()).ljust(5)+' coverage: '+str(round(self.coverage(), 2))
      if sec_pos_in_profile_1_based:       o+='    sec_position'+'s'*int(len(sec_pos_in_profile_1_based)>1)+': '+str(sec_pos_in_profile_1_based)
      o+='\n'    
    except: pass
    try:  o+='-ASI:           '+str(round(self.sequence_identity_with_profile(), 4))+'     (ignoring gaps: '+str(round(self.sequence_identity_with_profile(True), 4))+')\n'
    except: pass
    try:  o+='-AWSIc:         '+str(round(self.weighted_seq_identity_with_profile(with_coverage=True), 4))+ '     Z-score: '+str(round(self.z_score(with_coverage=True), 3))+'\n'
    except: raise
    try:  o+='-AWSIw:         '+str(round(self.weighted_seq_identity_with_profile(with_coverage=False), 4))+'     Z-score: '+str(round(self.z_score(with_coverage=False), 3))+'\n'
    except: raise  
    try:  o+="-State          "+str(self.filtered)
    except: pass
    o+='\n\n'    
    o+=" alignment ".center(25, '-') #+"\n" #not putting the \n, it's on the Query ali line
    o+=self.pretty_alignment(chars_per_line=chars_per_line)
    o+=' positions '.center(25, '-')+'\n'
    for exon_index in range(len( self.exons )):
      p_start, p_end= self.exons[exon_index]
      o+=('Exon '+str(exon_index+1)).ljust(10)+str(p_start).ljust(12)+' '+str(p_end)+'\n'
    o+='\n'
    o+=' features '.center(25, '-')+'\n'
    added_at_least_feature=False
    for f in self.features:
      if f.included_in_output:         
        o+=f.output()+'\n'
        added_at_least_feature=True
    if not added_at_least_feature:        o+='None\n'
    try:
      o+=" 3' seq ".center(25, '-')+'\n'
      o+='Total sequence length available downstream '
      seq_until_stop, translation, length_of_three_prime = self.three_prime_coding_tail()
      if length_of_three_prime >=get_MMlib_var('three_prime_length'):     o+='>= '+str(get_MMlib_var('three_prime_length'))+'\n'
      else:                                o+='= '+str(length_of_three_prime)+'\n'
      o+='Sequence until first stop codon: \n'+seq_until_stop+'\n'
      o+=' '+join(list(translation), '  ')+' \n'
    except:
      o+="Cannot get the sequence downstream! "
    return o

  def get_misc_sequence_data(self):
    m=[f for f in self.features if f.__class__ == misc_sequence_data]
    if len(m): 
      if len(m)>1: raise Exception, "ERROR two misc_sequence_data_feature found in object: "+self.id
      return m[0]
    else: return self.compute_misc_sequence_data()

  def remove_misc_sequence_data(self):
    """ Remove sequence data from memory. useful if the coordinates of the prediction changed in anyway. Use compute_misc_sequence_data to recompute the data """
    for index_to_remove in [ index     for index, f in enumerate(self.features) if f.__class__ == misc_sequence_data ][::-1]:        self.features.pop(index_to_remove)

  def compute_misc_sequence_data(self):
    """ This function loads some sequence data for lazy computing, and stores it into a genomic feature of this object. calling it is equivalent to force recomputing it, if data is already found"""
    ## removing misc_sequence_data object if already present.
    self.remove_misc_sequence_data()
    ## preparing to fetch just a single time all necessary sequence
    chromosome_length= get_MMlib_var('chromosome_lengths')[self.chromosome]
    g=misc_sequence_data(chromosome=self.chromosome, strand=self.strand);    #minimal copy of self object, like a gene_boundaries
    self_boundaries=self.boundaries(); g.exons=[ [self_boundaries[0], self_boundaries[1]] ]; length_before_extension= g.length()
    if    self.strand=='+':     g.extend( right= get_MMlib_var('three_prime_length'), left= get_MMlib_var('five_prime_length'),   inplace=True )
    elif  self.strand=='-':     g.extend( left=  get_MMlib_var('three_prime_length'), right=get_MMlib_var('five_prime_length'),   inplace=True )
    g.check_boundaries(chromosome_length)
    g.full_seq= upper( g.fasta_sequence()[1] )
    if   self.strand=='+':     
      g.extension_upstream=     self.boundaries()[0] - g.boundaries()[0]
      g.extension_downstream=   g.boundaries()[1]    - self.boundaries()[1]
    elif self.strand=='-':     
      g.extension_upstream=     g.boundaries()[1]    - self.boundaries()[1] 
      g.extension_downstream=   self.boundaries()[0] - g.boundaries()[0]
    g.parent=self
    self.features.append(g)
    return self.features[-1]

  def subsequence(self, start, length, include_introns=False):
    """ Lazy computed subsequence. It should be used to get the sequence of the gene or any subportion, such as CDS, three prime etc. 
    Start is 1-based, and everything is nucleotide based. So for example, .subsequence(1, 3) reports the first codon in the prediction.  
    Negative start positions can be used to get the five prime sequence. start==0 is the nt right before the first coding position, etx. So for example, .subsequence(-99, 100) reports 100 nt upstream of the p2g. Large values of length will result in getting the sequence at 3' UTR. 
    IMPORTANT NOTE: the function will crash raising an exception if the sequence which you're trying to fetch is not available in the misc_sequence_data, which is controlled by three_prime_length and five_prime_length parameters. So if five_prime_length is 0, you cannot use this function to get the upstream sequence.
  Only the coding sequence is considered normally, if include_introns is active, the positions are relative to the whole gene structures, from the first coding letter to the last. """    
    subseq_upstream='';    subseq_downstream='';   subseq_middle=''
    subseq_length=length;  subseq_start=start
    g = self 
    if include_introns:         g= self.boundaries_gene(minimal=True)
    if subseq_start<1:
      subseq_upstream=   self.sequence_of_internal_gene(   self.upstream(0, -subseq_start+1)   )
      if len( subseq_upstream ) > subseq_length:  subseq_upstream= subseq_upstream[:subseq_length]        
      subseq_length=    max(0, subseq_length - len(subseq_upstream) ) 
      subseq_start=1
    max_subseq_length=  g.length() - subseq_start +1
    if subseq_length > max_subseq_length:
      subseq_downstream= self.sequence_of_internal_gene(   self.downstream(0, subseq_length-max_subseq_length )   )
      subseq_length= max_subseq_length
    if subseq_length>0:
      s=g.subseq(subseq_start, subseq_length, minimal=True)
      subseq_middle=self.sequence_of_internal_gene(s)
    subseq =   subseq_upstream+      subseq_middle   +subseq_downstream
    return  subseq
    
  def sequence_of_internal_gene(self, s, trying_recomputing_misc_data=False):
    """ Get the sequence for a gene object s, that must be included within the boundaries of the misc_sequence_data available """
    g=self.get_misc_sequence_data()
    if s.boundaries()[0] < g.boundaries()[0] or s.boundaries()[1] > g.boundaries()[1]: 
      if not trying_recomputing_misc_data:
        self.compute_misc_sequence_data()
        return self.sequence_of_internal_gene(s, trying_recomputing_misc_data=True)
      else:  #we're trying for the second time already, no way... 
        raise lazycomputationException, "ERROR sequence_of_internal_gene  the gene is not internal to the available misc_sequence_data! you may try extending three_prime_length / five_prime_length" 
    rel_coord = s.relative_coordinates(    *g.boundaries() , reverse=self.strand=='-'   )  # also gene obj
    out_seq=''
    for exon_start, exon_end in rel_coord.exons:      out_seq+= g.full_seq[ exon_start-1: exon_end ]     
    return out_seq

  def splice_site_sequences(self):
    """ Returns the splice sites for this prediction. A list of elements is reported, one for each intron (bigger than 5 bp), each element is 4 letters long and is composed by the first 2 bp of the intron plus the last 2bp."""
    out_list=[]
    introns_g= self.introns(minimal=True)    
    for intron_index in range(len(introns_g.exons)):
      intron_g= introns_g.get_exon(intron_index, minimal=True)
      if intron_g.length()>5:            #real intron
        intron_sequence= self.sequence_of_internal_gene(  intron_g  )
        out_list.append(  intron_sequence[:2]+intron_sequence[-2:] )
    return out_list

  def frameshift_sequences(self):    
    """ Returns the sequences of the insertions causing frameshifts for this prediction. A list of sequence elements is reported, one for each "intron" smaller than 6 bp."""
    out_list=[]
    introns_g= self.introns(minimal=True)    
    ### splice sites and frameshifts
    for intron_index in range(len(introns_g.exons)):
      intron_g= introns_g.get_exon(intron_index, minimal=True)
      if intron_g.length()<=5:            #frameshift
        f_sequence= self.sequence_of_internal_gene(  intron_g  )
        out_list.append(  f_sequence )
    return out_list
      
  def three_prime_coding_tail(self):
    """ Returns [nt_seq, aa_translation, length_of_three_prime] of the three_prime sequence until the first stop codon included.  For a complete prediction, it would return a stop codon and an "*" as translation. """
    t=self.three_prime()
    if t is None: return [ None, None, 0 ]
    seq=upper(join(t.split('\n')[1:], ''))    
    seq_until_stop=''; index_codon=0; codon='XXXXX'
    while codon:
      codon= seq[index_codon*3:index_codon*3+3]
      seq_until_stop+=codon
      if codon in ['TGA', 'TAG', 'TAA']: codon=''
      index_codon+=1                    
    translation=transl(seq_until_stop)
    return [seq_until_stop, translation, len(seq)]

  def remove_terminal_uga(self, distance_from_end_of_prediction=0):
    """ This function serves to correct a unwanted behaviour which is a consequences of the scoring schemes used for selenoproteins. Genewise (I don't yet if exonerate does it as well) sometimes includes a terminal UGA in the prediction, aligned to a amino acid which is not U in the query. This function detects and remove this codon from the prediction. Only one codon can be removed. The variable distance_from_end_of_prediction determines the maximum distance of the removed UGA from the end of coding sequence predicted (0= can be only the last codon of CDS)
    
    GPGGKTFVSLLDMPRPFTKLKALDVEALAEEVLAALKE
    |||  /||/ |||||  |/ ||//  | /|/ ||| 
    VHGGKEVLSLVGMPRPFKSLRELDMDEAATQVVEALK*
    gcggaggctcggaccctaaccgcgaggggacggggtat
    taggaattcttgtcgctagtgatataacccattactag
    ttacagcacagtgttgcatttagtgtgtttagggagga

    ----> 

    GPGGKTFVSLLDMPRPFTKLKALDVEALAEEVLAALK
    |||  /||/ |||||  |/ ||//  | /|/ ||| 
    VHGGKEVLSLVGMPRPFKSLRELDMDEAATQVVEALK
    gcggaggctcggaccctaaccgcgaggggacggggta
    taggaattcttgtcgctagtgatataacccattacta
    ttacagcacagtgttgcatttagtgtgtttagggagg     """

    for index_ali in range( self.alignment.length()-1, self.alignment.length()-2-distance_from_end_of_prediction, -1):
      aa_query =self.alignment.seq_of('q')[index_ali]
      aa_target=self.alignment.seq_of('t')[index_ali]      
      if aa_query!='U' and aa_target=='*' and upper(self.aligned_cds()[index_ali*3:index_ali*3+3])=='TGA':
        #trying to extend the region to remove, checking for gaps in query or target right before the TGA
        while self.alignment.seq_of('q')[index_ali-1] =='-' or self.alignment.seq_of('t')[index_ali-1] =='-': index_ali-=1 #NB increasing index_ali! this works only if we never go to the next element in the for loop
        query_seq_to_remove         = self.alignment.seq_of('q')[index_ali:]
        number_of_query_aas_removed = len(query_seq_to_remove) - query_seq_to_remove.count('-')
        target_seq_to_remove         = self.alignment.seq_of('t')[index_ali:]
        number_of_target_aas_removed = len(target_seq_to_remove) - target_seq_to_remove.count('-')
        #correcting alignment
        self.alignment.set_sequence('q',  self.alignment.seq_of('q')[:index_ali] )
        self.alignment.set_sequence('t',  self.alignment.seq_of('t')[:index_ali] )
        #correcting query positions
        query_start, query_end= self.query.exons.pop(-1)
        query_end=  query_end - number_of_query_aas_removed
        self.query.add_exon(query_start, query_end)
        #correcting target positions
        target_start, target_end= self.exons.pop(-1)
        if self.strand=='+':  target_end=   target_end   - number_of_target_aas_removed*3
        if self.strand=='-':  target_start= target_start + number_of_target_aas_removed*3
        self.add_exon(target_start, target_end)        
        return number_of_target_aas_removed           #returning a value while inside the for loop  =  breaking it on the first aa to remove

  def join_adjacent_exons(self):
    """This method is necessary because very rarely when you join exonerate predictions you end up with two exons which are adjacent and this may screw the functioning of methods such as .introns(). This full join such exons in a single exon (so we don't leave a 0 bp intron). Also, it checks wheter the joined exons are too close together and out of frame, which  imply frameshifts are present."""
    for exon_index in range(1, len(self.exons)):
      if self.strand=='+':
        intron_length=self.exons[exon_index][0]-self.exons[exon_index-1][1]-1
        if intron_length==0:
          self.exons[exon_index-1][1]=self.exons[exon_index][1]
          self.remove_exon(exon_index)
          self.join_adjacent_exons()
          return
        elif intron_length<6 and intron_length!=3:
          position_in_the_cds_1_based= sum( [ end-st+1 for st, end in self.exons[:exon_index] ]  )            
      if self.strand=='-':
        intron_length=self.exons[exon_index-1][0]-self.exons[exon_index][1]-1
        if intron_length==0:
          self.exons[exon_index-1][0]=self.exons[exon_index][0]
          self.remove_exon(exon_index)
          self.join_adjacent_exons()
          return
        elif intron_length<6 and intron_length!=3:
          position_in_the_cds_1_based= sum( [ end-st+1 for st, end in self.exons[:exon_index] ]  )

  def clean_inframe_stop_codons(self, max_codons_removed=10, silent=False):
    """This methods analyze the prediction and, if there's stop codons in the target (predicted Sec, which are coded as Us, are not considered here), it tries to clean it from the prediction: so if they are close to a intron boundary, or to the start or end of prediction, the prediction is cut to remove the stop codons. """
    n_stop_codons_removed=0
    if self.frameshifts():         return ## function gives problem wwhen prediction has frameshifts
    if "*" in self.alignment.seq_of('t'): 
      old_p2g=self.copy()
      try:
        n_gaps_in_target=0
        for pos in range( self.alignment.length() ):
          if self.alignment.seq_of('t')[pos]=='-':          n_gaps_in_target+=1
          elif self.alignment.seq_of('t')[pos]=='*':
            this_aa_codon=self.subseq( 3*(pos-n_gaps_in_target)+1, 3)
          #now searching the nearest intron 
            introns_index_and_distances_and_direction=[] #each element is: [intron_index, distance_to_this_aa_codon, "upstream"or"downstream"]
            introns=self.introns()
            # adding pseudo introns standing for the non coding region at the two extremes of the gene.
            introns.add_exon(0, self.boundaries()[0]-1)  
            introns.add_exon(self.boundaries()[1]+1, sys.maxint)
          
            for intron_index in range(len(introns.exons)):
              if not introns.exons[intron_index][1] - introns.exons[intron_index][0] +1 <= 5: #excluding frameshifts
              
                intron=introns.get_exon(intron_index)
                d= intron.is_downstream_of(this_aa_codon)
                u= intron.is_upstream_of(this_aa_codon)
                b= this_aa_codon.boundaries_gene().overlaps_with(intron, phase=False )   # the stop codon maybe right on the intron boundary  

                if   b: introns_index_and_distances_and_direction.append([intron_index, -1, 'boundary']) #-1 is to ensure it will be first after sorting
                elif d: introns_index_and_distances_and_direction.append([intron_index, d-1, 'downstream'])
                elif u: introns_index_and_distances_and_direction.append([intron_index, u-1, 'upstream'  ])
          
            introns_index_and_distances_and_direction.sort(key=lambda x:x[1]) #sorting: first value refers now to closest intron
            nearest_intron_index, distance, direction =  introns_index_and_distances_and_direction[0]

            if distance <= max_codons_removed*3: 
              n_stop_codons_removed+=1
              if direction =='downstream':
                gene_portion_to_remove = self.subseq( 3*(pos-n_gaps_in_target)+1, distance+3)
                n_aminoacids_to_remove = distance / 3
                n_nucleotides_left_at_intron_boundary = distance % 3
                if n_nucleotides_left_at_intron_boundary:  #in this case there must be another exon downstream which carries the rest of the codon, which we want to remove. 
                  n_aminoacids_to_remove+=1
                
                  exon_to_remodel=self.exons[nearest_intron_index]
                  self.remove_exon( nearest_intron_index )  
                  if self.strand=='+':    exon_to_remodel[0]= exon_to_remodel[0]+(3-n_nucleotides_left_at_intron_boundary)
                  elif self.strand=='-':  exon_to_remodel[1]= exon_to_remodel[1]-(3-n_nucleotides_left_at_intron_boundary)
                  if exon_to_remodel[1]>= exon_to_remodel[0]: # to avoid bug when we have a single nt exon (which is a prediction bug by itself, but anyway)
                    self.add_exon(exon_to_remodel[0], exon_to_remodel[1])
                positions_in_alignment_to_remove=0; n_gaps=0   #not counting the positions of the stop
              
                while positions_in_alignment_to_remove != n_aminoacids_to_remove+n_gaps:
                  positions_in_alignment_to_remove+=1
                  if self.alignment.seq_of('t')[pos+positions_in_alignment_to_remove]=='-': n_gaps+=1
                
              elif direction =='upstream':
                gene_portion_to_remove = self.subseq( 3*(pos-n_gaps_in_target)+1-distance, distance+3)
                n_aminoacids_to_remove = distance / 3
                n_nucleotides_left_at_intron_boundary = distance % 3
                if n_nucleotides_left_at_intron_boundary:  #in this case there must be another exon downstream which carries the rest of the codon, which we want to remove. 
                  n_aminoacids_to_remove+=1
                  exon_to_remodel=self.exons[nearest_intron_index-1]
                  self.remove_exon( nearest_intron_index-1 )
                  if self.strand=='+':    exon_to_remodel[1]= exon_to_remodel[1]-(3-n_nucleotides_left_at_intron_boundary)
                  elif self.strand=='-':  exon_to_remodel[0]= exon_to_remodel[0]+(3-n_nucleotides_left_at_intron_boundary)
                
                  if exon_to_remodel[1]>= exon_to_remodel[0]: # to avoid bug when we have a single nt exon (which is a prediction bug by itself, but anyway)
                    self.add_exon(exon_to_remodel[0], exon_to_remodel[1])

                positions_in_alignment_to_remove=0; n_gaps=0   #not counting the positions of the stop
                while positions_in_alignment_to_remove != n_aminoacids_to_remove+n_gaps:
                  positions_in_alignment_to_remove+=1
                  if self.alignment.seq_of('t')[pos-positions_in_alignment_to_remove]=='-': n_gaps+=1

              elif direction =='boundary':
                gene_portion_to_remove = this_aa_codon  # --> this will do all the work
                n_aminoacids_to_remove = 0
                positions_in_alignment_to_remove =0
              #NOTE! positions_in_alignment_to_remove does not include the stop codon being cut off

              a= self.subtracted_of( gene_portion_to_remove  )
              if direction in ['downstream', 'boundary']: a.alignment.set_sequence('t',    a.alignment.seq_of('t')[:pos] + '-'*(1+positions_in_alignment_to_remove) + a.alignment.seq_of('t')[pos+positions_in_alignment_to_remove+1:])  
              elif direction=='upstream':                 a.alignment.set_sequence('t',    a.alignment.seq_of('t')[:pos-positions_in_alignment_to_remove] + '-'*(1+positions_in_alignment_to_remove) + a.alignment.seq_of('t')[pos+1:])    
              self.__dict__ = a.__dict__
              self.reset_derived_data()
              if not silent: write((self.profile.name+'.'+self.id+' ['+self.prediction_program()+'] ').ljust(22) +' stop codon removed from CDS, '+str(n_aminoacids_to_remove+1)+' codon'+('s'*bool(n_aminoacids_to_remove>0) )+' removed from the alignment', 1)

              self.alignment.remove_useless_gaps()
              self.clean_unaligned_tails()
              return self.clean_inframe_stop_codons(max_codons_removed=max_codons_removed, silent=silent) #instead of continuing the cycle through the positions, since it is difficult to take into account the reduction of the alignment, we just rerun the same function to see if it is necessary to cut off more codons from the prediction
      except:
        self.__dict__ = old_p2g.__dict__
        printerr('clean_inframe_stop_codons WARNING can\'t process '+self.profile.name+'.'+self.id+' ('+self.prediction_program()+') : skipping !', 1 )
        #if opt['debug']: raise

  def clean_unaligned_tails(self):
    """This function checks if there are unaligned portion of query or target at the N or C terminal and clean the prediction by removing them. """
    removed_something=False
    index=0 #end of n-.terminal tail of unaliged query
    while self.alignment.seq_of('q')[index] !='-' and self.alignment.seq_of('t')[index]=='-':          index+=1
    if index:
      self.query.exons[0][0]+=index
      self.alignment.set_sequence('q', self.alignment.seq_of('q')[index:])
      self.alignment.set_sequence('t', self.alignment.seq_of('t')[index:])
      removed_something=True
    index=-1 #index of c.-terminal tail. -2 means there's a 1 amionoacid long c-terminal tail
    while self.alignment.seq_of('q')[index]!='-' and self.alignment.seq_of('t')[index]=='-': index-=1
    if index!=-1:
      self.query.exons[0][1]+=(index+1)
      self.alignment.set_sequence('q', self.alignment.seq_of('q')[:index+1])
      self.alignment.set_sequence('t', self.alignment.seq_of('t')[:index+1])
      removed_something=True
    #checking N-terminal unaligned target
    index=0
    while self.alignment.seq_of('q')[index]=='-' and self.alignment.seq_of('t')[index]!='-':        index+=1
    if index: 
      gene_portion_to_remove=self.subseq(1, index*3)
      a=self.subtracted_of( gene_portion_to_remove  )
      a.alignment.set_sequence('q', self.alignment.seq_of('q')[index:])
      a.alignment.set_sequence('t', self.alignment.seq_of('t')[index:])
      self.__dict__ = a.__dict__
      self.reset_derived_data()
      removed_something=True
    #checking C-terminal unaligned target
    index=-1 
    while self.alignment.seq_of('q')[index]=='-' and self.alignment.seq_of('t')[index]!='-': index-=1
    if index!=-1:
      gene_portion_to_remove=self.reverse_subseq(1, -(index+1)*3 ) ###         -(index+1) = length of tail
      a=self.subtracted_of( gene_portion_to_remove  )
      a.alignment.set_sequence('q', self.alignment.seq_of('q')[:index+1])
      a.alignment.set_sequence('t', self.alignment.seq_of('t')[:index+1])
      self.__dict__ = a.__dict__
      self.reset_derived_data()
      removed_something=True
    #rerunning the function in case we have removed something, since our procedure may have left some new tails
    if removed_something: return self.clean_unaligned_tails()

  def exclude_large_introns(self, max_intron_length=140000, silent=False):
    """ Examine the prediction for the presence of large introns. If any is present, remove them by keeping the side with the largest codin sequence predicted """
    introns=self.introns(minimal=True)
    original_cols_alignment=self.alignment.length()
    for intron_index, (intron_start, intron_end) in enumerate( introns.exons ):
      intron_length= intron_end - intron_start +1
      if intron_length  > max_intron_length:
        ### removing this intron!
       
        length_left=  sum([ end-start+1 for index, (start, end) in enumerate(self.exons) if  index<= intron_index   ])   # total nt in cds prediction on the left  of this intron
        length_right= sum([ end-start+1 for index, (start, end) in enumerate(self.exons) if  index> intron_index   ])    # total nt in cds prediction on the right of this intron
        if length_right >= length_left :
          ### removing alignment on the left of this intron
          parsed_positions_target=0; parsed_positions_query=0; 
          for ali_index in range( self.alignment.length() ):
            seq_target= self.alignment.seq_of('t')[ali_index]
            seq_query=  self.alignment.seq_of('q')[ali_index]
            if seq_target!='-':   parsed_positions_target+=3
            if seq_query!='-':    parsed_positions_query+=1
            if parsed_positions_target>= length_left:
              break
          # now ali_index is the last position to drop
          ## modifying self.alignment
          self.alignment.columns(ali_index+1,    self.alignment.length() - ali_index-1,    inplace=True)
          ## modifying target coordinates in self.exons
          g=self.subseq(  parsed_positions_target+1, self.length()-parsed_positions_target, minimal=True ) ## using parsed_positions_target because it's %3  =  0 
          self.exons= g.exons
          ## modifying query coordinates in self.query.exons
          g=self.query.subseq(  parsed_positions_query+1, self.query.length()-parsed_positions_query , minimal=True ) ## using parsed_positions_query because it's %3  =  0 
          self.query.exons= g.exons

        elif length_left > length_right:
          ### removing alignment on the right of this intron
          parsed_positions_target=0; parsed_positions_query=0
          for ali_index in range( self.alignment.length() )[::-1]:  ## parsing backward 
            seq_target= self.alignment.seq_of('t')[ali_index]
            seq_query=  self.alignment.seq_of('q')[ali_index]
            if seq_target!='-':   parsed_positions_target+=3
            if seq_query!='-':    parsed_positions_query+=1
            if parsed_positions_target>= length_right:
              break
          # now ali_index is the last position to drop
          ## modifying self.alignment
          self.alignment.columns(0, ali_index, inplace=True)
          ## modifying target coordinates in self.exons
          g=self.subseq(1, self.length()-parsed_positions_target, minimal=True) ## using parsed_positions_target because it's %3  =  0 
          self.exons= g.exons
          ## modifying query coordinates in self.query.exons
          g=self.query.subseq(  1, self.query.length()-parsed_positions_query, minimal=True  ) ## using parsed_positions_query because it's %3  =  0 
          self.query.exons= g.exons

        if not silent:   write((self.profile.name+'.'+self.id+' ['+self.prediction_program()+'] ').ljust(22) +' removed intron too large: '+str(intron_length) +' nt; '+str(original_cols_alignment - self.alignment.length())+' columns removed from pairwise alignment', 1)
        self.alignment.remove_useless_gaps()
        self.clean_unaligned_tails()
        self.reset_derived_data()
        return self.exclude_large_introns(max_intron_length=max_intron_length, silent=silent) #rerunning function to remove other large introns, if present. not finishing the loop though
        
  def complete_at_three_prime(self, max_extension=10, max_query_unaligned=30, stops=None, silent=False):
    """ function to modify inplace the prediction to try adding some sequence at the 3' to the coding sequence prediction, if the number of aminoacids that would be added is <= max_extension, and if there's no strange character (anything different from ACGT); also, there must be at the most max_query_left aminoacids unaligned on the right side of the query   """
    if stops is None: 
      codon_table=get_genetic_code_table()      
      stops=set([codon  for codon in codon_table if codon_table[codon]=='*'  ])
    q_bounds=self.query.boundaries(); q_length=len(nogap(self.profile.seq_of(self.query_full_name())));    
    if q_length - q_bounds[1] >max_query_unaligned: return 
    tp= self.three_prime(length=max_extension*3+3)
    if tp is None: return   #out of boundaries exception
    splt=tp.split('\n')
    title= splt[0]
    three_prime_seq= join(splt[1:], '')
    for codon_index in range(0, len(three_prime_seq)/3):
      codon= upper(three_prime_seq [codon_index*3:codon_index*3+3])
      if not all([ lett in 'ACGT'  for lett in codon] ):     return 
      if codon in stops:
        ### ALRIGHT!
        if codon_index!=0:  
          aa_extension=    transl(three_prime_seq[:codon_index*3])
          self.alignment.set_sequence( 't',  self.alignment.seq_of('t')+aa_extension )
          self.alignment.set_sequence( 'q',  self.alignment.seq_of('q')+'-'*len(aa_extension) )
          if    self.strand=='+':            self.exons[-1]= [  self.exons[-1][0],                self.exons[-1][1]+codon_index*3  ]
          elif  self.strand=='-':            self.exons[-1]= [  self.exons[-1][0]-codon_index*3,  self.exons[-1][1]                ]
          if not silent:   write((self.profile.name+'.'+self.id+' ['+self.prediction_program()+'] ').ljust(22) +' completed prediction at 3\' until the next downstream stop, '+str(codon_index)+' codons added', 1)
          self.reset_derived_data() 
        return

  def complete_at_five_prime(self, max_extension=15, max_query_unaligned=30, full=False, stops=None, starts=None,  silent=False):
    """function to modify inplace the prediction to try adding some coding sequence at the 5', until the first methionine. This is done only if the extension is <= max_extension, and if there's no strange character (anything different from ACGT), and no stop codon is found before Met, and if there's at the most max_query_left aminoacid unaligned on the left side of the query. 
    If full==True, the function choose the leftmost possible ATG instead of the first one  """
    if stops is None: 
      codon_table=get_genetic_code_table()      
      stops=set([codon  for codon in codon_table if codon_table[codon]=='*'  ])
    if starts is None: starts=set(['ATG'])
    q_bounds=self.query.boundaries() #q_length=len(nogap(self.profile.seq_of(self.query_full_name())));    
    if self.protein()[0]=='M' and not full: return 
    if q_bounds[0] > max_query_unaligned: return 
    fp= self.five_prime(length=max_extension*3)
    if fp is None: return   #out of boundaries exception
    five_prime_seq= join(fp.split('\n')[1:], '')
    if len(five_prime_seq)%3:       five_prime_seq= five_prime_seq[len(five_prime_seq)%3:]   ## if not enough seq is available, we make sure it's in-frame
    atg_found=False
    for codon_index in range(0, len(five_prime_seq)/3)[::-1]:   ## reading backwards    12, 11, 10, .., 0
      codon= upper(five_prime_seq [codon_index*3:codon_index*3+3])
      if not all([ lett in 'ACGT'  for lett in codon] ):     break
      if codon in stops:                     break
      if codon in starts:
        atg_found=True ; best_codon_index=codon_index
        if not full: break 
    if atg_found:     
      aa_extension=    transl(five_prime_seq[best_codon_index*3:])
      self.alignment.set_sequence( 't',           aa_extension+self.alignment.seq_of('t') )
      self.alignment.set_sequence( 'q',  '-'*len(aa_extension)+self.alignment.seq_of('q') )
      if    self.strand=='+':            self.exons[0]=  [  self.exons[0][0]- len(aa_extension)*3,  self.exons[0][1]                      ]
      elif  self.strand=='-':            self.exons[0]=  [  self.exons[0][0],                       self.exons[0][1]+ len(aa_extension)*3 ]
      if not silent:   write((self.profile.name+'.'+self.id+' ['+self.prediction_program()+'] ').ljust(22) +' completed prediction at 5\' until the next upstream start, '+str(len(aa_extension))+' codons added', 1)
      self.reset_derived_data()       
      return 

def load_p2g(input_file, profile_ali=False, profiles_filenames_hash=False, silent=False):
  """Loads an p2g output of selenoprofiles into a p2ghit object identical to the one that generated it.  inputfile can be a filehandler, a string with the path to the file to load, or a string with the file content.  It is necessary to have the profile that generated this prediction loaded before running this function. variable profiles_hash must be defined and include the profile for this prediction or an exception will be raised. Two alternative ways to provide the profile are: directly with the variable profile_ali, or with the hash profiles_filenames_hash, which has as keys the filenames loaded and as values the profile objects.
  NB: it is also necessary to define the variables chromosome_lengths (hash of: chromosome_name -> integer length) and three_prime_length (integer) using the function set_MMlib_var :
  set_MMlib_var('chromosome_lengths', chromosome_lengths)
  set_MMlib_var('three_prime_length', three_prime_length)
  silent==True avoids warnings being printed
  OLD FUNCTION! MAYBE IT DOESN'T WORK
  """
  try:        profiles_hash=get_MMlib_var('profiles_hash')
  except:     profiles_hash={}
  if type(input_file) == file:     lines=input_file   #type: filehandler
  if type(input_file) == str:      
    if is_file(input_file):        lines=open(input_file, 'r')    #type: string with file path
    else:                          lines=input_file.split('\n')   #type: string with content of a file

  old_lines=[]
  aligned_seq_query=''; aligned_seq_target=''; next_line_is=''
  for cline in lines:
    old_lines.append(cline)
    cline=cline.rstrip() #chomping
    if cline.startswith('Output_id:'):
      output_id=del_white(cline.split('Output_id:')[1])
      family_name= output_id.split('.')[0]
      numeric_id=  output_id.split('.')[1]
      label=       output_id.split('.')[2]

    elif cline.startswith("-Species"):
      species_taxid=0
      species_name=del_white(   cline.split("-Species")[1].split('-Taxid')[0] )
      try:        species_taxid= int(del_white(cline.split('-Taxid')[1]))
      except:     
        if not silent: printerr('selenoprofiles WARNING loading p2ghit from the .p2g output '+input_file+' cannot determine taxid. Setting it to 0', 1 ) 

    elif cline.startswith("-Target"):
      abs_path_to_target=del_white(cline.split("-Target")[1])
      
    elif cline.startswith('-Chromosome'):  #-Chromosome (+)   gi|jkasdba .....
      strand=     cline.split('(')[1].split(')')[0]
      chromosome= del_white(join(cline.split(')')[1:], ')'))
    
    elif cline.startswith("-Program"):
      prediction_program=    del_white( cline.split('-Program')[1] )
      if   prediction_program=='blast':            x=blasthit()
      elif prediction_program=='exonerate':        x=exoneratehit()
      elif prediction_program=='genewise':         x=genewisehit()
      else: raise Exception, "selenoprofiles->load_p2g ERROR program not recognized: "+prediction_program

      x.species=species(species_name);        x.species.taxid=species_taxid
      x.target=abs_path_to_target
      x.strand=strand;                        x.chromosome=chromosome
      if not profile_ali and (not profiles_filenames_hash or not profiles_filenames_hash.has_key(input_file)) and (not profiles_hash or not profiles_hash.has_key(family_name)): raise Exception, "selenoprofiles->load_p2g ERROR profile not found in memory: "+family_name
      if profile_ali:                      x.profile=profile_ali
      elif profiles_hash:                  x.profile=profiles_hash[family_name]
      elif profiles_filenames_hash:        x.profile=profiles_filenames_hash[input_file]
      
      x.id=numeric_id
      x.label=label

    elif cline.startswith("-Query name"):
      query_name=   del_white(cline.split("-Query name")[1])
      x.query.chromosome=query_name 
      x.query.id=numeric_id+'_query'

    elif cline.startswith("-Query range"):
      query_positions_string=   cline.split("-Query range")[1].split()[0]
      query_start=              int(query_positions_string.split('-')[0])
      query_end=                int(query_positions_string.split('-')[1])
      x.query.add_exon(query_start, query_end)
    
    elif cline.startswith("-Profile range"):
      pass
    
    elif cline.startswith("-State"):
      filtered= del_white(cline.split("-State")[1])
      x.filtered=filtered

    elif cline.startswith("Query"):
      for piece in cline.split()[1:]:
        is_ok=True #checking if this is a piece of the alignment. It may also be a intron flag or frame shift flag, and we want to get rid of those
        if piece=='FRAME': is_ok=False
        else:
          for char in piece:
            if not char in AA_LETT+'-*X':
              is_ok = False; break
        if is_ok:  aligned_seq_query+=piece

    elif cline.startswith("Target"):
      for piece in cline.split()[1:]:
        is_ok=True #checking if this is a piece of the alignment. It may also be a intron flag or frame shift flag, and we want to get rid of those
        for char in piece:
          if not char in AA_LETT+'-*xX':
            is_ok = False; break
        if is_ok:  aligned_seq_target+=piece

    elif cline == '------- positions -------':
      ali=alignment()
      ali.add('q', aligned_seq_query)
      ali.add('t', aligned_seq_target)
      x.alignment=ali


    elif cline.startswith('Exon '):
      exon_start= int(cline.split()[2])
      exon_end=   int(cline.split()[3])
      x.add_exon(exon_start, exon_end)
          
  if type(lines)==file: lines.close()

  return x


class empty_p2g(p2ghit):
  """ Puppet class for certain cases."""

##########
#BLASTHIT

class blasthit(p2ghit):
  """ Class to keep two genes object: one is on the target (the self object), the other one is on the query (self.query). 
  Should be initiated like b=blasthit(blast_parsed_line, program='tblastn'); but it should read automatically any format.
  Also, it can accept as first argument a filename as well as a blast_parsed_line. 
  If you want to save memory and not keep the alignments in the object, you can initiate a blasthit with dont_keep_ali=1 as keyword arg (NOT TESTED!)  """
  def load_blaster(self, line, dont_keep_ali=None, id=None, keep_lengths=None):
    dont_keep_ali=dont_keep_ali or self['dont_keep_ali']
    splt=line.split()
    self.chromosome=splt[0]
    self.strand=splt[3] #not 100% sure ...
    start, stop= int(splt[1]), int(splt[2])
    if start >stop: #USELESS!
      start, stop= stop, start
      self.strand='-'
    self.add_exon(start, stop)
    #alignment    
    self.alignment=alignment()
    if not dont_keep_ali:
      self.alignment.add('q', splt[8] ) #not keeping the title here, I already have it in chromosome
      self.alignment.add('t', splt[7])
    if keep_lengths:
      self.target_length= int(splt[11])      
      self.query_length=  int(splt[12])

    self.evalue=e_v(splt[9])
    self.bits=  float(splt[10])
    if id:      self.id=id
    else:       self.id= str(uniq_id(self))
    #query
    start, stop= int(splt[5]), int(splt[6])
    self.query=gene(chromosome=splt[4], id=self.id+'_query')
    if start >stop:
      start, stop= stop, start
      self.query.strand='-'
    else:      self.query.strand='+'    
    self.query.add_exon(start, stop)
  load=load_blaster
  
  def summary(self, description='BLAST HIT', old_style=False, **keyargs):
    if self.alignment.nseq() and not old_style:      return self.pretty_summary()
    else:      return gene.summary(self, description, **keyargs)+'\n'+gene.summary(self.query, description+'_query', **keyargs)
    
  def pretty_summary(self, chars_per_line=60 ):
    '''returns a human, blast-like readable summary from a blaster_parser line referred to a single blast HSP    '''
    if not self.alignment.nseq():      raise Exception, "blasthit->pretty_summary ERROR the alignment was not saved, can't show summary for blasthit "+self.id
    q_seq, t_seq=self.alignment.seq_of(0), self.alignment.seq_of(1)
    q_start, q_end=self.query.boundaries()
    t_start, t_end=self.boundaries()
   #identifying type of blast favour (query, target--> protein? nucleotide?)
    if len(nogap(q_seq)) == abs(q_start-q_end)+1:      increase_pos_q=1
    elif 3*len(nogap(q_seq)) == abs(q_start-q_end)+1:  increase_pos_q=3
    else:        raise Exception, "blasthit->pretty_summary ERROR invalid positions/length of alignment query!\n"+self.summary(old_style=True)+'\n'+q_seq
    if len(nogap(t_seq)) == abs(t_start-t_end)+1:       increase_pos_t=1
    elif 3*len(nogap(t_seq)) == abs(t_start-t_end)+1:   increase_pos_t=3
    else:        raise Exception, "blasthit->pretty_summary ERROR invalid positions/length of alignment target!\n"+self.summary(old_style=True)+'\n'+t_seq
    strand_t,strand_q =self.strand, self.query.strand
    #computing what to put as Frame description
    frame_text='/'
    if increase_pos_q==3 and increase_pos_t==3:      frame_text=strand_q+{1:'1', 2:'2', 0:'3'}[q_start%3]+'/'+strand_t+{1:'1', 2:'2', 0:'3'}[t_start%3]
    elif increase_pos_q==3 :                         frame_text=strand_q+{1:'1', 2:'2', 0:'3'}[q_start%3]
    elif increase_pos_t==3:                          frame_text=strand_t+{1:'1', 2:'2', 0:'3'}[t_start%3]
    #begin to build output
    try: query_length=str(self.query_length)
    except: query_length='???'
    try: target_length=str(self.target_length)
    except: target_length='???'
    out=''
    out+="Query= "+self.query.chromosome+'\n         ('+query_length+' letters)\n'
    out+=">"+self.chromosome+'\n         Length = '+target_length+'\n\n'
    bits='???'
    if hasattr(self, 'bits'): bits=str(self.bits)
    out+=" Score = "+bits+" bits (???), Expect = "+str(self.evalue)+'\n'
    out+=" Identities = ???/??? ("+str(round(self.alignment.sequence_identity()*100, 1))+"%), Positives = ???/??? (??%)\n"
    out+=" Frame = "+frame_text+'\n\n'
    if strand_t=='-':      t_start, t_end=  t_end, t_start 
    if strand_q=='-':      q_start, q_end=  q_end, q_start 
    length_field = 8+max(len(str(t_start)),len(str(t_end)),len(str(q_start)),len(str(q_end)),)
    c=0
    q_line=  'Query: '+str(q_start)
    q_line+=' '*(length_field-len(q_line))
    mid_line=' '*length_field
    t_line=  'Sbjct: '+str(t_start)
    t_line+=' '*(length_field-len(t_line))
    for i in range(len(q_seq)):
      if i-chars_per_line*c >= chars_per_line:
        out+= q_line+' '+str(q_start-1*(int(strand_q+'1')))+'\n'
        out+= mid_line+'\n'
        out+= t_line+' '+str(t_start-1*(int(strand_t+'1')))+'\n'
        out+='\n'
        q_line=  'Query: '+str(q_start)
        q_line+=' '*(length_field-len(q_line))
        mid_line=' '*length_field
        t_line=  'Sbjct: '+str(t_start)
        t_line+=' '*(length_field-len(t_line))
        c+=1
      q_line+=q_seq[i]
      t_line+=t_seq[i]
      if t_seq[i]==q_seq[i]:                       mid_line+=q_seq[i]
      elif similar_aas(q_seq[i], t_seq[i]):        mid_line+='+'
      else:                                        mid_line+=' '
      if q_seq[i]!='-':        q_start+=increase_pos_q*(int(strand_q+'1'))
      if t_seq[i]!='-':        t_start+=increase_pos_t*(int(strand_t+'1'))
    out+= q_line+' '+str(q_start-1*(int(strand_q+'1')))+'\n'
    out+= mid_line+'\n'
    out+= t_line+' '+str(t_start-1*(int(strand_t+'1')))+'\n'
    return out
  
  def set_query_to_master(self, clusters_relative_positions, queries_alignment):
    """ This function changes the query attribute of the protein (being one the BLAST_QUERY of one of the cluster of the profile) to the sequence of the BLAST_QUERY_MASTER in the profile. It uses a hash of relative positions such as the one returned by the function clusters_relative_positions of profile_alignment"""
    new_query=gene()
    new_query.chromosome='BLAST_QUERY_MASTER'
    new_query_start=  clusters_relative_positions[ self.query.chromosome ][ self.query.boundaries()[0] ]
    new_query_end=    clusters_relative_positions[ self.query.chromosome ][ self.query.boundaries()[1] ]
    new_query.add_exon(new_query_start, new_query_end)
    new_query.id=    self.query.id
    new_query.strand=self.query.strand #always + actually
    self.query=new_query
    q_ali_cut = queries_alignment.columns(new_query_start-1, new_query_end-new_query_start+1) #this will work cause there are no gaps in BLAST_QUERY_MASTER seq in queries_alignment, for construnction
    ali = q_ali_cut.transfer_alignment( self.alignment, dont_shrink=True )
    self.alignment = alignment()
    self.alignment.add('q', ali.seq_of('BLAST_QUERY_MASTER'))
    self.alignment.add('t', ali.seq_of('t'))
    return 
  
  def place_selenocysteine(self):
    """ This function parses the blast output and put Us for selenocysteine: all * in query are replaced, and in the target all UGAs aligned to such positions are replaced """
    for pos in range(len(self.alignment.seq_of('q'))):
      if self.alignment.seq_of('q')[pos] in 'U*':
        self.alignment.set_sequence('q', self.alignment.seq_of('q')[:pos]+'U'+self.alignment.seq_of('q')[pos+1:] )
        if self.alignment.seq_of('t')[pos] == '*':
          codon= replace_chars( upper(self.aligned_cds()[pos*3:pos*3+3]), 'T', 'U') 
          if not codon in STOP_CODONS: raise Exception, "ERROR getting blast codon for aminoacid \"*\" . This is not a stop codon: "+codon +' '+str(self)
          if codon=='UGA':            self.alignment.set_sequence('t', self.alignment.seq_of('t')[:pos]+'U'+self.alignment.seq_of('t')[pos+1:] )

  def remove_internal_introns(self, min_length=18, silent=False):
    """ This function detects and remove the portion of blast hits which corresponds to an intron. In fact it happens often, for small introns and for flanking exons which are in the same frame respect to the full chromosome sequence , that the blast HSP extends over it. Such blast hits presents large portions of target not aligned to query. Often this target portion is full of stop codons as well. 
    The argument min_length is the minimum number of consecutive aligned nucleotides that will be called intron.
    """
    min_length_in_aa= min_length / 3
    seq_query=self.alignment.seq_of('q'); seq_target=self.alignment.seq_of('t')
    intron_start_in_ali=seq_query.find( min_length_in_aa*'-'  )
    while intron_start_in_ali != -1:
      intron_length=min_length_in_aa
      while intron_start_in_ali + intron_length <= len(seq_query) and seq_query[ intron_start_in_ali + intron_length]=='-':         intron_length+=1 #extending intron to max length possible
      number_of_gaps_in_target_before_intron =  seq_target[:intron_start_in_ali].count('-')
      gene_portion_to_remove = self.subseq(  1+(intron_start_in_ali-number_of_gaps_in_target_before_intron) *3, intron_length*3)
      a= self.subtracted_of( gene_portion_to_remove  ).copy()
      a.alignment.set_sequence('q', seq_query [:intron_start_in_ali] + seq_query [intron_start_in_ali+intron_length:])
      a.alignment.set_sequence('t', seq_target[:intron_start_in_ali] + seq_target[intron_start_in_ali+intron_length:])      
      self.__dict__ = a.__dict__
      self.reset_derived_data()
      if not silent: write( (self.profile.name+'.'+self.id+' [blast] ').ljust(23)+'internal intron removed, '+str(intron_length)+' codons removed from the alignment', 1)
      #updating these variables to check if there's another intron in the blasthit
      seq_query=self.alignment.seq_of('q'); seq_target=self.alignment.seq_of('t')
      intron_start_in_ali=seq_query.find( min_length_in_aa*'-'  )


class superblasthit(blasthit):
  """ This class is a particular blasthit, which is actually a collection of blasthits. It derives from the merge_by_colinearity method applied to blasthits.
  The purpose of this class is to manage the exceptions that arise from the fact that it does not have a single exon, it generally has more. 
  The coordinates of exons on both self (target) and query are produced by the union_with method. The alignment is produced adding "x" to fill the unaligned positions of query.
  Each superblasthit has by default a .merged list attribute storing all blasthits that are stored inside it. Delete this attribute if you don't plan to use it to save memory usage.
  """
  def __init__(self):
    blasthit.__init__(self) 
    self.merged=[] #contains the blast hits merged in this. to avoid using much memory, just delete this.
  def pretty_summary(self, chars_per_line=60, description='SUPEREXON' ):
    """ Overrides the function of the parent class to print the summary of all blast hits merged into this one."""
    o='##### '+description+' '+self.id+' ##### \n'; chars_header =len(o)
    for index, blast_h in enumerate(self.merged):
      o+=blast_h.pretty_summary(chars_per_line=chars_per_line)
      if not index==len(self.merged)-1:         o+='+'*chars_header+'\n'
    o+='#'*chars_header
    return o

class parse_blast_tab(parser):
  """ Read blast hits (without .alignment)"""
  def parse_next(self):
    if not self.last_line: self.stop()
    splt=self.last_line.split('\t')
    #target
    s, e, strand = int(splt[8]), int(splt[9]), '+' 
    if s > e: strand='-'; s,e=e,s
    g=blasthit(); g.chromosome=splt[1]; g.strand=strand  ## will return g
    g.alignment=None; #alignment(); g.alignment.add('q', 'X'); g.alignment.add('t', 'X');
    g.add_exon(  s, e  )
    #query
    s, e, strand = int(splt[6]), int(splt[7]), '+' 
    if s > e: strand='-'; s,e=e,s
    g.query.chromosome= splt[0]
    g.query.strand=strand; g.query.add_exon(s, e)
    g.evalue=e_v( splt[10] )
    g.bits=float( splt[11] )
    g.identity=float(splt[2])

    self.last_line=self.file.readline()        
    return g

class parse_blaster(parser):
  """ Parse a blaster_parser output file, which was used to scan a blast output. blasthit instances are returned on each next() call. Define dont_keep_ali=1 when calling the parser to ignore the alignments in the blast output."""
  def parse_next(self):
    g=blasthit()
    g.load(self.last_line,   dont_keep_ali=bool(self['dont_keep_ali']),   keep_lengths=bool(self['keep_lengths']) )
    self.last_line=self.file.readline()
    return g
  
class parse_blast(parser):
  """ Parse a ncbi blast output file, which is passed through blaster_parser to. blasthit instances are returned on each next() call. Define dont_keep_ali=1 when calling the parser to ignore the alignments in the blast output.
  When initialising, add keyargs full_target=1 or full_query=1 to have complete names instead of just the first word. You can use full=1 to have both
  """
  def load(self, filename=''):
    if not filename:
      filename=self.file.name
    if filename!='<fdopen>':   check_file_presence(filename, 'filename')

    add_option=''
    if self['full']:  
      add_option+=' -v FULL=1 '
      self.full_query=1; self.full_target=1
    else: 
      if self['full_query']:  add_option+=' -v FULL_QUERY=1 '
      if self['full_target']: add_option+=' -v FULL_TARGET=1 '
    if self['keep_lengths']:  add_option+=' -v QLENGTH=1 -v TLENGTH=1 '

    self.file=bash_pipe('blaster_parser.g '+add_option+' '+ filename)
    self.last_line=self.file.readline()
  def parse_next(self):
    g=blasthit()
    g.load(self.last_line,   dont_keep_ali=bool(self['dont_keep_ali']),   keep_lengths=bool(self['keep_lengths']))
    if self['full_query']:    
      g.query.chromosome=  replace_chars(g.query.chromosome, '%', ' ')

    if self['full_target']:   
      g.chromosome=  replace_chars(g.chromosome, '%', ' ')
    self.last_line=self.file.readline()          
    return g

##########
#EXONERATEHIT

class exoneratehit(p2ghit):
  """ This class handles the exonerate output predictions """

class superexoneratehit(exoneratehit):
  """ This class handles the superexonerate hits, meaning sets of exonerate prediction merged by colinearity into a single one. 
  These class is loaded typically from a exonerate output file which was generated by cyclic_exonerate, so it has a seed gene object which was used to initiate the procedure.
  Only the best scoring exoneratehit (in case, extended by the merging procedure) which overlaps the original seed is loaded. The original seed information is read from a specially formatted comment line, the first one of the file. 
  In case you don't have the comment line in the file, you can still specify the seed with the seed option in the load method.
  """

  def __init__(self, **keyargs):
    exoneratehit.__init__(self) #not tested!
    self.merged=[] #contains the exonerate hits merged in this. to avoid using much memory, just delete this.

  def load(self, filename, seed='', query_full_sequence='', merge_multiple=True):
    """ see doc string of the class"""
    self.__init__() #flushing data
    check_file_presence(filename, 'exonerate file')
    fileh= open(filename, 'r')
    first_line= fileh.readline();              second_line= fileh.readline()
    fileh.close()
    if not seed:      
      if not (first_line.startswith('#') and 'seed:' in first_line): raise Exception, "superexoneratehit->load ERROR no seed was specified and no comment line was detected in the file "
      seed=gene()
      try:    seed.load_from_header(del_white( join(first_line.split('seed:')[1:], 'seed:'))) 
      except: printerr("superexoneratehit->load ERROR trying to load the seed from the comment line", 1); raise
    all_exonerate_hits_in_file=   parse_exonerate(filename).all()
    if not all_exonerate_hits_in_file: return 'empty' #message for: file is empty.  the self object is empty too
    if not query_full_sequence:
      try:    query_full_sequence=del_white(second_line[:-1].split('query_full_seq:')[1].split(';')[0])
      except: printerr("superexoneratehit->load ERROR trying to load the query sequence from the comment line", 1); raise
    fake_ali=alignment()   #building fake alignment just to pass the full query sequence to merge_p2g_hits_by_colinearity function
    if not all(all_exonerate_hits_in_file): return 'output with gff bug'
    fake_ali.add(all_exonerate_hits_in_file[0].query.chromosome, replace_chars(query_full_sequence, '*', 'U')   ) 
    if merge_multiple:        
      merge_p2g_hits_by_colinearity(all_exonerate_hits_in_file, inplace=True, post_function=merge_p2g_hits_by_colinearity_post_function_exonerate_hits,  sequence_collection=fake_ali)
    for exonerate_or_superexoneratehit in sorted(all_exonerate_hits_in_file, key=lambda e:e.score, reverse=True):
      #navigating the exonerate hits in the merged list, taking first the best scoring ones. We will stop when we find the first one overlapping the seed.
      if exonerate_or_superexoneratehit.overlaps_with(seed, phase=False):
        if  exonerate_or_superexoneratehit.__class__.__name__=='superexoneratehit':
          self.__dict__=deepcopy(exonerate_or_superexoneratehit.__dict__)
        elif exonerate_or_superexoneratehit.__class__.__name__=='exoneratehit':
          for i in exonerate_or_superexoneratehit.__dict__:      self.__dict__[i]=exonerate_or_superexoneratehit.__dict__[i]
          self.query=exonerate_or_superexoneratehit.query
        self.join_adjacent_exons() #adjust the prediction for rare cases
        return 'ok'      
    #if no exonerate or superexonerate was found overlapping the seed: do nothing.. the exonerate object will be empty.    
    return 'none_overlapping'  #message for: no exonerate hit was found overlapping.  the self object is empty
  
  def cds(self, **keyargs): #**keyargs are not used. but we accept them to allow compatibility with the cds function in the upper class p2ghit 
    """ This returns the coding sequence of the entire prediction. Frameshifts nucleotide are excluded so that the translation of the cds is always equal to the predicted protein sequence"""
    #if self.merged:
    #  out=''
    #  for e in self.merged: 
    #    out+=p2ghit.cds(self)
    #  return out        
    #else: 
    return p2ghit.cds(self)

class parse_exonerate(parser):
  """ Parse an exonerate file and returns a exoneratehit object for each prediction inside. 
  The target names obtained by fastasubseq are recognized and set to absolute coordinates. Also the target names obtained through the method from gene class "fasta_sequence", with title set to "fasta_title", are recognized (they are used in selenoprofiles).
  """
  def load(self, filename=''):
    if not filename:  filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file = open(filename, 'r')
    self.last_line=self.file.readline()
    
  def parse_next(self):
    line=self.last_line;     cfile=self.file #necessary to recycle code: this below is the old exonerate parser
    # reading inputfile
    cont_b=0
    while line and line !='-- completed exonerate analysis\n':
      query_three_lett_seq='';   three_lett_seq='';   target_dna_seq='';   ali_target_seq='';   ali_query_seq=''; frameshifts_data=[]; # contains elements like [position_in_the_cds_1_based, length_of_insertion]
      while line and line !='-- completed exonerate analysis\n' and line.split(':')[0]!="  Target range":
        if line!='\n' and line.split()[0]=='Target:':   full_target_name=del_white(line[:-1].split('Target:')[1])
        line=cfile.readline()
      if line and line !='-- completed exonerate analysis\n':
        ### HERE WE ARE IN THE FIRST LINE (BEGINNING WITH Target range) OF A EXONERATE PREDICTION. THIS BLOCK OF CODE IS TO FILL THE QUERY_THREE_LETT_SEQ , THREE_LETT_SEQ AND TARGET_DNA_SEQ
        skip_intron=False
        passed_the_N_of_intron=False
        line=cfile.readline()
        line=cfile.readline()
        while line and line.split(':')[0]!='vulgar':
          align_block=[]        # reading 5 lines-alignment blocks
          for i in range(5):
            align_block.append(line)
            line=cfile.readline()
          ccc=1            #index of the current char in the line
          query_block=align_block[0].split(':')[1]+'}'
          piece_length= len(align_block[0].split(':')[0])+1 #lenght of string til ":" in query line
          target_block=align_block[2][piece_length:]
          target_dna_block=align_block[3][piece_length:]
          if align_block[3].find(':') != align_block[0].find(':'):             query_block= query_block.strip()  #fixing a rare bug in which query is not aligned with targert because of wrong positinos are printed (999 and 1000)
          
          while ccc<len(query_block):
            ### procedure: each splitted piece of text is tried to be added to the growing query_three_lett_seq. At each add, the seq is checked in its last codon (actualli the 3 lett code aminoacid). If it is not a valid codon, the last seq added is taken away. A special case is made for the ">" char.
            just_added=False
            char=query_block[ccc]
            if skip_intron:
              if passed_the_N_of_intron and symbols_after_intron_found==4 :              skip_intron=False
              elif char =='n':                                                          passed_the_N_of_intron=True
              elif char =='>' and passed_the_N_of_intron:                                symbols_after_intron_found +=1
            else:
              if char in ' >{#}+' or is_number(char) or (char =='-' and   ( (ccc>1  and  query_block[ccc-1]!='<') or  (ccc==1 and  query_three_lett_seq[-1]!='<') )) or  (char=='-' and target_block[ccc]=='#'):
                for i in range(cont_b):
                  query_three_lett_seq+=query_block[ccc-(cont_b-i)]
                  three_lett_seq+=target_block[ccc-(cont_b-i)]
                  target_dna_seq+=target_dna_block[ccc-(cont_b-i)]
                just_added=True
                if char=='-' and target_block[ccc]=='#': #annotating frameshifts
                  if frameshifts_data and frameshifts_data[-1][0]==len(nogap(target_dna_seq))+1:
                    previous_frameshift_length=frameshifts_data[-1][1]
                    frameshifts_data.pop(-1)
                    frameshifts_data.append([len(nogap(target_dna_seq))+1, previous_frameshift_length+1])
                  else: frameshifts_data.append([len(nogap(target_dna_seq))+1, 1])
              else:              cont_b+=1
              if char=='>':
                if query_three_lett_seq[-1]=='-' and query_three_lett_seq[-2]=='<':
                  query_three_lett_seq+=query_block[ccc]
                  three_lett_seq+=target_block[ccc]
                  target_dna_seq+=target_dna_block[ccc]
                else:
                  skip_intron=True
                  passed_the_N_of_intron=False
                  symbols_after_intron_found=0
              alright=False
              while not alright and len(query_three_lett_seq)>=3:
                a = len(query_three_lett_seq) % 3 
                if a==0:                  last_codon = query_three_lett_seq[len(query_three_lett_seq)-3:][:3]      
                else:                    last_codon = query_three_lett_seq[len(query_three_lett_seq)-3-a:-a][:3]  
                if not three_letter_codon_diz.has_key(last_codon) :
                  if cont_b==0:                  cont_b=1
                  query_three_lett_seq=query_three_lett_seq[:-cont_b] 
                  three_lett_seq = three_lett_seq[:-cont_b]
                  target_dna_seq = target_dna_seq[:-cont_b]
                else:                  alright=True
                if just_added:         cont_b=0
            ccc+=1

        ### END OF CODE BLOCK
        for codon_index in range(   len(three_lett_seq) / 3 ) :
          ali_target_seq+=  three_letter_codon_diz[  three_lett_seq[3*codon_index:3*codon_index+3 ]    ]    #translating 3 letters code to one letter
        for codon_index in range(   len(query_three_lett_seq) / 3 ) :
          ali_query_seq+=  three_letter_codon_diz[  query_three_lett_seq[3*codon_index:3*codon_index+3 ]    ]    #translating 3 letters code to one letter
        #now in vulgar line; we parse positions from here.
        vulgar_line=line
        qname=vulgar_line.split()[1] ;      qstart=int(vulgar_line.split()[2])    + 1   ;     qend=int(vulgar_line.split()[3])   
        tname=vulgar_line.split()[5];       real_tstart=int(vulgar_line.split()[6])    +1 ;    real_tend=int(vulgar_line.split()[7])
        if real_tstart> real_tend:    #negative frame
          real_tstart = real_tstart-1
          real_tend = real_tend+1
        raw_score=vulgar_line.split()[9]
        ali=join(vulgar_line.split()[9:], ' ')
        frameshifts=len(ali.split('F') )-1
        line=cfile.readline()
        # processing considering subseqing
        #now in GFF comment startline 
        if line=="# --- START OF GFF DUMP ---\n":
          while line[0]=='#':
            line=cfile.readline()
          #now on actual gff line
          gff_lines=''
          while line[0]!='#':
            gff_lines+=line
            line=cfile.readline()
        else:        raise Exception, 'ERROR no gff ouput in the exonerated file; please run exonerate with --showtargetgff option'
        while line and line !='-- completed exonerate analysis\n' and line !='C4 Alignment:\n':          line=cfile.readline()  
        #replacing Sec-TGAs with U.. replacing also all * in query with U
        
        for pos in range(len(ali_query_seq)):
          if ali_query_seq[pos]=='*':
            ali_query_seq=ali_query_seq[:pos]+'U'+ali_query_seq[pos+1:]
            if  ali_target_seq[pos]=='*' and upper(target_dna_seq[pos*3:pos*3+3])=='TGA':  ali_target_seq=ali_target_seq[:pos]+'U'+ali_target_seq[pos+1:]            
        e=exoneratehit()
        ex_mode= gff_lines.split('\t')[1] #reading from first line
        if 'protein2dna' in ex_mode:       gff_tag='similarity'
        elif 'protein2genome' in ex_mode:  gff_tag='cds'
        else:   raise Exception, 'ERROR I do not know how to read the gff! The exonerate mode was not recognized: '+ex_mode
        e.load_gff(gff_lines, tag=gff_tag)
        e.alignment=alignment()
        e.alignment.add('q', ali_query_seq );       e.alignment.add('t', ali_target_seq )
        e.cds_sequence=nogap(target_dna_seq)
        e.query=gene(chromosome=qname, strand='+')
        e.query.add_exon(qstart, qend)
        e.id=str(uniq_id(e))
        e.query.id=str(uniq_id(e))+'_query'
        e.score=int(raw_score)
        #e.frameshifts_data=frameshifts_data
        if ':subseq(' in tname:
          strand_subseq='+'
          if '[revcomp]' in full_target_name:              strand_subseq='-'
          start_subseq=int(tname.split(':subseq(')[1].split(',')[0]  ) # 0 based
          length_subseq=int(tname.split(':subseq(')[1].split(')')[0].split(',')[1]  )
          e.chromosome=tname.split(':subseq(')[0]
          subseq_gene=gene(chromosome=tname, strand=strand_subseq)
          subseq_gene.add_exon(start_subseq+1, start_subseq+length_subseq)
          e.restore_absolute_coordinates(parent_gene=subseq_gene)
        elif '[positions:' in tname:
          subseq_gene=gene()
          subseq_gene.load_from_fasta_title(tname)
          e.chromosome=tname.split('[positions:')[0]
          e.restore_absolute_coordinates(parent_gene=subseq_gene)
        if line=='-- completed exonerate analysis\n':  line=''
        self.last_line=line
        if frameshifts_data:  #introducing the frameshifts in the gff (self.exons) 
          for pos_start, insertion_length in frameshifts_data:
            portion_to_remove=e.subseq(pos_start, insertion_length)
            e=e.subtracted_of(portion_to_remove)
        
        if e.length() != len(e.cds_sequence):         
            #try to correct a bug by exonerate. in some occasions the gff and the visual alignment are not showing the same, with one codon missing from the gff. here i'm correcting a quite specific bug...
          try:
            difference_in_length=len(e.cds_sequence) - e.length() 
            for i in range(0, len(e.cds_sequence), 3):
              codon_in_cds_sequence=upper(e.cds_sequence[i:i+3])
              codon_as_in_gff= upper(e.subseq(i+1, 3).fasta_sequence()[1])
              if codon_in_cds_sequence != codon_as_in_gff:
                #trying different possibilities, the maximum seq missing at once is e.cds_sequence- e.length()
                for difference_in_length in range(  len(e.cds_sequence)- e.length(), 0, -3   ):
                  codon_just_upstream_obj = e.subseq(i+1, 3).upstream(0, difference_in_length) 
                  if not codon_just_upstream_obj.overlaps_with(e) and upper( codon_just_upstream_obj.fasta_sequence()[1] ) == upper(e.cds_sequence[i:i+difference_in_length]):
                    for exon_index in range(len(e.exons)): #getting which exon we're extending
                      if e.get_exon(exon_index).is_downstream_of(codon_just_upstream_obj):                      break
                    if   e.strand=='+':            e.exons[exon_index][0]= e.exons[exon_index][0]-difference_in_length
                    elif  e.strand=='-':            e.exons[exon_index][1]= e.exons[exon_index][1]+difference_in_length
                    break
          except: pass

        if e.length()!= len(e.cds_sequence):  
          printerr('WARNING error in the exonerate file: '+(self.file.name)+' ; the gff does not have a perfect correspondence with the alignment displayed. To avoid crash, returning None like if it were an empty exonerate prediction', 1)
          return None
        return e

    self.last_line='' #no output found. setting last line to '' to induce the .stop method in the .next method.
      
##########
#GENEWISE
class genewisehit(p2ghit):
  """ This class handles the genewise output prediction, as output by the command: genewise  -pretty -sum -gff """
  def __init__(self, filename=''):
    p2ghit.__init__(self)
    if filename:  self.load(filename)
  def load(self, filename):
    check_file_presence(filename, 'genewise file')
    g=parse_genewise(filename).all()[0]
    self.__dict__=g.__dict__.copy()
  def check_alignment(self):
    """ This functions checks whether the alignment shown by genewise have some wrong characteristics and it correct them.
    For example sometimes genewise outputs alignments ending like:
gi|160896338|re  135 TKFLVGRDGQVIRRYAPQDAPAKLSTDIEAALAL                
                     TK+LV  DG   +RY+    P  +  DI AAL                  
                     TKWLV-VDGTPTKRYSYDVKPEAIEADIAAAL--                
scaffold_6[posi10424 aatcg gggacaactttggacggagggagggt                  
                     cagtt tagcccagacaatacactacatccct                  
                     gggcc cccggggccccccgcgccgcccccgg                  

    With the last two positions which are not informative and actually are not even reported in the gff they provide
    """
    if self.alignment.titles(): #may also be empty, when no alignment overlapping the seed is found. In this case, we have nothing to do here
      title1, title2= self.alignment.titles()
      while self.alignment.seq_of(title2)[-1]=='-':
        self.alignment.set_sequence(title1,    self.alignment.seq_of(title1)[:-1]  )
        self.alignment.set_sequence(title2,    self.alignment.seq_of(title2)[:-1]  )
        self.query[-1][1]-=1
    
class parse_genewise(parser):
  """ Parse a genewise output produced with -pretty -gff -sum and return a genewise object with all useful information.
  The target names obtained through the method from gene class "fasta_sequence", with title set to "fasta_title", are recognized (they are used in selenoprofiles).
  NB: if a seed description is found in the first commented lines, this is loaded and its overlaps with the genewise hit is tested: if there is not overlap, an empty genewise object is returned, with an .error_message attribute filled with  "Not overlapping with seed" """

  def __init__(self, filename='',  **keyargs):
    self.skip_comments=False #necessary not to skip the first lines
    for key in keyargs:      self[key]=keyargs[key]
    if filename:             self.load(filename)
  def load(self, filename=''):
    if not filename:  filename=self.file.name
    check_file_presence(filename, 'genewise output filename')
    self.file = open(filename, 'r')
    self.last_line=self.file.readline()
    
  def parse_next(self):
    def remove_ending_white_spaces(stringg):  #and \n also. this small function is used only here.
      stringg=list(stringg)
      if stringg:
        while stringg and (stringg[-1] in ['\n', ' '] ):        stringg.pop(-1)
      return join(stringg, '')  
    line=self.last_line;     cfile=self.file #necessary to recycle code: this below is the old genewise parser
    ## reading comment lines with information about the seed and original range 
    original_range=False
    if line.startswith("#") and 'seeded' in line: #first comment line. Example: '#blast_seeded; Target_range:19 chromosome:15 strand:- positions:91774890-91775531
      seed_program=line.split('#')[1].split('_seeded;')[0]
      original_range_header=line[:-1].split('Target_range:')[1]
      original_range=gene()
      try:    original_range.load_from_header(original_range_header)
      except: printerr("parse_genewise ERROR trying to load the original genomic range from the comment line", 1); raise
      line=cfile.readline()
    full_query_name=False;  full_query_sequence=False;  profile_name=False
    if line.startswith("#") and 'profile_alignment:' in line:
      full_query_name= line.split('; query_name: ')[1].split(';')[0]
      full_query_sequence=line.split('; query_full_seq: ')[1].split()[0]
      profile_name=del_white(line.split('profile_alignment: ')[1].split(';')[0])
      if ' ; target_for_this_genewise_run: ' in line:   ### this is the norm; but we keep the 'if' for compatibility
        range_fasta_title=line.rstrip().split('; target_for_this_genewise_run: ')[1].split(';')[0]
      line=cfile.readline()
    while line and line!="See WWW help for more info\n":       line=cfile.readline()
    if not line:       raise Exception, "parse_genewise: some error occured parsing " +self.file.name+' ; empty file.'
    line=cfile.readline()
    currently_intron=0;    aa_seq='';     aa_query='';     dna_target_seq=''
    # parsing pieces of aligmnment of 8 lines. saving protein alignment and dna sequence. The parsing procedure is not easy to understand but well tested.
    line=cfile.readline()
    while line !='//\n' and line:
      alignment_line=[]
      while len(alignment_line)<6:
        if len(alignment_line)==3 and line[0]==' ': # TO AVOID BUG OF RANDOM INSERTION OF NEWLINE INTO THE ALIGNMENT BY GENEWISE
          done_corrected_bug=0
          l=0
          while l<3 and not done_corrected_bug:
            if (  not any( alignment_line[l][:-1].split() )  ):
              alignment_line.pop(l)
              done_corrected_bug=1
            l+=1
        alignment_line.append(line)
        line=cfile.readline()
      line=cfile.readline()
      line=cfile.readline()
      end_index= 21+len(remove_ending_white_spaces(alignment_line[0][21:]) )
      aa_seq_portion= alignment_line[2][21: ]
      aa_seq_query_portion=alignment_line[0][21: end_index]
      dna_target_lines=[alignment_line[3][21:], alignment_line[4][21:], alignment_line[5][21:]]
      index=0
      for s in aa_seq_portion.split():
        if not ':' in s:
          aa_seq+=s
          aa_query += aa_seq_query_portion.split()[index]
          index+=1
        else:
          aa_seq+=s.split(':')[1].split('[')[0]
          aa_query+=s.split(':')[0]
      for pos in range(len(dna_target_lines[0])):
        current_column=  dna_target_lines[0][pos]+dna_target_lines[1][pos]+dna_target_lines[2][pos]
        if "<" in current_column:          currently_intron=1
        if not currently_intron and del_white(join(current_column.split('\n'), '')):          dna_target_seq+=del_white(current_column.upper())
        if ">" in current_column:          currently_intron=0
    line=cfile.readline();     line=cfile.readline()
    ### now on the sum line like: Bits   Query         start end Target      start end   idels introns   ##NOT THE HEADINGS, THE ACTUAL LINE
    found=False  
    if line:
      splitted=line.split()
      score=splitted[0]
      query=splitted[1]
      query_range_start, query_range_end = int(splitted[2]), int(splitted[3])
      target=splitted[4]
      target_range_start, target_range_end = int(splitted[5]), int(splitted[6])
      ## genewise weird bug: the target end is WRONG in this line. rescuing this
      if len(no_gap(aa_query) ) != (query_range_end-query_range_start+1):         query_range_end = len(no_gap(aa_query))+ query_range_start -1
      found=True  
    line=cfile.readline();     line=cfile.readline()
    ### now on the very first line of the gff output
    if line:
      out_gff=''
      while line!='//\n':
        out_gff+=line
        line=cfile.readline()
    ## preparing frameshifts data
    aa_seq=join(aa_seq.split('!'), '-')
    frameshifts_data=[] # contains elements like [position_in_the_cds_1_based, length_of_insertion]
    for p in range(len(dna_target_seq)):
      pos=p-len(frameshifts_data)
      if dna_target_seq[pos] in  '12345':
        frameshifts_data.append( [pos+1, int(dna_target_seq[pos])]  )
        dna_target_seq=dna_target_seq[:pos]+dna_target_seq[pos+1:]
    #correcting dna_target_seq to be aligned.
    all_pos_of_gaps=[]
    f=find(aa_seq, '-')
    while f!=-1:
      all_pos_of_gaps.append(f)
      f=find(aa_seq, '-', f+1)
    for pos in all_pos_of_gaps:      dna_target_seq=dna_target_seq[:pos*3]+'---'+dna_target_seq[pos*3:]
    #replacing U with * in unaligned TGAs and replacing also other stop codons to *
    for pos in range(len(aa_seq)):
      if    aa_seq[pos]=='U' and not aa_query[pos]=='U':                                                                 aa_seq=aa_seq[:pos]+'*'+aa_seq[pos+1:]
      elif  aa_seq[pos]=='X' and replace_chars(upper(dna_target_seq[pos*3:pos*3+3]), 'T', 'U') in STOP_CODONS:           aa_seq=aa_seq[:pos]+'*'+aa_seq[pos+1:]
    ## closing file
    while line:         line=cfile.readline()
    self.last_line=line
    ### g is the genewise object that is returned
    g=genewisehit()
##### 
#  ## expanding this:  g.load_gff(out_gff[:-1])    to check for yet another genewise bug
#               ---IPK 
#                  IPK 
#               EY!IPK 
# Intron 6   AAGgt2aca      --> small exon before frameshift will have weird coordinates
#7080 : 7509]-0>aa tca 
#               at taa 
    try:       g.load_gff(out_gff[:-1])
    except:           
      d=genewisehit()
      d.error_message='gff not valid. WARNING this is a bug of genewise'
      return d
    # restoring absolute coordinates
    if target=='gw_target':   ## this is the norm from version 3.2e
      target=range_fasta_title

    if '[positions:' in target:
      subseq_gene=gene()
      subseq_gene.load_from_fasta_title(target)
      g.restore_absolute_coordinates(parent_gene=subseq_gene)
      target=target.split('[positions:')[0]
    if original_range: 
      if not original_range.overlaps_with(g, phase=False): 
        g=genewisehit(); g.error_message='Not overlapping with seed'
        return g
    g.id=str(uniq_id(g)) #setting id of genewisehit to a uniq numeric id (used as string)  taken from the python environment
    g.chromosome=target
    g.query=gene(chromosome=query, strand='+')
    g.query.add_exon(query_range_start, query_range_end)
    g.query.id=g.id+'_query'
    #setting alignment property of g
    g.alignment=alignment()
    g.alignment.add('q', aa_query)
    g.alignment.add('t', aa_seq)
    g.cds_sequence=nogap(dna_target_seq)
    #g.frameshifts_data= frameshifts_data   
    if frameshifts_data:
      # correcting possible presence of - to - aligned positions
      g.alignment.remove_useless_gaps()
    # genewise output routine has definitely some problems. the gff doens't always correspond to the aligment output. Here I correct this
    if g.length() > len(g.protein()*3):         g=g.subseq(  1, len(g.protein()*3)  )
    if frameshifts_data and not len(frameshifts_data) != len([1   for st, end in g.introns().exons if end-st+1<7]):
      # correcting genewise bug! extending introns for certain frameshifts. anyway it must not interfer with another bug correction (below) so I put the condition above to tell which one is screwing us now.
      current_frameshift_index=-1
      introns= g.introns()
      for intron_index in range(len(introns.exons)):
        intron_length= introns.exons[intron_index][1] - introns.exons[intron_index][0] + 1
        if intron_length<7:
          current_frameshift_index+=1
          if intron_length != frameshifts_data[current_frameshift_index][1]:
            nt_need_to_add= intron_length-frameshifts_data[current_frameshift_index][1]
            if g.strand=='+':    g.exons[intron_index][1]+=nt_need_to_add
            elif g.strand=='-':  g.exons[intron_index][0]-=nt_need_to_add

    if g.length()!= 3*len(g.protein()): 
      try:
        current_exon_index=0; target_non_gaps=0; total_exon_length=g.exons[0][1]-g.exons[0][0]+1
        for ali_pos in range(g.alignment.length()):
          aa_target=g.alignment.seq_of('t')[ali_pos]
          if aa_target!='-': 
          
            codon_obj=g.subseq(target_non_gaps*3+1, 3)
            codon=upper(codon_obj.fasta_sequence()[1])
            if upper(g.cds_sequence[target_non_gaps*3:target_non_gaps*3+3])!=codon:
              if current_exon_index!=0:  raise
              if upper( g.get_exon(current_exon_index).downstream(0, 3).fasta_sequence()[1] ) == upper(g.cds_sequence[target_non_gaps*3:target_non_gaps*3+3]):
                if g.strand=='+':   g.exons[current_exon_index][1]+=3
                elif g.strand=='-': g.exons[current_exon_index][0]-=3
                total_exon_length+=3
                if g.length()== 3*len(g.protein()): break
            target_non_gaps+=1
          if target_non_gaps*3 > total_exon_length:
            current_exon_index+=1
            total_exon_length+=g.exons[current_exon_index][1]-g.exons[current_exon_index][0]+1
        assert g.length()== 3*len(g.protein())
      except:
        #printerr('WARNING there are some bugs in the genewise output: ' +self.file.name+' ; can\'t parse it. Returning an empty object. ', 1)
        d=genewisehit()
        d.error_message='alignment and gff don\'t correspond. WARNING this is a bug of genewise'
        return d
    g.remove_terminal_uga()    
    return g
    
##########
#SECIS or generic features
class p2g_feature(gene):
  """ Generic class for a feature of a prediction: it can be a label for a signal sequence, or a secondary structure, anything. It is a subclass of gene so you can check overlaps and stuff, if you use the .exons attributes. Otherwise it can be used even in simpler way, just to label predictions. Important methods that must be defined:
  - dump_text()            --> returns a text containing all information relevant to the object, which is used to store the obj in the sqlite database
  - load_dumped_text(txt, [p2g])  --> load that text back into the self object. To be used on an empty instance of this class. Note that a .parent attribute is already added when loading fromt he database, pointing to the relevant p2ghit
  - output()               --> used to output this object in the p2g output (only if included_in_output, see below)
  - gff( **keyargs)        --> used to output the feature in the gff output (only if included_in_gff, see below); if the feature attributes inherited from the gene class are used (.exons, .chromosome, .strand etc), this is not even necessary, as a native gff method is defined. The **keyargs are needed for compatibility
  Important class (not instances!) attributes:
  - included_in_output boolean, if true, the output() method is used to include this in p2g output. if this is false, there's no need to define the method output
  - included_in_gff    boolean, if true, the feature will be included in the gff file, with a line for each feature instance
  - program            string,     shown in the gff tag for this features in gff output
  - gff_add            list of strings.  Each element is an attribute name (must be present in all feature objects) that you want to be added in the last field of the gff for this kind of feature. You can use also methods, in this case, use M:method_name
  """

class protein_motif(p2g_feature):
  """ protein motif is an example of a p2g_feature, to annotate the positions of a certain motif defined as a perl-style regexp.  The motif is defined in the line following this, as a class attribute. In the example, the redox box (CXXC) is the motif. 
      Attributes:
      - start      start of the protein motif in the protein sequence (1-based, included)
      - end        end of protein motif in the protein sequence (1-based, included)
      - sequence   motif sequence 
  """
  motif=re.compile( 'C..C' )
  included_in_output=True
  included_in_gff=   True  

  def dump_text(self):    
    """ Returns a string with all the information for this feature. This string is stored in the sqlite database. """
    return str(self.start)+':'+str(self.end)+':'+self.sequence

  def load_dumped_text(self, txt):
    """ Reverse the dump_text method: gets a string as input, and loads the self object with the information found in that string. """
    start, end, sequence= txt.split(':')
    self.start= int(start);    self.end=int(end);    self.sequence=sequence

  def output(self):           
    """ Returns a string. This will be added to the p2g output of the prediction to which this feature is linked -- if class attribute included_in_output is True"""
    return 'Motif: '+self.sequence+' Start: '+str(self.start)+' End: '+str(self.end)

  def gff(self, **keyargs): 
    """This must return a gff-like tab-separated string. In this case, we are exploiting and overriding the gff method of the gene class, which is a parent class for p2g_feature"""
    ## getting a gene object with the genomic coordinates of the protein motif. we use the gene method subseq, which returns a subsequence of the parent gene. Indexes are adjusted for protein-nucleotide conversion
    motif_gene_object= self.parent.subseq(    start_subseq=(self.start-1)*3 +1,      length_subseq=(self.end-self.start+1)*3,      minimal=True    )
    #now motif_gene_object has a .exons attributes with the genomic coordinates of the protein motif. now we can use the native gff method of the obtained gene object
    return gene.gff(motif_gene_object,    **keyargs)

  def reset(self):
    """ This method is called when the linked prediction is modified, to allow to recompute some or all attributes of the feature. In this case, we are removing all features of this class, and annotating them again with the same method used to add them in first place: annotate_protein_motif"""
    ##removing instances of this class
    for index_to_remove in [ index     for index, f in enumerate(self.parent.features) if f.__class__ == protein_motif ][::-1]:        self.parent.features.pop(index_to_remove)      
    #reannotating
    annotate_protein_motif( self.parent, silent=True )

def annotate_protein_motif(p, silent=False):
  """p is a p2ghit. This is an example of method to annotate the p2g_feature protein_motif. To use, add this to the main configuration file:  
  ACTION.post_filtering.annotate_motif =    if x.filtered == 'kept':  annotate_protein_motif(x)   
  """
  s= protein_motif.motif.search(  p.protein() )   ##using search method of re.RegexObject  --  protein_motif.motif is such an object
  while s:
    protein_motif_instance=         protein_motif()
    protein_motif_instance.start=   s.start()+1   #making 1 based
    protein_motif_instance.end=     s.end()       #making 1 based and included, so it'd be +1-1 
    protein_motif_instance.sequence=      p.protein()   [ protein_motif_instance.start-1 : protein_motif_instance.end ]
    p.features.append(protein_motif_instance)     ## adding feature to p2g object
    if not silent:    printerr('annotate_protein_motif found a motif: '+protein_motif_instance.output()+' in prediction: '+p.output_id(), 1)
    s=protein_motif.motif.search(  p.protein(), pos= s.start()+1 ) ## searching again, starting from just right of the previously found position


class misc_sequence_data(p2g_feature):
  """ Hidden feature to store the sequence information about the sequence at the 3' of the prediction, until the first stop codon. It is loaded by method compute_misc_sequence_data """
  included_in_output=False
  included_in_gff   =False 
  # self.full_seq                     string
  # self.extension_downstream         int
  # self.extension_upstream           int
  # self.parent                       p2g

  def dump_text(self):    return str(self.extension_upstream)+':'+str(self.extension_downstream)+':'+self.full_seq

  def load_dumped_text(self, txt):
    extension_upstream, extension_downstream, full_seq=   txt.split(':')
    self.extension_upstream=   int(extension_upstream)
    self.extension_downstream= int(extension_downstream)
    self.full_seq=             full_seq
    ## reading parent attribute to set coordinates
    self.chromosome= self.parent.chromosome;     self.strand= self.parent.strand;      self.exons= list(self.parent.exons)
    if self.strand=='+':   self.extend( left=self.extension_upstream,  right=self.extension_downstream,    inplace=True)
    elif self.strand=='-': self.extend( right=self.extension_upstream, left=self.extension_downstream,     inplace=True)
    
  def reset(self):    
    ## resetting without recomputing, since it will be recomputed only if sequence data is requested. in this way, consecutive calls to reset will not causing repeating the computation
    self.parent.remove_misc_sequence_data()

class bsecis(p2g_feature):
  """ Class for SECISearch3/Seblastian secis predictions. Instanciated just with method load_secis (see below) """
  included_in_output=True
  included_in_gff   =True
  program='bSeblastian'
  gff_add=[ 'model_name', 'infernal_score'] 

  def seq(self):    return nogap(   self.ali.seq_of('t'))
  def ss(self):     return nogap(self.ss_ali.seq_of('t'))
  def pretty_alignment(self):
    o= '         '+self.ss_ali.seq_of('q')+'\n'
    o+=' model   '+self.ali.seq_of('q')+'\n'
    o+='         '+self.ss_ali.seq_of('t')+'\n'
    o+=' target  '+self.ali.seq_of('t')
    #o+='         '+self.annotation_line
    return o
  def summary(self):
    o='-SECIS: '+self.id+'\n'
    o+=' Positions:'+self.positions_summary()+' '
    #d=self.distance_from_cds() 
    #if d: o+='    CDS-Distance:'+str(d)+' '
    #d=self.distance_from_sec_uga() 
    #if d: o+='    Sec-Distance:'+str(d)+' '
    #o+='\n'+self.found_by_line  # has terminal \n
    o+=' Model = '+self.model_name+'  Infernal score = '+ str(self.infernal_score)+'\n'
    o+=' Apical loop length = '+str(self.apical_loop_length)+'  Upper stem length = '+str(self.upper_stem_length)+' UGA-loop spacing = '+str(self.uga_loop_spacing)+'\n'
    o+=' Free Energy: '+str(self.energy)+'\n\n'
    o+=self.pretty_alignment()
    return o
  __str__=  summary
  __repr__= summary
  output=summary

  def dump_text(self):
    """see above, function to save the object in a db """
    o= self.header()
    o+='#'+self.model_name+' '+str(self.infernal_score)
    o+='#'+self.ali.seq_of('q')+' '+self.ali.seq_of('t')+' '+self.ss_ali.seq_of('q')+' '+self.ss_ali.seq_of('t')
    o+='#'+str(self.energy)+' '+ str(self.apical_loop_length)+' '+str(self.upper_stem_length)+' '+str(self.uga_loop_spacing)
    return o

  def load_dumped_text(self, txt):
    """see above, functino to load object from a db"""
    splt=txt.split('#')
    self.load_from_header(splt[0])
    self.model_name= splt[1].split()[0]; self.infernal_score=float(splt[1].split()[1])
    splt_ali= splt[2].split()
    self.ali=alignment();     self.ali.add('q', splt_ali[0]    );   self.ali.add('t', splt_ali[1]); 
    self.ss_ali=alignment();  self.ss_ali.add('q', splt_ali[2] );   self.ss_ali.add('t', splt_ali[3]); 
    splt_rest=splt[3].split()
    self.energy = float(splt_rest[0]);   self.apical_loop_length= int(splt_rest[1]);  self.upper_stem_length= int(splt_rest[2]);  self.uga_loop_spacing= int(splt_rest[3]);

def load_bsecis(bseblastian_outfile):
  """ Loads all secis instances found in a seblastian_outfile  (the one normally with extension output.all_secis) """
  out_list=[]; current_id=1
  file_h=open(bseblastian_outfile)
  line=file_h.readline()
  while line:
    while line and not line.startswith('------- bSECIS element'):        line=file_h.readline()
    if line:       
      #grade= line.split('grade:')[1][0]    #------- SECIS element -- id:1 -- grade:B ----
      line=file_h.readline();       splt=line.strip().split()   # Chromosome = chrom Strand = - Positions = 271 - 195
      chrom= splt[2];     strand= splt[5]
      pos1=int(splt[8]);  pos2=int(splt[10])
      line=file_h.readline(); splt=line.strip().split()        # Model = z090-4   Infernal score = 19.64     
      model_name=splt[2];           infernal_score= float(splt[-1])
      line=file_h.readline(); splt=line.strip().split()        #  Apical loop length = 3   Upper stem length = 7   UGA-loop spacing = 19
      apical_loop_length= int(splt[4]); upper_stem_length=int(splt[9]); uga_loop_spacing=splt[13]
      line=file_h.readline();              # Free Energy = -13.66
      energy= float(line.strip().split()[-1])
## example of next lines:
#
#         ..(((((((.......(((((((((((((((..................)))))))))))))))...........)))))))..
# model   uucggucaugcgUuaAUGAcGgccuggcccuAAAcC---cuuuuuggggcgggccaggcCuGAUGuu--u-uuucaugaccggc
#         .(((((((((......{[[{(((({(((((------((..........)))))))}))))}]]}.......-..))))))))).
# target  UGCCGAAUUAUUUAAAUGAAGAUCAAGUGA------UGGAGUUGGUCGCAUCGUUCGGUCAGAUCAAGGCG-UUUAAUUUGGUC
#                         ****                                        ****                    
      line=file_h.readline();                         line=file_h.readline();  
      model_aligned_ss =line.strip();                 line=file_h.readline();  
      model_aligned_seq=line.strip().split()[-1] ;    line=file_h.readline();  
      aligned_ss =line.strip();                       line=file_h.readline();  
      aligned_seq=line.strip().split()[-1] ;          line=file_h.readline();  
      annotation_line= line[9:-1];                    line=file_h.readline()
      s= bsecis()
      s.chromosome=chrom;   s.strand=strand;   s.energy=energy; 
      s.model_name=model_name; s.infernal_score=infernal_score
      s.apical_loop_length=apical_loop_length; s.upper_stem_length=upper_stem_length; s.uga_loop_spacing=uga_loop_spacing
      
      s.add_exon( min([pos1, pos2]),   max([pos1, pos2])  )
      s.ali=alignment();    s.ali.add('q', model_aligned_seq);     s.ali.add('t', aligned_seq); 
      s.ss_ali=alignment(); s.ss_ali.add('q', model_aligned_ss);   s.ss_ali.add('t', aligned_ss); 
      out_list.append(s)
    line=file_h.readline()    
  return out_list
    
class secis(p2g_feature):
  """ Class for SECISearch3/Seblastian secis predictions. Instanciated just with method load_secis (see below) """
  included_in_output=True
  included_in_gff   =True
  program='SECISearch3'
  gff_add=['grade', 'M:infernal_score', 'energy', 'type'] #'M:covels_score',

  def seq(self):    return nogap(   self.ali.seq_of('t'))
  def ss(self):     return nogap(self.ss_ali.seq_of('t'))
  def pretty_alignment(self):
    o= '         '+self.ss_ali.seq_of('q')+'\n'
    o+=' model   '+self.ali.seq_of('q')+'\n'
    o+='         '+self.ss_ali.seq_of('t')+'\n'
    o+=' target  '+self.ali.seq_of('t')+'\n'
    o+='         '+self.annotation_line
    return o
  def summary(self):
    o='-SECIS: '+self.id+'     grade:'+self.grade+'\n'
    o+=' Chromosome:'+self.chromosome+' Strand:'+self.strand+'\n'
    o+=' Positions:'+self.positions_summary()+' '
    d=self.distance_from_cds() 
    if d: o+='    CDS-Distance:'+str(d)+' '
    d=self.distance_from_sec_uga() 
    if d: o+='    Sec-Distance:'+str(d)+' '
    o+='\n'+self.found_by_line  # has terminal \n
    o+=' Free Energy: '+str(self.energy)+'\n'
    if self.type: o+=' Type: '+str(self.type)+'\n'
    o+='\n'+self.pretty_alignment()
    return o
  __str__=  summary
  __repr__= summary
  output=summary

  def distance_from_sec_uga(self):
    """ Return the distance from the (last) sec uga (not counting predicted introns) """
    try: 
      self.parent
      u_pos_in_target_seq = self.parent.protein().rfind('U')
      if u_pos_in_target_seq!=-1:
        return ( len(self.parent.protein())-u_pos_in_target_seq-1 )*3 +self.distance_from_cds()
    except: sys.exc_clear()
  def distance_from_cds(self):
    """ Return the distance from the last codon of the coding sequence predicted by selenoprofiles. 0 would mean that the secis begins right after the last codon. """
    try:
      self.parent
      return self.parent.is_upstream_of(self)-1
    except: sys.exc_clear()

  def dump_text(self):
    """see above, function to save the object in a db """
    o= self.header()
    o+='#'+self.grade+' '+str(self.energy)+' '+str(self.type)
    o+='#'+self.ali.seq_of('q')+' '+self.ali.seq_of('t')+' '+self.ss_ali.seq_of('q')+' '+self.ss_ali.seq_of('t')
    o+='#'+self.found_by_line+'#'+self.annotation_line
    return o

  def load_dumped_text(self, txt):
    """see above, functino to load object from a db"""
    splt=txt.split('#')
    self.load_from_header(splt[0])
    splt_1=splt[1].split()
    self.grade=         splt_1[0] 
    self.energy=  float(splt_1[1]) 
    if len( splt_1 )> 2: self.type=splt_1[2]
    else:                self.type='None'

    splt_ali= splt[2].split()
    self.ali=alignment();     self.ali.add('q', splt_ali[0]    );   self.ali.add('t', splt_ali[1]); 
    self.ss_ali=alignment();  self.ss_ali.add('q', splt_ali[2] );   self.ss_ali.add('t', splt_ali[3]); 
    self.found_by_line=splt[3]
    self.annotation_line=splt[4]
    
  def infernal_score(self): 
    """ Returns a float with the infernal score of the prediction. Note: it will return None if the score is not found in the .found_by_line attribute"""
    if "Infernal" in  self.found_by_line: 
      score= float(   self.found_by_line.split('Infernal')[1].split(')')[0].split('(')[1].split('=')[1]   ) # if crash: wrong "Found by" line in prediction!
    else: score=None
    return score  

  def covels_score(self): 
    """ Returns a float with the covels score of the prediction. Note: it will return None if the score is not found in the .found_by_line attribute"""
    if "Covels" in  self.found_by_line: 
      score= float(   self.found_by_line.split('Covels')[1].split(')')[0].split('(')[1].split('=')[1]   ) # if crash: wrong "Found by" line in prediction!
    else: score=None
    return score  
    
def load_secis(seblastian_outfile):
  """ Loads all secis instances found in a seblastian_outfile  (the one normally with extension output.all_secis) """
  out_list=[]; current_id=1
  file_h=open(seblastian_outfile)
  line=file_h.readline()
  while line:
    while line and not line.startswith('------- SECIS element'):        line=file_h.readline()
    if line:       
      grade= line.split('grade:')[1][0]    #------- SECIS element -- id:1 -- grade:B ----
      line=file_h.readline()               # Chromosome = chrom Strand = - Positions = 271 - 195
      chrom= line.split()[2];     strand= line.split()[5]
      pos1=int(line.split()[8]);       pos2=int(line.strip().split()[10])
      line=file_h.readline();     found_by_line=line  # Found by:  Infernal (score=13.72) Covels (score=4.0)
      line=file_h.readline();              # Stem2 length = 9   Apical loop length = 14   Insertion in stem2 = 0, 0
      line=file_h.readline();              # Free Energy = -13.66
      energy= float(line.strip().split()[-1])
      line=file_h.readline();
      if line.strip():                     # Type = 2   BUT this line may be absent from old output files
        secis_type=line.strip().split()[-1]
        line=file_h.readline();
      else: secis_type='None'
## example of next lines:
#
#         ..(((((((.......(((((((((((((((..................)))))))))))))))...........)))))))..
# model   uucggucaugcgUuaAUGAcGgccuggcccuAAAcC---cuuuuuggggcgggccaggcCuGAUGuu--u-uuucaugaccggc
#         .(((((((((......{[[{(((({(((((------((..........)))))))}))))}]]}.......-..))))))))).
# target  UGCCGAAUUAUUUAAAUGAAGAUCAAGUGA------UGGAGUUGGUCGCAUCGUUCGGUCAGAUCAAGGCG-UUUAAUUUGGUC
#                         ****                                        ****                    
      line=file_h.readline(); 
      model_aligned_ss =line.strip();                 line=file_h.readline();  
      model_aligned_seq=line.strip().split()[-1] ;    line=file_h.readline();  
      aligned_ss =line.strip();                       line=file_h.readline();  
      aligned_seq=line.strip().split()[-1] ;          line=file_h.readline();  
      annotation_line= line[9:-1];                    line=file_h.readline()
      s= secis()
      s.chromosome=chrom;   s.strand=strand;  s.found_by_line=found_by_line;   s.energy=energy; s.grade=grade
      s.add_exon( min([pos1, pos2]),   max([pos1, pos2])  )
      s.ali=alignment();    s.ali.add('q', model_aligned_seq);     s.ali.add('t', aligned_seq); 
      s.ss_ali=alignment(); s.ss_ali.add('q', model_aligned_ss);   s.ss_ali.add('t', aligned_ss); 
      s.annotation_line=annotation_line
      s.type=secis_type
      out_list.append(s)
    line=file_h.readline()    
  return out_list

def Secisearch3(p2g, three_prime_length=-1, silent=False, full=False):
  """ Performs a complete secisearch3 with seblastian. Parse secis and restore their coordinates cosindering their parent gene so that they are absolute. Add the secis to the .features list of the p2g object. If silent!=True, it prints a message for every SECIS found.   If three_prime_length is -1, the value specified in the main config file (or in the command line) is used.
  If full==True, all methods are run with seblastian, to ensure maximal sensitivity """
  #cutting 3' UTR
  three_prime_file=temp_folder+'three_prime_for_secisearch.fa'
  three_prime_text=p2g.three_prime(length=three_prime_length)  
  
  if three_prime_text:
    three_prime_region_width=int( three_prime_text.split('_bp_downstream')[0].split('region_of_')[1] )
    three_prime_region=p2g.downstream(distance=0, region_length=three_prime_region_width)
    three_prime_text= ">three_prime_of_id:"+str(p2g.id)+'\n'+ join( three_prime_text.split('\n')[1:], '\n')
    write_to_file(three_prime_text, three_prime_file)
    #running SS3
    command= 'Seblastian.py '+three_prime_file+' -SS -infernal_no_mpi -type -c -temp '+temp_folder+'secisearch3'
    if full: command+=' -m all '
    try:         b=bash(command); assert not b[0]
    except:      raise notracebackException, "ERROR SECISearch3 (Seblastian) is not installed or didn't work for some reason, here below you find its error message. If you wish not to run SECISearch3, comment with a starting \"#\" the line with ACTION.post_filtering.secisearch in the selenoprofiles configuration file.\nCOMMAND: "+command+' ; ERROR: '+b[1]
    ss3_out_file=temp_folder+'three_prime_for_secisearch.output.all_secis'
    # loading results
    all_secises= load_secis(ss3_out_file)
    if all_secises:
      for secis_element_index in range(len(all_secises)): 
        all_secises[secis_element_index].restore_absolute_coordinates(parent_gene=three_prime_region)
        all_secises[secis_element_index].species=p2g.species;  all_secises[secis_element_index].target=p2g.target
        all_secises[secis_element_index].id=p2g.output_id()+'.secis'+str(secis_element_index+1)
        all_secises[secis_element_index].parent=p2g
        p2g.features.append(all_secises[secis_element_index])
        if not silent:   write('SECISearch3 found a grade '+all_secises[secis_element_index].grade+' eukaryotic SECIS for '+p2g.output_id(), 1)
        #all_euk_secises[secis_element_index].parent=p2g
    bbash('rm -r '+temp_folder+'three_prime_for_secisearch.* '+temp_folder+'secisearch3') # * is to remove secis files as well  
    if all_secises: return all_secises
    else: return []
  return None

def bSecisearch(p2g, silent=False, full=False):
  """ Performs a complete bsecisearch searche with bseblastian (crash if not installed). Parse secis and restore their coordinates cosindering their parent gene so that they are absolute. Add the secis to the .features list of the p2g object. If silent!=True, it prints a message for every bSECIS found.    """
  #cutting 3' UTR
  cds_seq= p2g.cds(); prot_seq=p2g.protein()  
  sec_ugas=[]  #pos codon based, 0 based
  for codon_index in range( len(cds_seq)/3 ):
    if prot_seq[codon_index]=='U': # and cds_seq[codon_index*3:codon_index*3+3]   =='TGA'
      sec_ugas.append(  codon_index  )
  if sec_ugas:
    ### cutting from 5 nts before the first uga, to 100 nts after the last one  (generally there's only one)
    pos_start_subsequence =   sec_ugas[0]*3  +1   -5     ## 1 based, nt based
    pos_end_subsequence   =   sec_ugas[-1]*3  +1   +2 +100    ## 1 based, nt based. Including 100 nts after the UGA codon 
    length_subsequence    =   pos_end_subsequence - pos_start_subsequence +1 

    constrained_start   = max([pos_start_subsequence, 1])
    available_length_in_p2g = p2g.length() - constrained_start + 1
    constrained_length = min ([available_length_in_p2g, length_subsequence])
    subseq_gene_obj =p2g.subseq( constrained_start, constrained_length, minimal=True )
    remainder_up  = constrained_start  -  pos_start_subsequence
    remainder_down= length_subsequence -  constrained_length
    if   subseq_gene_obj.strand=='+':     subseq_gene_obj.extend(left=remainder_up,   right=remainder_down, inplace=True)
    elif subseq_gene_obj.strand=='-':     subseq_gene_obj.extend(left=remainder_down, right=remainder_up, inplace=True)

    chromosome_length= get_MMlib_var('chromosome_lengths')[p2g.chromosome]
    prev_boundaries=subseq_gene_obj.boundaries()
    if subseq_gene_obj.check_boundaries(chromosome_length):  #this checks and modify inplace                                                                     
      # out of bounds                                                                                                                                            
      if p2g.strand=='+':
        pos_start_subsequence+=  subseq_gene_obj.boundaries()[0] - prev_boundaries[0]
        length_subsequence-=     (prev_boundaries[1]-subseq_gene_obj.boundaries()[1] +subseq_gene_obj.boundaries()[0]-prev_boundaries[0] )
      elif p2g.strand=='-':
        pos_start_subsequence+=  prev_boundaries[1] - subseq_gene_obj.boundaries()[1]
        length_subsequence-=     (subseq_gene_obj.boundaries()[0] - prev_boundaries[0] +prev_boundaries[1] - subseq_gene_obj.boundaries()[1])

    subseq=p2g.subsequence( pos_start_subsequence, length_subsequence)
    
    target_temp_file=temp_folder+'sequence_for_bsecisearch.fa'
    write_to_file(">region_from:"+str(p2g.id)+' '+subseq_gene_obj.header(no_id=True) +'\n'+subseq,      target_temp_file )
    command= 'bSeblastian.py '+target_temp_file+' -SS -infernal_no_mpi -c -temp '+temp_folder+'bsecisearch'
    try:         b=bash(command); assert not b[0]
    except: raise notracebackException, "ERROR bSECISearch (bSeblastian) is not installed or didn't work for some reason, here below you find its error message. If you wish not to run bSECISearch\
3, comment with a starting \"#\" the line with ACTION.post_filtering.bsecisearch in the selenoprofiles configuration file.\nCOMMAND: "+command+' ; ERROR: '+b[1]
    #raw_input('...'+temp_folder)
    bsecis_out_file=temp_folder+'sequence_for_bsecisearch.output.all_secis'
    all_secises= load_bsecis(bsecis_out_file)
    if all_secises:
      for secis_element_index in range(len(all_secises)): 
        all_secises[secis_element_index].restore_absolute_coordinates(parent_gene=subseq_gene_obj)
        all_secises[secis_element_index].species=p2g.species;  all_secises[secis_element_index].target=p2g.target
        all_secises[secis_element_index].id=p2g.output_id()+'.bsecis'+str(secis_element_index+1)
        all_secises[secis_element_index].parent=p2g
        p2g.features.append(all_secises[secis_element_index])
        if not silent:   write('bSeblastian found a bacterial SECIS for '+p2g.output_id(), 1)
        #all_euk_secises[secis_element_index].parent=p2g
    bbash('rm -r '+temp_folder+'sequence_for_bsecisearch.* '+temp_folder+'bsecisearch') # * is to remove secis files as well  
    if all_secises: return all_secises
    else: return []
  return None

def gene_ontology(obo_file=''):
  """ Returns a GeneOntologyNX class, loaded from the obo_file provided or from the one defined in the configuration file as "GO_obo_file".
  This class is from http://gitorious.org/annotation/annotation/trees/master
  The two methods used here are: 
  get_term_by_id(go_id)         -> returns a GOTerm class object 
  get_all_parents_ids(go_id)    -> returns a list of strings (GO_ids)   
  """

  if not 'gene_ontology_data' in globals():
    #checking if necessary module is loaded
    try:          oboparser
    except:       raise notracebackException, "gene_ontology ERROR cannot find the modules necessary to use the gene ontology features. Please see the installation script."
    #checking if obo_file is defined (and valid) either as argument of this function or in the main config file
    if opt['GO_obo_file'] and is_file(opt['GO_obo_file']):      obo_file=opt['GO_obo_file']
    else:      raise notracebackException, "gene_ontology ERROR cannot load obo file, not defined or not found. Add a line like ' GO_obo_file  = PATH_TO_FILE '   in the main configuration file."
    global gene_ontology_data
    parserO = oboparser.Parser (open(obo_file))
    gene_ontology_data = parserO.parse()
  return gene_ontology_data

def set_selenoprofiles_var(varname, value):
  """ Utility to set a variable inside the selenoprofiles module from an external python program which is importing it"""
  globals()[varname]=value
def get_selenoprofiles_var(varname):
  """ Utility to get a variable inside the selenoprofiles module from an external python program which is importing it"""
  return globals()[varname]

def uniq_id_for_file(ffile):  return replace( fileid_for_temp_folder(ffile), '.', '_')

####
class notracebackException(Exception):
  """ When this exception are raised, the traceback is not printed, just the message it contains"""

class lazycomputationException(Exception):
  """ This exception is raised when lazy computation cannot be applied, as a flag for: compute it again!"""

class skipprofileException(Exception):
  """ When this exception is raised, the pipelines print the error and shifts to the next profile"""

#######################################################################################################################################
def close_program():
  """ Utility to perform operation at the end of the computation, even in case the pipeline crashes"""
  global allowed_output_formats
  if 'temp_folder' in globals():
    if 'opt' in globals() and opt['debug']: 
      if sys.exc_info()[1]!= None:
        printerr('error preview: '+str(sys.exc_info()[1]), 1)
        raw_input('some error occured. check temp folder '+temp_folder)
    try:     bash('rm -r '+temp_folder)
    except:  pass
  for keyword in allowed_output_formats:
    handler_opened='output_'+keyword+'_file_h'
    if handler_opened in globals():
      if not globals()[handler_opened].closed: globals()[handler_opened].close()

  if sys.exc_info()[0]: #an exception was raised. let's print the information using printerr, which puts it in the logfile as well, if any.   
    if issubclass(sys.exc_info()[0], notracebackException):      printerr( sys.exc_info()[1], 1)
    elif issubclass(sys.exc_info()[0], SystemExit):      pass
    else:                                                                  printerr('ERROR '+ traceback.format_exc( sys.exc_info()[2]) , 1)

  if 'log_file' in globals(): log_file.close()
  try:
    if get_MMlib_var('printed_rchar'):       printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:    pass

if __name__ == "__main__":
  try:
    main()
    close_program()  
  except:
    close_program()


