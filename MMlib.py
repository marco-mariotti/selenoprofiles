import os
import re
import sys
from string import *
import commands
import subprocess
try:      import cPickle as pickle
except:   import pickle
from copy import copy, deepcopy
from math import log as math_log, sqrt
import random
NUMBERS='0123456789'  
DNA_LETT='ACGT'
SS_LETT='0SETBH'
AA_LETT='ACDEFGHIKLMNPQRSTVWYUOB'  #  including selenocysteine
AA_LETT_STRICT='ACDEFGHIKLMNPQRSTVWY'
RNA_LETT='ACGU'
three_letter_codon_diz={ 'Ala':'A' ,  'Cys':'C' ,  'Asp':'D' ,  'Glu':'E' ,  'Phe':'F' ,  'Gly':'G' ,  'His':'H' ,  'Ile':'I' ,  'Lys':'K' ,  'Leu':'L' ,  'Met':'M' ,  'Asn':'N' ,  'Pro':'P' ,  'Gln':'Q' ,  'Arg':'R' ,  'Ser':'S' ,  'Thr':'T' ,  'Val':'V' ,  'Trp':'W' ,  'Tyr':'Y', '***':'*', 'Unk':'X' , '<->':'-', '---':'-','Asx':'B', 'SeC':'U', 'Zed':'Z', 'SeC(e)':'U' }
one_to_three_letter_codon_diz={ 'A':'Ala' ,  'C':'Cys' ,  'D':'Asp' ,  'E':'Glu' ,  'F':'Phe' ,  'G':'Gly' ,  'H':'His' ,  'I':'Ile' ,  'K':'Lys' ,  'L':'Leu' ,  'M':'Met' ,  'N':'Asn' ,  'P':'Pro' ,  'Q':'Gln' ,  'R':'Arg' ,  'S':'Ser' ,  'T':'Thr' ,  'V':'Val' ,  'W':'Trp' ,  'Y':'Tyr', '*':'***', 'X':'Unk' , '-':'---' , 'B':'Asx' }
STOP_CODONS=    ['UAA', 'UAG', 'UGA', 'URA', 'UAR'] 
STOP_CODONS_DNA=['TAA', 'TAG', 'TGA', 'TRA', 'TAR']
one_letter_to_three_letter_aa_diz={'A':'alanine', 
'R':'arginine', 'N':'asparagine','D':'aspartic_acid','C':'cysteine','E':'glutamic_acid','Q':'glutamine','G':'glycine','H':'histidine','I':'isoleucine','L':'leucine','K':'lysine','M':'methionine','F':'phenylalanine','P':'proline','S':'serine','T':'threonine','W':'tryptophan','Y':'tyrosine','V':'valine','U':'selenocysteine'}



try: from numpy import average, std as std_deviation
except: 
  sys.exc_clear()
  def average(ilist):    return sum(ilist)/float(len(ilist))
  def std_deviation(ilist):
    a=average(ilist)
    return   sqrt(  sum(   [pow( v-a, 2)  for v in ilist]  )/len(ilist)  )

def set_temp_folder(folder):  set_MMlib_var('temp_folder', folder)
def set_split_folder(folder): set_MMlib_var('split_folder', folder)
def get_temp_folder():        return get_MMlib_var('temp_folder')
def get_split_folder():       return get_MMlib_var('split_folder')
 
def set_local_folders(temp='/tmp'):
  """ Used in ipython to quickly set the environment for fetching chromosomes and other stuff that required temp files"""
  try: assert is_directory(opt['temp'])
  except:  opt['temp']=temp
  temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_temp_folder(temp_folder)        
  split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_split_folder(split_folder)                                              

def mute(also_stderr=False):
  """ Turns off any output to stdout (to stderr as well if option is True). To go back to normal , use unmute()"""
  sys.stdout = open(os.devnull, "w")
  if also_stderr:   sys.stderr = open(os.devnull, "w")

def unmute():
  sys.stdout = sys.__stdout__;  sys.stderr = sys.__stderr__

def phylome_connector():
  """  ete2 connector to PhylomeDB. allows retrieving data such as trees, alignments etc.   Examples: (calling p the object returned by this function)
  p.get_algs("Phy0005QIK_DROME", 8)
  p.get_tree("Phy0005QIK_DROME", 8, best_tree = True)

  """
  import ete2
  p = ete2.PhylomeDB3Connector(host = "phylomedb.org", user = "phyReader", passwd = "phyd10.-Reader", db = 'phylomedb_3')
  p._algs = "alignment"
  p._trees = "tree"
  p._phylomes = "phylome"
  p._phy_content = "phylome_content"
  return p

def bash(command, print_it=0):
  """Utility to run bash commands. a tuple (exit_status, message) is returned, where message includes both std input and stderr. If argument print_it==1 or the variable print_commands is defined in MMlib, the command is printed before execution. If variable bin_folder is defined in MMlib, this folder is added to bash $PATH before running the command. """
  if 'print_commands' in globals() and print_commands: print_it=1
  if print_it:    write(command, 1)
  if 'bin_folder' in globals(): 
    if not bin_folder  == os.environ['PATH'].split(':')[0]:      os.environ['PATH']=str(bin_folder)+':'+os.environ['PATH']
  b1, b2= commands.getstatusoutput(command)
  return [b1, b2]
  
def bbash(command, print_it=0, dont_die=0):
  """Utility to run bash commands. A string is returned, where including both std input and stderr. If the exit status of the command is different than 0, it is assumed something went wrong, so an exception is raised indicating the command and the output. If argument dont_die==1, no exception is raised and output is returned as normal.
   If argument print_it==1 or the variable print_commands is defined in MMlib, the command is printed before execution. If variable bin_folder is defined in MMlib, this folder is added to bash $PATH before running the command. """
  if 'print_commands' in globals() and print_commands: print_it=1
  if print_it:    write(command, 1)
  cmnd=command
  if 'bin_folder' in globals(): 
    if not bin_folder  == os.environ['PATH'].split(':')[0]:      os.environ['PATH']=str(bin_folder)+':'+os.environ['PATH']
  bb=commands.getstatusoutput(cmnd)
  if bb[0]!=0 and not dont_die:     raise Exception, 'COMMAND: ' + command+' ERROR: "'+bb[1]+' "'
  else:                             return bb[1]
  
def bash_pipe(cmnd, print_it=0, return_popen=0, stdin=None):
  """ Open a filehandler which reads from a pipe opened in bash, given a command. Useful when you want to read one line at the time. If variable bin_folder is defined in MMlib, this folder is added to bash $PATH before running the command.
  stdin can be used to input large chunks of data, through a filehandler. valid formats are: stdin=existing_write_filehandler, or stdin='PIPE' (equivalent to subprocess.PIPE); in this last case you want to use return_popen to be able to access the handler (with this option, the popen object is returned instead of its stdout filehandler), like this:
  p=bash_pipe('command', stdin='PIPE', return_popen=True)
  print >> p.stdin, 'input_lines!'   #repeat as many times as needed
  p.stdin.close()   # important! as most programs wait for a EOF signal to give output, you need to close the input filehandler before reading the lines of output
  p.stdout.readline()  # --> now you can read the output lines from this handler
  
  """
  if stdin=='PIPE': stdin=subprocess.PIPE 
  if 'print_commands' in globals() and print_commands: print_it=1
  if 'bin_folder' in globals(): 
    if not bin_folder  == os.environ['PATH'].split(':')[0]:      os.environ['PATH']=str(bin_folder)+':'+os.environ['PATH']    
  if print_it: write(cmnd, 1)
  s=subprocess.Popen(cmnd.split(), stdout=subprocess.PIPE, stdin=stdin, env=os.environ)
  if return_popen: return s
  else:    return s.stdout

def md5_executable():
  b=bash('echo | md5sum')
  if not b[0]: return 'md5sum'  ## command found, no error
  b=bash('echo | md5')
  if not b[0]: return 'md5'  ## command found, no error
  else: raise Exception, "ERROR neither md5sum  or  md5  found on this system!" 

def checksum(ffile, is_text=False):
  """ Returns the checksum for the file in input """  
  if is_text:
    pipe=bash_pipe(md5_executable()+' ', return_popen=1, stdin='PIPE')
    print >> pipe.stdin,  ffile
    pipe.stdin.close()
    m=    pipe.stdout.readline().split()[0]
  else:
    m=bbash(md5_executable()+' '+ffile).split()[0] 
  return m

    
def Folder(string):
  if not string:
    return string
  ff=string+'/'*int(string[-1]!='/')
  cmnd='mkdir '+ff
  bb=bash(cmnd)
  return ff
  
def random_folder(parent_folder='', dont_create_it=0):
  if parent_folder:
    parent_folder = Folder(parent_folder) #checking or creating parent folder. IF CRASHED: do you have writing privileges here?
  a=parent_folder+ bash("date +%F%T%N | "+md5_executable()+" | cut -c 1-32")[1] #creating a random named folder inside the parent folder
  if dont_create_it:
    return a+'/'
  a=Folder(a)
  if not bash('cd '+a+' ; cd ..')[0] :
    return a
  else:
    return 'some_error_creating_random_folder.hopefully_there_is_noooo_file_named_like_this_when_you_try_to_delete_it'

temp_folder=random_folder('/tmp', 1)
split_folder=Folder('/tmp')

def set_MMlib_var(varname, value):  globals()[varname]=value
def get_MMlib_var(varname):         return globals()[varname]

def is_significant(pvalue):         return pvalue<opt['alpha']

printed_rchar=0

# trans={};
# trans ['GCA'] = "A";     trans ['GCC'] = "A";     trans ['GCG'] = "A";     trans ['GCT'] = "A";     trans ['TGC'] = "C";     trans ['TGT'] = "C";
# trans ['GAC'] = "D";     trans ['GAT'] = "D";     trans ['GAA'] = "E";     trans ['GAG'] = "E";     trans ['TTC'] = "F";     trans ['TTT'] = "F"; 
# trans ['GGA'] = "G";     trans ['GGC'] = "G";     trans ['GGG'] = "G";     trans ['GGT'] = "G";     trans ['CAC'] = "H";     trans ['CAT'] = "H";  
# trans ['ATA'] = "I";     trans ['ATC'] = "I";     trans ['ATT'] = "I";     trans ['AAA'] = "K";     trans ['AAG'] = "K";     trans ['TTA'] = "L";   
# trans ['TTG'] = "L";     trans ['CTA'] = "L";     trans ['CTC'] = "L";     trans ['CTG'] = "L";     trans ['CTT'] = "L";     trans ['ATG'] = "M";   
# trans ['AAC'] = "N";     trans ['AAT'] = "N";     trans ['CCA'] = "P";     trans ['CCC'] = "P";     trans ['CCG'] = "P";     trans ['CCT'] = "P";  
# trans ['CAA'] = "Q";     trans ['CAG'] = "Q";     trans ['AGA'] = "R";     trans ['AGG'] = "R";     trans ['CGA'] = "R";     trans ['CGC'] = "R";  
# trans ['CGG'] = "R";     trans ['CGT'] = "R";     trans ['AGC'] = "S";     trans ['AGT'] = "S";     trans ['TCA'] = "S";     trans ['TCC'] = "S";  
# trans ['TCG'] = "S";     trans ['TCT'] = "S";     trans ['ACA'] = "T";     trans ['ACC'] = "T";     trans ['ACG'] = "T";     trans ['ACT'] = "T";  
# trans ['GTA'] = "V";     trans ['GTC'] = "V";     trans ['GTG'] = "V";     trans ['GTT'] = "V";     trans ['TGG'] = "W";     trans ['TAC'] = "Y";   
# trans ['TAT'] = "Y";    # trans ['taa'] = "!";     trans ['tag'] = "#";     trans ['tga'] = "@";
# trans ['TAA'] = "*";     trans ['TAG'] = "*";     trans ['TGA'] = "*";
# trans ['---'] = "-";
## std: FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG

# build alternative genetic code translation tables based on NCBI codes
genetic_codes={}
genetic_codes_AAs={    1:'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                       2:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
                       3:'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                       4:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                       5:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
                       6:'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                       9:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
                      10:'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      11:'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      12:'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      13:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
                      14:'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
                      16:'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      21:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
                      22:'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      23:'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      24:'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
                      25:'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      26:'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      27:'FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      28:'FFLLSSSSYYQQCCWWLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      29:'FFLLSSSSYYYYCC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      30:'FFLLSSSSYYEECC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
                      31:'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'}
for gc_code in genetic_codes_AAs:
  i=-1
  genetic_codes[gc_code]={'---':'-'}
  for a in 'TCAG':
    for b in 'TCAG':
      for c in 'TCAG':
        i+=1
        genetic_codes[gc_code][ a+b+c ]=genetic_codes_AAs[gc_code][i]
trans=genetic_codes[1]

retrotrans={}
for codon in trans:  retrotrans.setdefault( trans[codon], []  ).append(codon)
for aa in retrotrans: retrotrans[aa].sort()

species_code_file="/home/mmariotti/software/selenoprofiles/libraries/species_codes.tab"
genome_config="/home/mmariotti/software/selenoprofiles/libraries/genome.config"

def contain_chars(string, to_check=uppercase+lowercase):
        for char in string:
            if char in to_check:
                return 1
        return 0
def is_number(string, mode='int'):
  if mode=='float':
    try:
      float(string)
      return True
    except ValueError:
      return False
  else:  #mode =int
    try:
      int(string)
      return True
    except ValueError:
      return False
def is_option(s):
        return (s[0]=='-' and contain_chars(s[1:]))
def option_value(value):
  """value=string; returns value changed into the appropriate type.
  """
#  if value.startswith('[') and value.endswith(']'):     
#    write(value, 1, how='yellow')
#    return [option_value(x) for x in value[1:-1].split(', ') ]
  if  is_number(value):                      return int(value)
  elif is_number(value, 'float'):            return float(value)
  elif value=='None':                        return None
  else:
    if value and value[0]==value[-1] and value[0] in ['"', "'"] and len(value)>=2:   value=value[1:-1]
  return value

def update_opt(new_opt, current_opt):  
  """sometimes it is useful to read options from command line, then manipulate them before using them. In this case, it is worth doing to update the opt object and sys.argv.
     Each option in new_opt (key) goes to replace the one in current_opt with the same name. Sys.argv is also updated.DONT KNOW ABOUT BOOLEAN VALUES
  """
  for k in new_opt:
    current_opt[k]=new_opt[k]
    for c in range(len(sys.argv)):
      if sys.argv[c]=='-'+k:
        sys.argv[c+1]=new_opt[k]
  return current_opt


def fill_option(default_opt, option, dont_add_to_argv=0):   #fill the keys of option (dictionary like {option: value}) that are not present, with default values specified in default_opt
  for key in default_opt:
    if not option.has_key(key):
      option[key]=default_opt[key]
      if not dont_add_to_argv:
        sys.argv.extend(['-'+key, str(default_opt[key])])
  return option

class options(dict):
  """ """
  def __init__(self, dict_in={}, synonyms={}):
    super(options, self).__init__()
    for k in dict_in:
      self[k]=dict_in[k]
    self.set_synonyms(synonyms)
    
  def __getitem__(self, name):
    if self.has_key(name):                        return super(options, self).__getitem__(name)
    elif self['__synonyms__'].has_key(name):      return self[ self['__synonyms__'][name] ]
    elif name in self['__synonyms__'].values():
      for k in self['__synonyms__']:
        if name== self['__synonyms__'][k]:
          return self[k]
    else:
      return None

  def set_synonyms(self, syno_hash):
    self['__synonyms__']=syno_hash
  def add_synonym(synonym, key):
    self['__synonyms__'][synonym]=key
  def synonyms(self):
    return self['__synonyms__']

def command_line(default_opt, command_line_help='Command line usage:...', default_order='io', synonyms={}, dont_die=0, silent=0, dont_add_to_argv=0, nowarning=0, strict=0, tolerated_regexp=[], advanced={}):
    import sys
    def command_line_option():
      """this is a utility that read sys.argv and returns a dictionary like {option: value}
          options in command line must be preceded by -, and option can't consist only of number ;
          if after an option there is no value (for example there is another option), the value is set to 1.
          Numberic values are converted to their types; floating point numbers are discriminated from integers if they contains a dot '.' 
          ex:
          [python] prog_name.py -c 12 -verbose -gap -0.45    ----> returns {'c': 12, 'g': -0.45, 'verbose':1}

    The default_opt contains all default values for compulsory arguments. These are taken from here if they are not find in the commandline. Also, sys.argv is enriched with these value, therefore if the command_line function is called again (even with specifying def_opt), the resulting returned opt  will be the same.
      """
      option=options()
      llist=sys.argv[1:]+['-EnD!!']
#    write_to_file(str( sys.argv), 'btemp')
      c=0
    #flushing empty values in argv
      while c<len(llist):
        if not llist[c]:
          llist.pop(c)
        else:
          c+=1
        while llist[0]!='-EnD!!':
          first=str(llist.pop(0))
              #option=option_value(value)
          if first[0]=='-':
            if  is_option(llist[0]): #the next object in the row  is a option. So this option is simply set to True (1)
              option[first[1:]]=1              
            else:
              value=llist.pop(0)
              option[first[1:]]=option_value(value)
        return option
    if silent: nowarning=1
    opt=command_line_option()
    if default_order=='*':
      obj_list=[]
      i=1
      while  i<len(sys.argv) and ( not  is_option(sys.argv[i]) )  :   #reading as syntax is  --> program_name  file1 file2 file3 .... fileN -o 1 -n 7
        obj_list.append(  sys.argv[i]  )
        i+=1
      if not obj_list and i<len(sys.argv): #it may be that the syntax used in the command line is also: program_name  -o 1 -n 7 file1 file2 file3 ....
        i=1
        while  i<len(sys.argv):
          if not is_option(sys.argv[i]) and  (not is_option(sys.argv[i-1])):
            obj_list.append(  sys.argv[i]  )
          i+=1
      opt['*']=obj_list
    else:
      i=1
      while  i<len(sys.argv) and i<=len(default_order) and ( not  is_option(sys.argv[i]) )  :
          opt[ default_order[i-1] ]=option_value(sys.argv[i])
          i+=1

    synonyms['help']='h';     synonyms['-help']='h'
    opt.set_synonyms(synonyms)
    special_options=['h', 'print_opt', 'print_option', 'print_options', 'config', 'bin_folder', '__synonyms__']
    for k in opt:
      if not default_opt.has_key(k) and not nowarning and not k in special_options and not k in opt.synonyms() and not  match_any_word(k, tolerated_regexp, ignore_case=0) :
        if strict!=0:
          if type(strict) in [int, bool]: e=Exception
          elif issubclass(strict, Exception): e=strict
          raise e, "command_line ERROR the option "+k+" is not present among the possible options. run with -print_opt to see the current option configuration"
        printerr('WARNING possible typo: the option '+k+' is not present among the possible options. run with -print_opt to see the current option configuration\n')

    opt=fill_option(default_opt, opt, dont_add_to_argv)
    for k in opt:
      if type(opt[k]) == str and opt[k].startswith('$('): opt[k]=bbash(opt[k][2:-1])

    #dealing with synonyms
    keyss=opt.keys()
    for k in keyss:
      if k in opt.synonyms():
        opt[opt.synonyms()[k]]=opt[k]
        del opt[k]

    #printing options to screen in case we have -print_opt
    if not silent and( opt.has_key('print_opt') or opt.has_key('print_options') or opt.has_key('print_option') ):
      write( "| ACTIVE OPTIONS:", 1)
      keys=opt.keys()
      keys.sort()
      for k in keys:
        a="| "+str(k)
        write( a+' '*(30-len(a))+': '+str(opt[k]), 1)
      write('', 1)
    #printing help message in case we have -h or --help
    if not silent and  (len(sys.argv)<2 or opt.has_key('h') ) :
      write(command_line_help, 1)
      if advanced and opt['h'] in advanced: 
        write(advanced[opt['h']], 1)
      if not dont_die:
        sys.exit()
        
    return opt

uniq_id=id

def lineage_string_to_abstract(lineage):
  """ lineage is a string which is usually returned by this program. This function condensate it keeping the most interesting classes. """
  splt=lineage.split('; ')
  if "Bacteria; " in lineage  :    return 'B; '+join(splt[2:min(5, len(splt))], '; ')
  elif 'Archaea; ' in lineage :    return 'A; '+join(splt[2:min(5, len(splt))], '; ')
  elif 'Eukaryota; ' in lineage:   
    out='E; '
    if 'Metazoa; ' in lineage:
      out+='M; '
      if 'Deuterostomia; ' in lineage:
        out+='Deuterostomia; '
        if 'Vertebrata; ' in lineage:
          out+='Vertebrata; '      
          if 'Mammalia; ' in lineage:           out+='Mammalia; '
          elif 'Sauropsida; ' in lineage:       out+='Sauropsida; '
          elif 'Amphibia; ' in lineage:         out+='Amphibia; '
          elif 'Actinopterygii; ' in lineage:   out+='Actinopterygii; '
          elif 'Chondrichthyes; ' in lineage:   out+='Chondrichthyes; '
        elif 'Tunicata; ' in lineage:          
          out+='Tunicata; '
          if 'Ascidiacea; ' in lineage:           out+='Ascidiacea; '         
        elif 'Branchiostomidae; ' in lineage:  out+='Branchiostomidae; '
        elif 'Echinodermata; ' in lineage:     out+='Echinodermata; '
      elif 'Protostomia; ' in lineage:
        out+='Protostomia; '
        if  'Arthropoda; ' in lineage:      
          out+='Arthropoda; '
          if 'Insecta; ' in lineage:        out+='Insecta; '
          elif 'Crustacea; ' in lineage:    out+='Crustacea; '
          elif 'Myrapoda; ' in lineage:     out+='Myrapoda; '
          elif 'Arachnida; ' in lineage:     out+='Arachnida; '
          elif 'Merostomata; ' in lineage:     out+='Merostomata; '
        elif 'Nematoda; ' in lineage:     out+='Nematoda; '
        elif 'Mollusca; ' in lineage:           
          out+='Mollusca; '
          if 'Gastropoda; ' in lineage: out+='Gastropoda; '
          elif 'Bivalvia; ' in lineage: out+='Bivalvia; '
        elif 'Annelida; ' in lineage:           out+='Annelida; '
        else:      out+=     lineage.split('Protostomia; ')[1].split(';')[0]+'; '
      else: #basal metazoan
        if 'Cnidaria; ' in lineage:        out+='Cnidaria; '
        elif 'Porifera; ' in lineage:      out+='Porifera; '          
        elif 'Ctenophora; ' in lineage:    out+='Ctenophora; '
        elif 'Placozoa; ' in lineage:      out+='Placozoa; '
        elif 'Platyhelminthes; ' in lineage: out+='Platyhelminthes; '

    else:      out+= join(splt[2:min(4, len(splt))], '; ')+'; '
    return out[:-2]
  else:      return join(splt[0:min(4, len(splt))], '; ')


def get_species_fullname(species_name):
    b=bash('egrep -w "'+species_name+'" '+species_code_file)
    if b[0]:
      raise Exception, "get_species_fullname ERROR: "+species_name+' not found'
    else:
      return b[1].split('\t')[1]
def get_species_code(species_name):
    b=bash('egrep -w "'+species_name+'" '+species_code_file)
    if b[0]:
      raise Exception, "get_species_code ERROR: "+species_name+' not found'
    else:
      return b[1].split('\t')[0]

def get_genome_file(species_fullname):
    b=bash('egrep "'+species_fullname+'.*=" '+genome_config)
    if b[0]:
      raise Exception, "get_genome_file ERROR: "+species_fullname+' not found'
    else:
      return del_white(b[1].split('=')[1])

def second_max(alist):
  """returns the second biggest number in a list. if it has only one element, that is returned. If it has none, it returns False
  """
  current_max='init'
  current_second_max='init'
  for item in alist:
    if current_max=='init' or item>current_max:
      current_second_max=current_max
      current_max=item
    elif current_second_max=='init' or item>current_second_max:
      current_second_max=item
  if current_second_max=='init':
    if current_max=='init':
      return False
    return current_max
  return current_second_max
      
blosum62_matrix={}
def load_blosum(from_file="/home/mmariotti/software/selenoprofiles/libraries/BLOSUM62sel"):
  try: assert blosum62_matrix
  except:
  #ncbi format
    ordered_aminoacids=[]
    main_diz={}
    cfile=open(from_file, 'r')
    cline=cfile.readline()
    index_line=0 #fake, check below
    while cline:
      if cline[0]=='#':
        cline=cfile.readline()
      elif not ordered_aminoacids:
       ordered_aminoacids=cline.split()
       index_line=-1
      else:
        for index_aa, num in enumerate(cline.split()[1:]):
          if not main_diz.has_key(ordered_aminoacids[index_line]):
            main_diz[ordered_aminoacids[index_line]]={}
          main_diz[ordered_aminoacids[index_line]][ordered_aminoacids[index_aa]]=int(num)

      cline=cfile.readline()
      index_line+=1
    cfile.close()
    main_diz['U']=main_diz['*']
    for k in main_diz:
      main_diz[k]['U']=main_diz[k]['*']

    blosum62_matrix=main_diz
  return blosum62_matrix

def blosum(aa1, aa2, matrix={}):
  if not matrix:
    matrix=load_blosum()
  if 'x' in [aa1, aa2]: #for my exonerate parser. actual Xs  (uppercase) are treated as in the blosum
    return 0
  #for my recoding of stop codons
  if aa1 in 'JZ':
    aa1='*'
  if aa2 in 'JZ':
    aa2='*'
  if matrix.has_key(aa1) and matrix[aa1].has_key(aa2):
    return matrix[aa1][aa2]
  else:
    raise Exception, "ERROR blosum score not defined for:"+aa1+' '+aa2

def similar_aas(aa1, aa2):
  '''returns True is the two aas are similar, false if not. They are defined as similar if they have a positive value in the blosum62 matrix
  '''
  similar_diz={"A":["S"], "R":["Q","K"], "N":["D","H","S","B"], "D":["N","E","B","Z"], "C":[], "Q":["R","E","K","Z"], "E":["D","Q","K","B","Z"], "G":[], "H":["N","Y"], "I":["L","M","V"], "L":["I","M","V"], "K":["R","Q","E","Z"], "M":["I","L","V"], "F":["W","Y"], "P":[], "S":["A","N","T"], "T":["S"], "W":["F","Y"], "Y":["H","F","W"], "V":["I","L","M"]}
  if similar_diz.has_key(aa1):
    if similar_diz.has_key(aa2):
      return   (aa2 in similar_diz[aa1])
  return False

def parsed_blast_to_summary(pline, program='tblastn', chars_per_line=60 ):
  return blasthit(pline).pretty_summary()

def all_chars_in(astring):
  """ This function returns the list of characters contained in the input string. The characters are in the order of first appearance"""
  outlist=[]; chars_hash={}
  for c in astring:
    if not chars_hash.has_key(c):
      outlist.append(c)
      chars_hash[c]=1
  return outlist

def find_all(substring, sstring):
  """ Find all arg1 occunreces in arg2 , and returns their indices. work with overlapping occurrences """
  l=len(substring); out=[]
  for pos in range(len(sstring)+1-l):
    if sstring[pos:pos+l]==substring: out.append(pos)
  return out

default_genetic_code=1
def set_genetic_code(code):      
  """ Set the default translation table to this code; Input is numerical, follows NCBI standards (e.g. 1 is standard, 6 is ciliate).
  This affects later calls of transl(seq)  """
  return set_MMlib_var('default_genetic_code', code)

def get_genetic_code():          
  """ Get the default translation table code; numerical, follows NCBI standards"""
  return get_MMlib_var('default_genetic_code')

def get_genetic_code_table(code=None):    
  """ Returns a dictionary codon->aminoacid for the given code (numerical, NCBI standard)
  If code is not provided, the default set in MMlib is used"""
  if code is None:   code=get_genetic_code()
  return genetic_codes[code]

def transl(cds_seq, include_selenocysteine=False, gaps_to=None, code=None):
  '''translate a nucleotide sequence in aminoacids.
  Use code=X to give a integer identifying the genetic code to be used, as NCBI codes (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
  Use include_selenocysteine=1 to use U for TGA or UGA
  Use gaps_to to force a certain char for gaps codons (---); normally translated as -
  '''
  if code is None: code=default_genetic_code
  out=''
  i=0
  codon_table=genetic_codes[code]
  while i < len(cds_seq):
    codon= replace_chars(upper(cds_seq[i:i+3]), 'U', 'T')
    if include_selenocysteine and codon=='TGA':      out+='U' 
    elif gaps_to and codon=='---':                   out+=gaps_to
    elif codon_table.has_key(codon):                 out+=codon_table[codon]
    else:                                            out+='X'
    i+=3
  return out

def retrotransl(aa_seq, gaps_to='', codon_hash={}):
  """translate an aminoacid sequence back to coding sequence. The first codon in alphabetical order is considered, unless a different codon_hash is provided. Argument gaps_to can be used to provide the character for gaps. Notice that a three character argument should be provided"""
  if not codon_hash: codon_hash=retrotrans
  out=''
  if len(gaps_to)==1: gaps_to*=3
  if gaps_to and len(gaps_to)!=3: raise Exception, "retrotransl ERROR gaps_to should be a string with length 3 or 1! gaps_to="+str(gaps_to)
  for aa in aa_seq:
    aa=upper(aa)
    if gaps_to and aa=='-':          out+=gaps_to
    elif codon_hash.has_key(aa):     
      codon=codon_hash[aa]
      if type(codon)==list: codon=codon[0] #taking first in alphabetical order
      out+=codon
    else: out+='NNN'
  return out

class e_v:
  """  Class for coping with evalues without going crazy because of out of memory
  """

  def __init__(self, anything):
    self.string=str(anything)
    self.value='!'
    try:
      self.value=float(self.string)
    except ValueError:
      "nothing"


  def exponent(self):
    if 'e' in self.string:
      return int(self.string.split('e')[1])
    else:
      if self.value =='!':
        raise ValueError
      if self.value==0:
        return -1000
      elif self.value<1:

        a= str(self.value).split('.')[1]
        exp=-1
        while a and a[0]=='0':
          exp-=1
          a=a[1:]
        return exp
      elif self.value>=1:
        return len(str(self.value).split('.')[0])-1

  def is_minor_than(self, other_e_v, or_equal=0):
    if type(other_e_v).__name__ in ["str",'int','float']:     #converting if necessary
      other_e_v=e_v(str(other_e_v))
    if self.exponent() < other_e_v.exponent():      return True
    elif self.exponent() > other_e_v.exponent():    return False
    else:
      if self.string == other_e_v.string:        return or_equal
      else:
        if self.value== "!" or other_e_v.value== "!":
          raise ValueError, "ERROR in evalue class! can't compare the two evalues: "+self.string+' AND '+other_e_v.string
        else:
          if or_equal:        return self.value <= other_e_v.value
          else:              return self.value < other_e_v.value

  def __lt__(self, other):    return self.is_minor_than(other)
  def __le__(self, other):    return self.is_minor_than(other, or_equal=1)
  def __gt__(self, other):    return not self.is_minor_than(other, or_equal=1)
  def __ge__(self, other):    return not self.is_minor_than(other)

  def __repr__(self):    return self.string
  def __str__(self):     return self.string
  
  def __float__(self):   return self.value
  def __int__(self):     return int(self.value)
    
  def __eq__(self, other):    return self.is_minor_than(other, or_equal=1) and not self.is_minor_than(other)

def shortcut_log(evalue):
  """ utility for evalues e_v ... defining the log(0) as 200, since I think blast is not computing evalues < than e-200"""
  try:
    if 'e' in str(evalue):
      return -int(str(evalue).split('e')[1])
    elif evalue>=1:
      return 0
    elif evalue==0.0:
      return 200
    else:
      a= -int(round(math_log(evalue, 10)))
      if a>0:
        return a
      else: 
        return 0
  except:
    print "____________ERROR: ",[evalue], type(evalue)


def match_any_word(main_string, word_list, is_pattern=1, ignore_case=1):
  """ Given a string and a list of strings/perl_patterns, it returns 1 is any of them matches the string, 0 otherwise
  """
  for w in word_list:
    try: 
      if is_pattern:
        if ignore_case:
          pattern=re.compile(w, re.IGNORECASE)
        else:
          pattern=re.compile(w)
        if pattern.search(main_string):
          return 1          
      else:
         if ignore_case:
           if w.lower() in main_string.lower():
             return 1 
         else:
           if w in main_string:
             return 1          
    except:
      printerr('ERROR with pattern: '+w)
      raise
  return 0
    
def score_tags(blast_evalues_and_titles, positive=[], negative=[], neutral=[], verbose=0, max_n_titles=0):
  """Input: list of objects [evalue, title] coming from parsing a blast file; evalue should be a e_v class object. titles should be complete (not just gis). Other arguments are: positive, negative, neutral tags, to be searched in the titles. The lines matching a neutral tag are skipped. Lines matching a negative tag (or no positive tag) are score negatively. Lines matching a positive tag are scored positively. The score depends on the negative logarithm of the evalue.
  """

  score=0
  if max_n_titles==0:
    max_n_titles=len(blast_evalues_and_titles)
  for evalue, title in blast_evalues_and_titles[:max_n_titles]:
    score_of_this_tag=shortcut_log(evalue)

    if not  match_any_word(title, neutral):
      if match_any_word(title, negative):
        score-=score_of_this_tag
        if verbose:
          print 'NEGATIVE_M ',title, evalue, '*-'+str(score_of_this_tag)
      elif match_any_word(title, positive):
        score+=score_of_this_tag
        if verbose:
          print 'POSITIVE ',title, evalue, '*+'+str(score_of_this_tag)
      else:
        if verbose:
          print 'NEGATIVE ',title, evalue, '*-'+str(score_of_this_tag)
        score-=score_of_this_tag
    elif verbose:
      print 'NEUTRAL', title, evalue, '/'+str(score_of_this_tag)

  return score

  
def dict_to_py_obj(sstring):
  """Take as input a string, which is a dictionary (hash) as printed usually by python, like :    {'a':'asdasfas', 'b': 'asfdasf'} and return a copy of the original object that was run.
  """
  hash_out={}
  if sstring[0]=='{':
    sstring=sstring[1:]
  while sstring and not sstring[0]=='}':

    #key  
    if sstring[0] in ["'", '"']:
      #key is string
      key=sstring[1:].split(sstring[0]+  ":")[0]
      sstring=sstring[1+len(key)+2+1:]
    else:
      value_string=sstring.split(":")[0]
      key=option_value(value_string)
      sstring=sstring[len(value_string)+1:]

    #value
    if sstring[0] in ["'", '"']:
      #value is a string
      if sstring[1:].split(sstring[0]+",")[0]==sstring[1:]:
        #last value
        value=sstring[1:].split(sstring[0]+"}")[0]
        sstring=''

      else:
        value=sstring[1:].split(sstring[0]+",")[0]
        sstring=sstring[1+len(value)+2+1:]
    else:
      value_string =sstring.split(',') [0] 
      if value_string ==sstring.split(',') [0] ==sstring:
        #last value
        value=option_value(sstring[1:].split("}")[0])
        sstring=''
      else:
        value= option_value(value_string )
        sstring=sstring[len(value_string)+2:]
    
    hash_out[key]=value
  return hash_out
def dict_to_config_file(diz, fileout=''):
  """Take as input a dictionary (typically a opt dictionary) and return a string corresponding to a configuration file that would be loaded in memory as the input dictionary.
  """
  out=''
  for k in diz:
    out+='\n'+str(k)+' = '+str(diz[k])
  if out:
    out=out[1:]
  if fileout:
    write_to_file(out, fileout)
  return out

def float_generalized(stringg):
  try:
    a=float(stringg)
    return a
  except ValueError:
    if stringg[0]=='e':
      return float('1'+stringg)
def replace_chars(astring, chars_list, replace_to_this=''):
  return ''.join([c if not c in chars_list else  replace_to_this      for c in astring])
  # out=''
  # for c in astring:
  #   if not c in chars_list:
  #     out+=c
  #   else:
  #     out+=replace_to_this
  # return out

def debug(msg):
#  if opt.has_key('v') and opt['v']:
    'nothing'
#    print msg

opt=options()
colored_keywords={}
def printerr(msg, put_newline=0, how='', keywords={}, is_service=False):
  global printed_rchar
  if not keywords and colored_keywords: keywords=colored_keywords
  msg=str(msg)
  if put_newline:    msg=msg+'\n'
  no_color_msg=msg
  if printed_rchar:
    sys.stderr.write('\r'+printed_rchar*' '+'\r' )
    printed_rchar=0
  if sys.stdout.isatty() and not opt['no_colors']:
    if how:
      for c in how.split(','): 
        if not terminal_codes.has_key(c): raise Exception, "ERROR option 'how' for write was not recognized: "+str(c)+' ; possible values are: '+join([i for i in terminal_codes.keys() if i], ',')
        msg=terminal_codes[c]+msg+terminal_codes['']
    for word in keywords:
      code=''
      for c in keywords[word].split(','): code+=terminal_codes[c]
      msg= replace(msg, word, code+word+terminal_codes[''])
  sys.stderr.write(str(msg))
  if not is_service and 'log_file' in globals(): print >> log_file, str(no_color_msg),
  
def service(msg):
  """ see write function"""
  msg=str(msg)
  global printed_rchar, opt
  if sys.stderr.isatty() and  not opt['Q']:
    if printed_rchar:
      printerr('\r'+printed_rchar*' ', is_service=True )
    printerr( "\r"+msg, is_service=True)
    printed_rchar=len(msg)
  #if 'log_file' in globals(): print >> log_file, str(msg+'\n') #putting a newline

def verbose(msg, put_newline=0):
  global opt
  if put_newline:    msg=str(msg)+'\n'  
  if opt['v']:
    write( msg )
    if 'log_file' in globals(): print >> log_file, str(msg),
    
terminal_codes={'':'\033[0m', 'red':'\033[31m', 'green':'\033[32m', 'black':'\033[30m', 'yellow':'\033[33m', 'blue':'\033[34m', 'magenta':'\033[35m', 'cyan':'\033[36m', 'white':'\033[37m', 'bright':'\033[1m', 'dim':'\033[2m', 'underscore':'\033[4m', 'blink':'\033[5m', 'reverse':'\033[7m', 'hidden':'\033[8m'}
def write(msg, put_newline=0, how='', keywords={}):
  """ Function to extend the functionalities of the standard 'print'. First argument (put_newline) when set to 1 put a newline after the string passed, as print would normally do. The argument "how" can be given a color to write the message in that color (only for atty terminals). This is prevented if opt['no_colors'] is active.  The function write is coupled with function "service" which prints service message which are deleted when another service message is printed, or another message is printed with the write function. If you use service, you should only print things with "write".
Argument keywords allows to use certain colors (or other "how" arguments) for certain keywords. The argument is a hash of keywords and correspoding how arguments. for example if you want to higlight all "ERROR" in red, pass keywords={'ERROR':'red'} 
"""
  msg=str(msg)
  global printed_rchar, opt
  if not keywords and colored_keywords: keywords=colored_keywords
  if put_newline:     msg=msg+'\n'
  no_color_msg=msg
  if sys.stdout.isatty() and not opt['no_colors']:
    if how:
      for c in how.split(','): 
        if not terminal_codes.has_key(c): raise Exception, "ERROR option 'how' for write was not recognized: "+str(c)+' ; possible values are: '+join([i for i in terminal_codes.keys() if i], ',')
        msg=terminal_codes[c]+msg+terminal_codes['']
    for word in keywords:
      code=''
      for c in keywords[word].split(','): code+=terminal_codes[c]
      msg= replace(msg, word, code+word+terminal_codes[''])
  if printed_rchar:
    sys.stderr.write('\r'+printed_rchar*' '+'\r' )
    printed_rchar=0
  sys.stdout.write(msg)
  if 'log_file' in globals(): print >> log_file, no_color_msg, 
warning=write

def is_empty_file(filename):
  return is_file(filename) and os.stat(filename)[6]==0

def is_valid_blast_output(filename):
  b=bash('tail -30 '+filename)
  if not b[0] and 'Lambda' in b[1]:
    return True
  return False

def framesBlastToTranseq( Blastframe, seq_length  ):
### translate the blast notation of frames to the ebi transeq notation
  if Blastframe>0:
    return Blastframe
  else:
    return -( 1+ ( seq_length+ Blastframe+1  )%3  )

def overlapping( range1, range2):    
  if max(range1)<=min(range2)   or  max(range2) <= min(range1)   : 
    return False
  else:  
    return True


def positions_to_frame(pos_start, strand='+', chr_length=None):
  """ Tells the frame (in blast notation) of a genomic range. we just need the starting position (which is > the end if strand is negative).  The chr_length is needed only if the strand is negative.
  Blast frame notation: 
  
  +  positions:        123456789    123456789   123456789
  if starts here:      |__|__        |__|__       |__|__
  then frame:           +1             +2           +3

  -  positions:        123456789    123456789   123456789
  if starts here:         __|__|      __|__|     __|__|
  then frame:             -1            -2           -3 

  so the trick for neg strand  is:   pstart%3 - chrlength%3 -1    --> frame   ; it doesn't work only if pstart%3 - chrlength%3 -1 > 0, which is only whe pstart%3 is 2  and chrlength%3 is 0 --> supposed result 1, must be -2.  I correct it manually. more elegant solutions would be more expensive anyway, so let's do it this way.  """

  if strand=='+':
    frame=pos_start%3
    if frame==0: frame=3
  elif strand=='-':
    if chr_length is None: raise Exception, "ERROR positions_to_frame: if strand is negative, the chromosome length must be provided to know the frame! "
    frame =  pos_start%3 - chr_length%3 -1 
    if frame == 1 : frame=-2 
  else:   
    raise Exception, "ERROR positions_to_frame: strand not recognized: '"+strand+"'"    
  return frame

def are_overlapping_ranges( range1, range2):    #[min, max], [min, max]  ; return False or the new region containing both; 
  #look into picture (lab notes) for case listing. #n refers to the picture
  #taking mininum and maximum of the 2 ranges###################################################  ehiiii
  if range1[0]<range1[1]:
    r1_0, r1_1 =range1[0], range1[1]
  else:
    r1_1, r1_0 =range1[0], range1[1]
  if range2[0]<range2[1]:
    r2_0, r2_1 =range2[0], range2[1]
  else:
    r2_1, r2_0 =range2[0], range2[1]
  ################################################

  if r1_1<r2_0 or r2_1 < r1_0:
    return False      # 2 or 3

  elif r2_0<r1_0:
    if r2_1 < r1_1:
      return   [ r2_0, r1_1 ]      # 4
    else:
      return    [ r2_0, r2_1 ]      # 6
  elif r2_1  < r1_1 :
    return    [ r1_0, r1_1 ]    #5
  else:
    return    [ r1_0, r2_1 ]    #1




to_comment=['^', '$', '.', '[', ']', '|', '(', ')', '*', '+', '?', "\\" ]
def comment_for_gawk(string):  # comment characthers so that the string can be found with gawk
  out=''
  for char in string:
    if char in to_comment:
      out+='\\'
    out+=char
  return out
  
def write_to_file(string, filename):
  filett =open(filename, 'w')
  print >> filett, string
  filett.close()


def del_white(string):
   """str-> str with only one white char between words
   """
   string=' '+string+'  '
   c=1
   while len(string)!=c:
       while string[c-1]==' ' and string[c]==' ' and len(string)!=c+1:
           string=string[:c-1]+string[c:]
       c=c+1
   return string[1:-2]


def configuration_file(filename, just_these_keys={}):
  """ Utility to read configuration files and return an opt with the parsed information, having as keys the option names, and reporting the values. These are automatically converted to the "minimal" type. If the value is an integere, it is cast to integer, then float is tried, otherwise string. It the argument just_these_keys is provided, only a subset of keys are reported, those present as keys of the hash just_these_keys.
 
  Example of configuration file format:

  temp= /tmp/
  profiles_folder = /users/rg/mmariotti/profiles
  keep_blast=1

  # commented lines like this one are not read. empty lines are also not read.
  blast_options.DEFAULT       =    -b 5000 -F F 
  blast_options.DEEP       =       -b 10000 -F F 


  # For dotted keys (what is before the "="), hashes are returned. In the two last lines, a hash is created as value of the output opt, corresponding to the key "blast_options". This nested hash will have two keys: DEFAULT and DEEP, and the corresponding values will be the strings reported in the config file.

  exonerate_options.example.DEEP       =       prova
  
  # when multiple dots are present in the key, a more complex nested structure of hashes is created. If the configuration file had only the line above, the reported opt would be: {'exonerate_options':{'example':{'DEEP':'prova'}} }
 
  
  """
  opt={}
  for line in open(filename, 'r'):
    try:
      if line.split() and line[0]!='#':
        corrected_line=line
        if corrected_line[-1]=='\n': corrected_line=corrected_line[:-1]
        value = join(corrected_line.split('=')[1:], '=')
        key=del_white(line.split('=')[0])
        
        if '.' in key and not just_these_keys or just_these_keys.has_key(key.split('.')[0]):
          main_key=key.split('.')[0]
          secondary_keys=join(key.split('.')[1:], '.')
          if not opt.has_key(main_key):           opt[main_key]={}
          current_targeted_hash=opt[main_key]
          while '.' in secondary_keys: # key can be e.g.    set.firstfield.secondfield  = 3  in this case a key 'set' is created in opt, its value is a new hash. Then a key in this has is added (firstfield) its value being an empty hash, then a key (secondfield) is added to this hash and the value 3 is added to it.
            secondary_key=secondary_keys.split('.')[0]; secondary_keys=join(secondary_keys.split('.')[1:], '.')
            if not current_targeted_hash.has_key(secondary_key): current_targeted_hash[secondary_key]={}
            current_targeted_hash=current_targeted_hash[secondary_key]
          secondary_key=secondary_keys #now no . is left
          if not current_targeted_hash.has_key(secondary_key): current_targeted_hash[secondary_key]={}
          current_targeted_hash[secondary_key]=option_value( del_white(value) )
        
        if not just_these_keys or just_these_keys.has_key(key):
          opt[key] = option_value( del_white(value) )
    except Exception, e:
      print "ERROR reading configuration file " +filename+ " reading line: "+line[:-1]
      raise
  return opt
  
def no_gap(seq):
    """return the seq without gaps
    """
    out=''
    for char in seq:
        if char not in ['-', '.']:
            out+=char
    return out
def nogap(seq):
  return no_gap(seq)

def alignment_relative_pos(group_align, global_align, neutral_chars='UX*'):
    """given two alignment class objects or two dictionaries {prot:sequence} like an alignment.diz, this function maps each position of the first alignment into a position of the second one.
        it returns a list long as the sequence(s) in first alignment, where the i-th element is the position of i-th aminoacid in the second alignment
        NB each protein in group_align should be present in global_align;
        NB2: for each position, there should be at least one protein which does not contain '-' in that position.
    """
    def find_prot_no_gap(prot_diz, pos): #returns the name of a protein in the alignment which has not a gap in position pos
        for prot in prot_diz:
            if prot_diz[prot][pos]!='-':
                return prot

    if type(group_align)==dict:
        diz=group_align
        lengh=len(diz.values()[0])
    else:
        diz=group_align.diz
        lengh=group_align.length()

    
    if type(global_align)==dict:
        glodiz=global_align
    else:
        glodiz=global_align.diz

    output=range(lengh)
    for i in output:
        output[i]=-1
    
    prot=diz.keys()[0]    
    for pos in range(lengh):
        if diz[prot][pos]=='-':
            prot=find_prot_no_gap(diz, pos)

        aa=diz[prot][pos]
        pos_global=output[pos-1]+1
        while output[pos]==-1:
            aa_global=glodiz[prot][pos_global]
            if aa_global!='-':
                if lower(aa_global)==lower(aa) or (aa_global in neutral_chars) or (aa in neutral_chars):
                    output[pos]=pos_global
#        print aa
                else:
                    print 'alignment_relative_pos ERROR AminoAcids dont correspond:\n>GROUP: '+prot+'\n'+no_gap(diz[prot])+'\n>GLOBAL: '+prot+'\n'+no_gap(glodiz[prot])
    

                    return False
            else:
                pos_global+=1
    return output

def mapping_alignments(group_align, global_align, neutral_chars='UX*'):
    """given two alignment class objects or two dictionaries {prot:sequence} like an alignment.diz, this function maps each position of the first alignment into a position of the second one.
        it returns a list long as the sequence(s) in first alignment, where the i-th element is the position of i-th aminoacid in the second alignment
        NB2: for each position, there should be at least one protein which does not contain '-' in that position.
    """
    def find_prot_no_gap(prot_diz, pos, allowed_names={}): #returns the name of a protein in the alignment which has not a gap in position pos
      for prot in prot_diz:
        if (not allowed_names or allowed_names.has_key(prot)) and prot_diz[prot][pos]!='-':
          return prot

    if type(group_align)==dict:
        diz=group_align
        lengh=len(diz.values()[0])
    else:
        diz=group_align.diz
        lengh=group_align.length()

    
    if type(global_align)==dict:
        glodiz=global_align
    else:
        glodiz=global_align.diz

    output=range(lengh)
    for i in output:
        output[i]=-1
    
    prot=diz.keys()[0]    
    for pos in range(lengh):
      impossibile_position=False
  #printerr(str(pos))
      if diz[prot][pos]=='-' or not glodiz.has_key(prot):
        prot=find_prot_no_gap(diz, pos, glodiz)
        if not prot:
      #in this case, this is a position that exist only in group_align, since all common proteins have a gap in group_align in this position. So we put the same number computed in the previuos position, or the automata will stop, and we correct just before outputing putting -1 in these positions
      #print "alignment_relative_pos ERROR no common proteins between the alignments"
    
          output[pos]=output[pos-1]
          impossibile_position=True
          prot=diz.keys()[0]    
      if not impossibile_position:

        aa=diz[prot][pos]
        pos_global=output[pos-1]+1
        while output[pos]==-1:
            aa_global=glodiz[prot][pos_global]
            if aa_global!='-':
                if lower(aa_global)==lower(aa) or (aa_global in neutral_chars) or (aa in neutral_chars):
                    output[pos]=pos_global
      #        print aa
                else:
                    print 'alignment_relative_pos ERROR AminoAcids dont correspond:\n>GROUP: '+prot+'\n'+no_gap(diz[prot])+'\n>GLOBAL: '+prot+'\n'+no_gap(glodiz[prot])
  

                    return False
            else:
                pos_global+=1
                    
    #correcting impossibile positions
    for p in range(1, len(output)):
        if output[p]==output[p-1]:
          output[p]=-1

    return output
        

def transfer_alignment(group_align_or_diz, global_align_or_diz, neutral_chars='UX*', dont_shrink=False):
    """given two alignment class objects or two dictionaries {prot:sequence} like an alignment.diz, this function return a new alignment including all sequences of both alignments.
       you need to provide more common sequences as possible.
    
    """

    def find_prot_no_gap(names_diz, pos, prot_diz): #returns the name of a protein in the alignment which has not a gap in position pos
        for prot in names_diz:
            if prot_diz[prot][pos]!='-':
                return prot
        return False

##############
#####    to deal with both alignment() objects or simple dictionaries {protname:seq}
####
    diz={}
    if type(group_align_or_diz)==dict:
      for k in group_align_or_diz:
        diz[k]=group_align_or_diz[k]
      lengh=len(diz.values()[0])
    else:
    #alignment class
      for k in group_align_or_diz.diz:
        diz[k]=group_align_or_diz.diz[k]
      lengh=group_align_or_diz.length()

    glodiz={}
    if type(global_align_or_diz)==dict:
     for k in global_align_or_diz:
       glodiz[k]=global_align_or_diz[k]
    else:
      for k in global_align_or_diz.diz:
        glodiz[k]=global_align_or_diz.diz[k]

#############
    diz_common={}  #contains the names of the common prot between group_diz and global
    for prot in diz:
      c=complete_word(prot, glodiz) #word, or -1, or False
      if c and c!=-1:
        diz_common[c]= True
        if prot!=c:
          diz[c] = diz[prot]
          del diz[prot]
 
    if len(diz_common)==0:
      outt='transfer_aligment ERROR: no common proteins between the two alignments. check the headers!\n>global:\n'
      for k in glodiz:
        outt+=k+'\t'
      outt+="\n>group:\n"
      for k in diz:
        outt+=k+'\t'
      raise Exception( outt )
     
#    print diz_common
    
    
    output=range(lengh)
    for i in output:
        output[i]=-1
    
    prot=diz_common.keys()[0]    
    pos=0
    while pos <lengh:
        #print 'pos: '+str(pos)
        if diz[prot][pos]=='-':            prot=find_prot_no_gap(diz_common, pos, diz)
    #    print prot
        if not prot:    #none of the common proteins has something different from a gap in position: pos.
    #        output.append(-1)
        #    lengh+=1
            pos_global=output[pos-1]+1   #becomes 0 if pos is 0, for construction
            output[pos] = pos_global
            for p_name in glodiz:
                glodiz[p_name]= glodiz[p_name][:pos_global]+'-'+glodiz[p_name][pos_global:]
        
            prot=diz_common.keys()[0]
    #increment iwht a gap each global seq.     set output[pos] = output[pos-1]+1
    
    
        else:    
            aa=diz[prot][pos]
            pos_global=output[pos-1]+1     #becomes 0 if pos is 0, for construction
            while output[pos]==-1:
                aa_global=glodiz[prot][pos_global]
                if aa_global!='-':
                  if lower(aa_global)==lower(aa) or (aa_global in neutral_chars) or (aa in neutral_chars):
                          output[pos]=pos_global
                  else:
                    print 'transfer_alignment ERROR AminoAcids dont correspond in pos '+str(pos)+'('+aa_global+' != '+aa+') for prot:'+prot
                    print '>cluster'
                    print diz[prot]
                    print '>global'
                    print glodiz[prot]
                    return
                else:
                            pos_global+=1
        pos+=1
    
    #over:
    

    for p_name in diz:                                    # ->up 
        seq=''    
        for pos in range(len(output)):
          if pos==0:
            for i in range(output[0]):        #adding initial gaps
                seq+='-'

          else:
            for i in range(output[pos]-output[pos-1]-1):
                seq+='-'
          seq+=diz[p_name][pos]
        for i in range(len(glodiz[glodiz.keys()[0]])-1 -output[-1]):        #adding final gaps
          seq+='-'
        
        if diz_common.has_key(p_name):    #just checking. if everything is ok, I should bring this IF ->up; if common_diz.has_key(p_name): don't do
          'nothing'        
#        if seq!=glodiz[p_name]:
#            print "transfer_alignment ERROR: with "+p_name+'; seq_group = '+seq +' seq_glodiz = '+glodiz[p_name]
            
        else:
          glodiz[p_name]=seq

    #correcting desert columns defect:

    if not dont_shrink:
      glodiz.shrink()

    return glodiz



class simmetrical_hash(dict):
	""" ........... old and bad implementation; see symmetrical_dict instead """
	def get(self, k1, k2):
		try:		  	return self[k1][k2]
		except:			return self[k2][k1]
symmetrical_hash=simmetrical_hash		


class symmetrical_dict(dict):
  """ Symmetrical dictionary. Usage:

  h=symmetrical_dict()
  h['a']['b']= 'something'     # this and     h['b']['a']='something'    have the same effect
  
  print h['a']['b']  --> 'something'
  print h['b']['a']  --> 'something'
  print h['a']['x']  --> None          # NOTE THIS!
  
  h.has_keys('a', 'c') -> False     
  """
  
  #def get_value(self, a, b):           a, b=sorted([a, b]);    return self[a][b] if a in self and b in self[a] else None
  # def set_value(self, a, b, value):    
  #     a, b=sorted([a, b]);    
  #     if not a in self:        self[a]=self.subdict(parent=self, mainkey=a)
  #     self[a][b]=value
  
  def __getitem__(self, key):
    if not key in self:      self[key]=self.subdict(parent=self, mainkey=key)
    return dict.__getitem__(self, key)

  class subdict(dict):
    def __init__(self, parent, mainkey, *args, **kargs):
      dict.__init__(self)
      dict.__setitem__(self, '__parent__',  parent )
      dict.__setitem__(self, '__mainkey__', mainkey)
            
    def __getitem__(self, key):
      if key < dict.__getitem__( self, '__mainkey__' ): 
        return dict.__getitem__(      #parent dict--> get the right subdict  #                        # index it with the mainkey of this subdict #
                                              dict.__getitem__(   dict.__getitem__(self, '__parent__'), key),     dict.__getitem__(self, '__mainkey__')  )
      else:                         return dict.__getitem__(self, key)   if key in self else None
    def __setitem__(self, key, value):
      if key=='__parent__' or key=='__mainkey__':    
        dict.__setitem__(self, key, value)
      elif key < dict.__getitem__(self, '__mainkey__'): 
        dict.__getitem__(self, '__parent__')[key][dict.__getitem__(self, '__mainkey__')] = value
      else:          
        dict.__setitem__(self, key, value) 
	
  def has_keys(self, a, b):
    if b<a: a,b=b,a
    return  dict.has_key(self, a) and b in self[a]

  def all_keys(self):
    out=set(dict.keys(self))
    for x in dict.keys(self):    
      for x in dict.keys( dict.__getitem__(self, x) ):
        if x!='__mainkey__' and x!='__parent__': out.add(x)
    return list(out)

class AliError(Exception):
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class alignment:
  """objects:
        diz= {protname: seq_with_gaps}
        order=[protname1, protname2, ...protnameN]
 OLD: [  ss= {protname: ss_with_gap}
         consensus= sequence long as the alignment (type =string)  . see fill_consensus()
         consensus_ss= secondary_structure long as the alignment (type =string)  . see fill_consensus() 
       ]        
       methods:
        init:  alignment(alignment_file_path, format='fasta', short_titles=False)       --->loads each prot seq in self.diz  ={protname:seq}
        titles():   return a (ordered) list of the protnames of the alignment
        seq_of(protname):   return the sequence of an entry in the alignment. The fill_title function is eventually activated for this
        fill_title(title):  given a uncomplete title (as cut by some programs), it return the complete title instead, if a unique title was found.
        add(protname, seq): adds protname:seq to self.diz; it's considered the last in order
        remove(protname): remove chosen protein from alignment. after that, remove_useless_gaps() is run
        remove_useless_gaps(): if all sequences have a '-' in a certain position, it's removed from all sequences
        nseq(): returns number of seq in alignment
        length(): checks if all sequences have same length, and returns it.
        display(): print alignment in fasta simple format
        fasta(): returns a string with a fasta text of alignment, removing gaps in sequences
        check_conservation(pos, lett): returns the conservation percentage of LETTER in global POSition in the alignment.
        conservation(pos, permitted_lett=AA_LETT+'-', threshold=0.0): returns a dictionary like {letter: conservation} at POS, if conservation>threshold
        conservation_map(thresholdC=0.0): returns a list of dictionary like {letter: conservation}, one per position 
        check_conservation_ss(pos, lett): same as check_conservation but for secondary structure
        conservation_ss( pos, permitted_lett=SS_LETT+'-', threshold=0.0) : same as conservation but for secondary structure
        fill_ss(xff_folder, suffix=''): fill the ss object with secondary structures read from files in xff_folder es (using _h as suffix): ../xff_folder/vdac_drome_h.xff
        fill_consensus(): calculate a consensus for sequences taking the most present aa in each position, and put it in self.consensus. The same is done for ss (if present); consensus for ss is put in self.consensus_ss
        transfer_consensus(xff_folder='xff', suffix='', outfold='consensus'): it extracts secondary structure information from files (with chosen suffix) in xff_folder, calculates the consensus ss and trasfer it to each protein: a xff file per protein is opened and written in the output folder (outfold)
        order_by_similarity_with(title, matrix): it returns a list of protnames, which are ordered by similarity with the protein named "title", which must be present in the alignment. This protein is not returned in the list.
        all_positions_of(aa):  it returns a list of positions (starting from 0) corresponding to occurences of the aminoacid specified. This can occur in any of the seqs of the alignment
        identity_matrix():  it returns a list of n elements (where n is the length of the alignment) of 0 or 1-> if all sequences have the same aminoacid at this position
        sequence_identity():    it returns a float from 0 to 1, counting all perfect match positions / the lenght of the alignment . Makes sense for two seqs alignments
        order_with(o_function):   the o_function must be a cmp like function taking as input two [title, seq] objects. The self alignment is ordered accordingly
         NEED UPDATE!
  """
  def __init__(self, alignfile='', format='fasta', short_titles=False, check_uniq_titles=False):
    self.reset_derived_data()
    self.order=[]
    self.diz={}
    if type(alignfile)==dict:       #you can initialise the alignment object even using directly the diz object. In this case, random order is assigned
      for k in alignfile:
        self.diz[k]=alignfile[k]
        self.order.append(del_white(k))
    #      self.display()

    elif alignfile!='':
      self.load_file(alignfile, format=format, short_titles=short_titles, check_uniq_titles=check_uniq_titles)
          
    self.ss={}
    self.length()       # check if all sequences have same length
  def titles(self):
    return list(self.order)

  def load_file(self, alignfile, format='fasta', short_titles=False, check_uniq_titles=False):
    """ loading a fasta or clustalw format alignment"""
    if format=='fasta':
      if check_uniq_titles: all_short_titles={}
      for title,seq  in parse_fasta(alignfile):
        if short_titles:
          title=title.split()[0]
        if check_uniq_titles:
          if all_short_titles.has_key(title.split()[0]):
            raise Exception, "alignment-> load ERROR short title in the alignment is not uniq: "+title.split()[0]
          all_short_titles[title.split()[0]]=1

        if self.diz.has_key(title):
          printerr('alignment class WARNING: the title "'+title+'" is present in more than one sequence. only the last one will be stored in the alignment object.', 1 )
        else:
              self.order.append(title)
        self.diz[title]=seq
    elif format=='stockholm': self.load_stockholm(alignfile)
    else: #clustal
      for [title, seq] in getfastaaln(open(alignfile, 'r'), order=1):
        if short_titles:
          title=title.split()[0]
        if self.diz.has_key(title):
                printerr('alignment class WARNING: the title "'+title+'" is present in more than one sequence. only the last one will be caught in the alignment object.' )
        else:
                self.order.append(title)
    
        self.diz[title]=seq

  def load_stockholm(self, filename):
    """Loads a stockholm file into the alignment. ss is stored into attribute .ss if present, same thing for .rf """
    self.order=[]; self.diz={}; self.ss=''; self.rf=''
    fileh=open(filename, 'r')
    line=fileh.readline()
    while line and (line[0]=='#' or not line.split()):     line=fileh.readline()      
    #now in first alignment line
    while line and line[0]!='#' and line!='//\n' and line.split():      
      self.add(line.split()[0], line.split()[1]) #parsing first block
      line=fileh.readline()
    while line and (line[0]=='#' or not line.split()):
      if line.startswith('#=GC SS_cons'): self.ss=line.split()[-1]
      if line.startswith('#=GC RF'): self.rf=line.split()[-1]
      line=fileh.readline()
    #now in beginning of new block, or last line
    while line and  line!='//\n':
      if line.startswith('#=GC SS_cons'): self.ss+=line.split()[-1]
      elif line.startswith('#=GC RF'): self.rf+=line.split()[-1]
      elif line.split():        
        if self.has_title( line.split()[0]  ):
          self.set_sequence(  line.split()[0], self.seq_of(line.split()[0])+ line.split()[1]     )
        else: self.add( line.split()[0],  line.split()[1] )
      line=fileh.readline()
    fileh.close()
    if line!='//\n': raise Exception, "ERROR loading stockholm file "+filename+' : \\ was not found at the end of file'
    self.convert_sequences(replace, '.', '-')

  def load_cmalign_out(self, filename):
    """ load an alignment from a cmalign out. this is almost identical to a stockholm format, but with 1 line which disrupts it (I think it's a bug); it also set the cm_name attribute to the name of cm found in cmalignout file"""
    try: 
      self.load_stockholm(filename)
      return
    except: pass
    fileh=open(filename); tempfileh=open(temp_folder+'temp_cmalign_out.stk', 'w')
    line=fileh.readline()
    while line and not line.startswith('# cm name'):           line=fileh.readline()    
    line=fileh.readline(); line=fileh.readline()
    self.cm_name=line.split()[1]
    while line and not line.startswith('# seq idx'):           line=fileh.readline()
    line=fileh.readline(); line=fileh.readline(); line=fileh.readline() #skipping line which makes the stockholm parser crash
    while line:
      print >> tempfileh, line,
      line=fileh.readline()
    tempfileh.close()
    self.load_stockholm(temp_folder+'temp_cmalign_out.stk')

  def __repr__(self):
    return self.summary()

  def subset(self, titles):
    """ Return an alignment containing only certain sequences, defined by argument titles (hash or list, hash is faster), which must contain the complete titles or the first word. Order in outptu is preserved from self alignment """
    b=alignment()
    for title in self.titles():
      if title in titles or title.split()[0] in titles:
        b.add(title, self.seq_of(title))
    return b

  def fill_title(self, uncomplete_title, silent=False):
    """given a title, returns the title of the alignment from where this title where cut. If it is not unique, or none is matching, reFalse
    """
    if uncomplete_title in self.diz: return uncomplete_title
    title_that_match=complete_word(uncomplete_title, self.diz, word_match=True)
    #print title_that_match
    if not title_that_match:
      if len(uncomplete_title.split())>1:        return self.fill_title( uncomplete_title.split()[0], silent  )
      else:
        if not silent:         raise Exception,  "ERROR the aligment object have no sequence starting with the title: "+uncomplete_title
        return False    
    elif title_that_match==-1:
      if len(uncomplete_title.split())>1:        return self.fill_title( uncomplete_title.split()[0], silent  )
      else:
        if not silent:       raise Exception,    "ERROR the aligment object have no unique sequence starting with the title: "+uncomplete_title
        return False
    return title_that_match

  def seq_of(self, title):
    if self.diz.has_key(title):      return self.diz[title]
    elif type(title)==int:
      if title>=len(self.titles()):
        raise AliError("seq_of function ERROR the alignment object have no sequence with the id: "+str(title)     )
        return False
      return self.seq_of(self.titles()[title])
    else:
      #trying to recover the seq of a uncomplete title.  returning it only if it's unique
      complete_title= self.fill_title(title, True)
      if complete_title:        return self.diz[complete_title]
      #trying to recover from only the first word
      complete_title= self.fill_title(title.split()[0]+' ', True)
      if complete_title:
        return self.diz[complete_title]
      raise AliError("seq_of function ERROR the alignment object have no sequence with the title: "+title)
      return False

  def has_title(self, title, even_partial=False):
    if even_partial:      return bool(self.fill_title(title, silent=True))
    else:      return self.diz.has_key(title)

  def codeml_format(self):
    o='  '+str(len(self.titles()))+' '+str(self.length())+'\n'
    for title in self.titles(): 
      short_title=   title.split()[0]
      o+= '>'+short_title+'\n'+self.seq_of(title)    +'\n'
    return o[:-1]

  def phylip_format(self, chars_per_block=60):
    if not self.check_length(): raise Exception, "Sequences are not aligned! can't output phylip"
    short_titles= [ title.split()[0] for title in self.titles() ] 
    max_short_title_length= max( [ len(s) for s in  short_titles ] ) 
    o='  '+str(len(self.titles()))+' '+str(self.length())+'\n'
    for c in range(    (self.length()-1)/chars_per_block +1): #number of blocks
      for index, title in enumerate( self.titles() ):
        if c==0:  title_shown= short_titles[index].ljust(max_short_title_length)+' '
        else:     title_shown=''# ' '*max_short_title_length + ' '
        o+=title_shown +   self.seq_of(index)[  c*chars_per_block: (c+1)*chars_per_block     ]+'\n'
      o+='\n'
    return o
        
  def clustal_format(self, ali_chars=60):
    if not self.nseq(): return ''
    h={}
    for t in self.titles(): 
      if h.has_key( t.split()[0] ):  raise Exception, "ERROR clustalw_format: the first word of titles must be unique to call this function"
      h[t.split()[0]]=1
    max_length= max (  [len(t) for t in h]  )
    n_lines   = self.length() / ali_chars
    out='CLUSTAL W (x.xx) multiple sequence alignment\n'
    if self.length() % ali_chars: n_lines+=1
    for line_index in range(n_lines):
      for t in self.titles():        out+=  t.split()[0].ljust(max_length)+' '  +  self.seq_of(t) [ line_index*ali_chars : (line_index+1)*ali_chars ] + '\n'
      if line_index != n_lines-1:    out+='\n\n'            
    return out

  def fill_sides(self, source_ali, inplace=False, wild_chars=[]):
    """ This function is useful when you have an alignment coming for example from a genomic prediction program, which aligns only a portion of the query. 
    This function takes the full seq of the query from source_ali, finds the missing parts and add them.
    It also replaces * in the self query with U.
    If the sequences do not match, if specified, attempt to match them considering wild_chars. The sequence returned will have the sequence matching perfectly the one from source ali.
    NOTE! If the sequence in self contains characters considered special by the re module, the behavior of this function is unpredictable.
    """
    if self.nseq()!=2:   raise Exception, 'alignment->fill_sides ERROR this function can be applied only to alignments with 2 sequences'
    if inplace:  a=self;     self.reset_derived_data()
    else:        a=self.copy()
    common_titles={} #complete_title: title ; they can be different, in that case title is the one in the self while complete_title is the one in the source_ali
    for title in a.titles(): 
      complete_title= source_ali.fill_title(title, silent=True)   #is false if it doesn't have it
      if complete_title:      common_titles[complete_title]=title
    if len(common_titles)!=1:      raise Exception, 'alignment->fill_sides ERROR too few or too many common titles between the two alignment provided (must be 1, it is: '+str(len(common_titles))+')'
    
    title_source_ali=common_titles.keys()[0];   title_self=common_titles[title_source_ali]
    complete_query_seq=           nogap(source_ali.seq_of(title_source_ali))
    partial_aligned_query_seq=    a.seq_of(title_self)
    pos_start= complete_query_seq.find(nogap(partial_aligned_query_seq)) #0 BASED
    #wild chars match
    if pos_start==-1 and wild_chars:
      pattern=re.compile('('+replace_chars(nogap(partial_aligned_query_seq), wild_chars, '.')+')')       
      match=pattern.search(complete_query_seq)
      if match:     
        pos_start= match.start()
        s=match.groups()
        seq=s[0]
        #adding gaps to seq
        find_gap= partial_aligned_query_seq.find('-')
        while find_gap!=-1:
          seq=seq[:find_gap]+'-'+seq[find_gap:]
          find_gap= partial_aligned_query_seq.find('-', find_gap+1)
        partial_aligned_query_seq=seq
        a.set_sequence(title_self,  partial_aligned_query_seq)
    ### NOTE the following code block is almost surely useless, but I keep it just for security. I now resolve the wild char match using the re module above.
    if pos_start==-1 and wild_chars:   # here trying to match considering some possible wild chars. In these positions the aminoacid is inferred from the complete_sequence in source_ali
      wild_char='%'
      wild_char_pos= partial_aligned_query_seq.find(wild_char)  ; next_wild_char_pos=partial_aligned_query_seq.find(wild_char, wild_char_pos+1)
      while wild_char_pos!=-1:
        bit_before_wild_char= partial_aligned_query_seq[:wild_char_pos]
        bit_before_wild_char_pos=complete_query_seq.find(nogap(bit_before_wild_char))
        while bit_before_wild_char_pos!=-1:          
          if      next_wild_char_pos==-1:   seq_to_match_in_partial_query_seq =   nogap(  partial_aligned_query_seq [   wild_char_pos+1:]  )
          else:                             seq_to_match_in_partial_query_seq =   nogap(  partial_aligned_query_seq [   wild_char_pos+1: next_wild_char_pos ]  )
          seq1=   complete_query_seq[bit_before_wild_char_pos+1+len(nogap(bit_before_wild_char)):]
          if seq1.startswith(seq_to_match_in_partial_query_seq): # checking that the match found by the sequence just before the wild char is the correct one. to know this, I check the characters right after the wild char in the complete sequence: they must correspond with the sequence in the partial sequence after the wild char.
            aligned_to_wild_char= complete_query_seq[bit_before_wild_char_pos+len(nogap(bit_before_wild_char))]
            partial_aligned_query_seq= partial_aligned_query_seq[:wild_char_pos]+ aligned_to_wild_char +partial_aligned_query_seq[wild_char_pos+1:]            
          bit_before_wild_char_pos=complete_query_seq.find(nogap(bit_before_wild_char), bit_before_wild_char_pos+1)
        wild_char_pos= partial_aligned_query_seq.find(wild_char, wild_char_pos+1)
      a.set_sequence(title_self,  partial_aligned_query_seq)
      pos_start= complete_query_seq.find(nogap(partial_aligned_query_seq)) #0 BASED    
    if pos_start==-1: #rechecking if wild chars mode, or just checking to raise if no wild_chars
      raise Exception, 'alignment->fill_sides ERROR sequences don\'t match for title: '+title_self+' -> '+nogap(a.seq_of(title_self))+" can't be found in "+complete_query_seq
    pos_end  = pos_start+len(nogap(a.seq_of(title_self)))-1 
    a.diz[title_self]= complete_query_seq[:pos_start] + a.seq_of(title_self) + complete_query_seq[pos_end+1:]
    for t in a.titles():
      if t!=title_self:
        a.diz[t]='-'*pos_start+a.seq_of(t)+'-'*len(complete_query_seq[pos_end+1:])
    if not inplace: return a        

  def remove(self, protname, dontremoveuselessgaps=False):
    del self.diz[protname]    
    if type(self.ss) == dict and self.ss.has_key(protname):       del self.ss[protname]
    self.reset_derived_data()
    for i in range(len(self.order)):
      if self.order[i]==protname:
        self.order.pop(i)
        break
    if not dontremoveuselessgaps:      self.remove_useless_gaps()
    return i

  def change_title(self, title, new_title):
    """Change the title of a sequence in the alignment to new_title. the order of the sequence in the alignmetn is preserved """
    try:     index= self.titles().index(title)
    except:  raise Exception, "alignment->change_title() ERROR the title provided is not found among the titles of this alignment: "+title
    seq=   self.seq_of(title)
    self.remove(title, dontremoveuselessgaps=True)
    self.add(new_title, seq, index)

  def is_gap_column(self, pos, char='-'):
    """ Return true is all sequence have a gap (or another character in input) at position pos (0-based!)"""
    for title in self.titles():
      if self.seq_of(title)[pos]!=char:          return False
    return True

  def stockholm(self, block_length=150):      
    """ stockholm output for alignments. if a .ss attribute is defined, a #=GC line with SS_cons will also be included. Same thing for .rf attribute """      
    stockholm_out='# STOCKHOLM 1.0 \n\n'
    max_length_title=max([len(t.split()[0]) for t  in self.titles()]+[12] ) #12 is for length of #GC tag
    print_ss=False; print_rf=False                                                                      
    try:  assert self.ss; assert len(self.ss)==self.length(); print_ss=True
    except: pass
    try:  assert self.rf; assert len(self.rf)==self.length(); print_rf=True
    except: pass
    for i_block in range(  1+(self.length()-1)/block_length   ):
      for title in self.titles():            stockholm_out+=(title.split()[0]).ljust(max_length_title)+' '+replace(self.seq_of(title)[i_block*block_length:(i_block+1)*block_length], '-', '.')+'\n'
      if print_ss:   stockholm_out+='#=GC SS_cons'.ljust(max_length_title)+' '+self.ss[i_block*block_length:(i_block+1)*block_length]+'\n'
      if print_rf:   stockholm_out+='#=GC RF'.ljust(max_length_title)+' '+self.rf[i_block*block_length:(i_block+1)*block_length]+'\n'
      stockholm_out+='\n'
    stockholm_out+='//'
    return stockholm_out

  def remove_useless_gaps(self):
    """if all sequences have a '-' in a certain position, it's removed from all sequences
    """
    if not self.diz:      return
    columns_removed={}
    for  pos in range( self.length() ):
      if self.is_gap_column(pos):        columns_removed[pos]=1  #annotating: we have to remove this
      pos+=1
    old_length=self.length()
    for title in self.titles():
      old_seq=self.seq_of(title)
      new_seq=''
      for  pos in range(  old_length  ):
        if not pos in columns_removed:         new_seq += old_seq[pos]
      self.set_sequence(title, new_seq)
    if columns_removed: self.reset_derived_data()
    return sorted(columns_removed.keys())

  remove_empty_columns=remove_useless_gaps

  def add(self, protname, seq, index=None):
      if self.diz.has_key(del_white(protname)):
        if seq == self.seq_of(del_white(protname)):
          printerr('alignment class WARNING: the title "'+protname+'" is already present in one sequence. Only one will be kept.\n')
        else:
          printerr('alignment class WARNING: the title "'+protname+'" is already present in one sequence and the sequences DIFFER! only the last one added will be kept.\n')
        self.set_sequence(protname, seq)
      else:
        self.diz[del_white(protname)]=seq
        if index is None:
          index=len(self.order)
        self.order.insert(index, del_white(protname))
      self.reset_derived_data()

  def __nonzero__(self):
    return bool(self.order)

  def nseq(self):
    return len(self.diz)
  def length(self):
    if self:    return len(self.seq_of(0))
    else:       return 0
  def check_length(self):
    """ Parse the alignment and returns False if it is not true that all sequences have the same length, True if it is regular"""
    l=0
    for title in self.titles():
      if not l:       l=len(self.seq_of(title))
      elif           l!=len(self.seq_of(title)): return False
    return True

  def translate(self):
    """For coding sequence alignments. Returns a copy of the alignment with the same titles and translated sequences """
    a=alignment()
    for t in self.titles():
      a.add(t, transl(self.seq_of(t)))
    return a
  
  def fill_ss(self, xff_folder, suffix=''):
    """fill the ss object with secondary structures read from files in xff_folder
    es (using _h as suffix): ../xff_folder/vdac_drome_h.xff
    """
    for prot in self.diz:
        temp=getfastaxff(open(xff_folder+'/'+prot+suffix+'.xff', 'r').readlines(), 0)
        if len(temp)==1:
          sss=temp[temp.keys()[0]]                
        else:
          sss=temp[prot]
        ss_cont=0
        ss_to_add=''
    #            print prot                             
        for pos in range(self.length()):
          a=self.diz[prot][pos]
          if a!='-':              
              if sss[0][ss_cont]==a:
                ss_to_add+=sss[1][ss_cont]
              ss_cont+=1
          else:
            ss_to_add+='-'
        self.ss[prot]=ss_to_add                

  def display(self, to_file='', return_it=False, fasta=False):
    if to_file:
      newfile=open(to_file, 'w')
      for t in self.titles():     
        seq=self.seq_of(t)
        if fasta: seq=fasta(seq)
        print >> newfile, ">"+t+'\n'+seq
      newfile.close()
    else:
      out=''
      for t in self.titles():
        seq=self.seq_of(t)
        if fasta: seq=fasta(seq)
        out+= '\n>'+t+'\n'+seq
      if out:      out=out[1:]
      if return_it: return out 
      else:         print  out
      
  def summary(self, titles=[]):
    out=''
    for prot_index in range(len(self.order)):
        if titles:
          prot=titles[prot_index] 
        else:
          prot=self.order[prot_index]
        out+= '\n>'+prot+'\n'+self.diz[prot]
    if out:
      out=out[1:]
    return out
      
  def fasta(self, protein=''):
    out=''
    for prot in self.titles():
        if protein=='' or prot==protein :          out+='>'+prot+'\n'+no_gap(self.seq_of(prot))+'\n'
    return out

  def aligned_fasta(self):
    out=''
    for prot in self.titles(): out+='>'+prot+'\n'+self.seq_of(prot)+'\n'
    return out

  def check_conservation(self, pos, lett):
    """ alignment = list of sequences; returns the conservation percentage of LETTER in global POSition in the ALIGNMENT. From 0.0 to 1.0 DEPRECATED"""
    nmatch=0.0
    for prot in self.diz:
        if self.diz[prot][pos]==lett:
          nmatch+=1
    cons=float(nmatch/len(self.diz))
    if len(str(cons))>5:                    #approximation to 0.nnn         ????????
        cons=float(str(cons)[:5])
    return cons

  def conservation(self, pos, permitted_lett='', threshold=0.0, exclude={}):
    """ alignment = list of sequence; returns a dictionary like {letter: conservation} at POS, if conservation>=threshold; you can provide a hash of titles to be excluded. If permitted_lett=='', then all letters are permitted. Pos is 0 based!!!"""
    counts={}; titles_considered=0; out={}
    for title in self.titles():
      if not  exclude.has_key(title):
        aa=self.seq_of(title)[pos]
        counts[aa]=counts.setdefault(aa, 0)+1
        titles_considered+=1
    for k in counts:
      if not permitted_lett or k in permitted_lett:
        percent=counts[k]/float(titles_considered)
        if percent>= threshold:
          out[k]=percent
    return out

  def most_conserved(self, pos, permitted_lett=AA_LETT, threshold=0.0):
    """return a tuple (aa,conservation) . If only gaps are present at that position, ('', 0.0) is returned DEPRECATED """
    cons_profile=self.conservation(pos, permitted_lett, threshold)
    max_percent, max_aa =0.0, ''
    for aa in cons_profile:
      if cons_profile[aa] > max_percent:
        max_percent=cons_profile[aa]
        max_aa=aa
    return (max_aa, max_percent)

  def check_conservation_ss(self, pos, lett):
    """ alignment = list of sequences; returns the conservation percentage of LETTER in global POSition in the ALIGNMENT. From 0.0 to 1.0 DEPRECATED"""
    nmatch=0.0
    for prot in self.ss:
        if self.ss[prot][pos]==lett:
          nmatch+=1
    cons=float(nmatch/len(self.ss))
    if len(str(cons))>5:                    #approximation to 0.nnn         ????????
        cons=float(str(cons)[:5])
    return cons
  def conservation_ss(self, pos, permitted_lett=SS_LETT+'-', threshold=0.0):
    """ alignment = list of sequence; returns a dictionary like {letter: conservation} at POS, if conservation>threshold DEPRECATED"""
    out={}
    for lett in permitted_lett:
        n=self.check_conservation_ss(pos, lett)
        if n>threshold:
          out[lett]=n
    return out
  def fill_consensus(self):
    """calculate a consensus for sequences taking the most present aa in each position, and put it in self.consensus
       if secondary structure information is present, this is done even for ss; consensus for ss is put in self.consensus_ss
       DEPRECATED
    """
    def max_value(diz, exclude=[]): #returns the key with a higher value associated, excluding keys in exclude list
        max_v=None
        the_key=None
        for key in diz:
          if diz[key]>=max_v and not key in exclude:
              the_key=key
              max_v=diz[key]
    #            if the_key==None: return '-'
        return the_key
    self.consensus=''
    for pos in range(self.length()):
        self.consensus+=max_value(self.conservation(pos), ['-'])
    if self.ss!={}:
        self.consensus_ss=''
        for pos in range(self.length()):
          self.consensus_ss+=max_value(self.conservation_ss(pos), ['-'])

  def position_in_seq(self, title, position):
    """ this function computes the position relative to a single sequence in the alignment, given the position in the alignment. NB !! All positions are 1 based here.
        Also, if you ask for a position which is a gap in the desired sequence, the position returned refers actually to the last non gap position in the seq (0 if non existing)
    """
    seq=self.seq_of(title)[:position]
    n_gaps_until_position=0
    for s in seq:
      if s=='-':         n_gaps_until_position+=1
    return position-n_gaps_until_position

  def position_in_ali(self, title, position):
    """ this function computes the position in the alignment in which the sequence "title" has position "position". it is like the contrary of position_in_seq.
    NB positions are 1 based here!
    """
    seq=self.seq_of(title)
    pos_seq=0
    for p, aa in enumerate(seq):
      if aa!='-':        pos_seq+=1
      if pos_seq==position:        return p+1

  def conservation_map(self, thresholdC=0.0, exclude={}, dont_save=False):
    """returns a list of dictionary like {letter: conservation}, one per position. Titles in list exclude will not taken into account    """
    exclude_key=  join(sorted(exclude.keys()), '&') # to index self.conservation_map_data to keep the data
    if self.conservation_map_data is None:      self.conservation_map_data={}
    if  not self.conservation_map_data.has_key( exclude_key ):
      out=[]
      for pos in range(self.length()):         out.append(self.conservation(pos, threshold=thresholdC, exclude=exclude))
      if dont_save: return out
      self.conservation_map_data[exclude_key]=out
    return self.conservation_map_data[exclude_key]

  def conservation_quadratic_score(self, titles=None):
    """ Return a list of scores, one per position of the alignment, computed as the sum of proportion of any character found at that position (except for character -, which then are like scored negatively). Titles can be used to call the function on a subset of the alignment."""
    exclude = {}    if titles is None    else    dict.fromkeys(  [t for t in self.titles() if not t in titles]  )
    c_map=self.conservation_map( exclude=exclude )
    out=[] #one score per position
    for pos in range(self.length()):
      cons_dict_pos= c_map[pos]  #pos 0 based
      score=0.0
      for char in cons_dict_pos:    
        if char!='-': score+=  cons_dict_pos[char]**2
      out.append(score)
    return out

  def reset_derived_data(self):
    """Delete the attributes derived from computation of the alignment, as for example .conservation_map_data. This function must be run everytime some modification to the self alignment is done. """
    self.conservation_map_data=None

  def trim_columns(self, max_non_gaps=0.1, inplace=True, remove_empty_seqs=False):
    """Clean the alignment by removing the desert columns. Only the columns for which at least X sequences have something different than a gap are kept """
    #getting list of positions to remove
    desert_columns=  self.find_desert_columns(max_non_gaps=float(max_non_gaps), length_clusters=1)
    positions_to_remove={}
    for pos, lenght in desert_columns: 
      for i in range(lenght):         positions_to_remove[ pos+i ]=True      
    ## now building new ali skipping those positions
    a=alignment()
    for title in self.titles():
      seq=''
      for pos in range(self.length()):
        if not pos in positions_to_remove:
          seq+=   self.seq_of(title)[pos]
      a.add(title, seq)
    if not remove_empty_seqs: a.remove_empty_seqs()
    a.reset_derived_data()
    if not inplace:  return a
    else:          
      self.__dict__=deepcopy(a.__dict__); 
      return sorted(positions_to_remove.keys())

  def consensus_sequence(self, threshold=0.0, sec_char='', exclude={}):
    """This function computes a consensus sequences taking into account all sequences in the alignment (apart from the titles in the input hash exclude). Not all the columns are taken into account: only those having at maximum \
"threshold" gaps (in proportion). If sec_char is set to a non-False value, any column containing a "U" will return a sec_char, regardless of number of gaps in this column. 
If a master alignment is provided, the conservation threshold is checked with this alignment instead of with self. This provided that the master alignment must have column numbering identical to self.
"""
    seq=''
    conservation_map_data=self.conservation_map(exclude=exclude)
    for pos in range(self.length()):
     # looking at the conservation profile at this position. sorting the aminoacid (or gap character) according to their representation in this column
      sorted_keys=sorted(conservation_map_data[pos].keys(), key=conservation_map_data[pos].get, reverse=True)
      if sec_char and 'U' in sorted_keys:                                          seq+=sec_char
      elif len(sorted_keys)==1 or sorted_keys[0]!='-':                             seq+=sorted_keys[0]
      elif conservation_map_data[pos].setdefault('-', 0.0) <= threshold and len(sorted_keys)>1:    seq+=sorted_keys[1]       
      else:                                                                        seq+='-'
    return  seq

  def transfer_consensus(self, xff_folder='xff', suffix='', outfold='consensus'):
    """it extracts secondary structure information from files (with chosen suffix) in xff_folder, calculates the consensus ss and trasfer it to each protein:
    a xff file per protein is opened and written in the output folder (outfold) 

    """
    self.fill_ss(xff_folder, suffix)
    self.fill_consensus()
    for prot in self.diz:
    #    print prot
        seq=no_gap(self.diz[prot])
        ss_cons=''
        rel_pos=alignment_relative_pos({prot:seq}, self.diz)
        for p in range(len(rel_pos)):
          ss_cons+=self.consensus_ss[rel_pos[p]]
        newfile=open(outfold+'/'+prot+suffix+'C.xff', 'w')        
        print >>newfile, '>'+prot+'\n'+seq+'\n#'+ss_cons
        newfile.close()

  def order_by_similarity_with(self, title, matrix=''):
    """ returns an ordered list of protnames ordered by the similarity with the sequence named title, which must be in the alignment (this protname is not reported)
        a matrix can be specified for scoring, otherwise the identity matrix is used. matrix object in input must be an hash of hashes. e.g. identity matrix would be like {'aa1':{'aa1':1, 'aa2':0, 'aa3':0, ....}, 'aa2':{'aa1':0, 'aa2':1, 'aa3':0...}, ...}
    """
    work_list=[]
    score_diz={}
    for protname in self.titles():
      if protname!=title:
        work_list.append(protname)
        score_diz[protname]=0
        for pos in range(self.length()):
          #aa1=self.diz[title][pos]
          #aa2=self.diz[protname][pos]
          if not self.diz[title][pos]=='-':
            if matrix:
              score_diz[protname] += matrix[  self.diz[title][pos]   ][  self.diz[protname][pos]   ]    #YOU HAVE TO PROVIDE THE MATRIX AS AN HASH LIKE  {'aa1':{'aa1':1, 'aa2':0, 'aa3':0, ....}, 'aa2':{'aa1':0, 'aa2':1, 'aa3':0...}, ...}
            else:    
              #simulating identity matrix
              score_diz[protname] += int(self.diz[title][pos]==self.diz[protname][pos])    

    def compare_function(x, y):
      if score_diz[x]>score_diz[y]:
        return -1
      elif score_diz[x]==score_diz[y]:
        return 0
      else: 
        return 1
    #  print score_diz
    work_list.sort(compare_function)
    #  print work_list
    return work_list

  def sequence_identity_at_position(self, pos):
    """Return the simple sequence identity at position pos (0 based, ali based). The most conserved character at each position is considered, unless it is a gap. In that case, either the second char proportion is reported (if any is present), or 0.0 is returned."""
    cmap_this_pos=self.conservation_map()[pos]
    sorted_chars=sorted (    cmap_this_pos.keys(),   key=lambda x: cmap_this_pos[x]   ,  reverse=True   )
    the_char=sorted_chars[0]
    if the_char == '-':
      if len(sorted_chars)==1: return 0.0
      else:                    the_char=sorted_chars[1]
    return cmap_this_pos[the_char]    

  def sequence_identity_list(self):
    """Return a list of same length of alignment, with values computed at each position by function sequence_identity_at_position"""
    return [self.sequence_identity_at_position(pos)   for pos in range(self.length())]      

  def sequence_identity(self, count_terminal_gaps=True):
    """ returns the seq identity of the alignment. if count_terminal_gaps=False, then terminal gaps of the sequence having the longest one are excluded from the count of the total length of the alignment"""
    if len(self.titles())==1:
      return 1.0
    if len(self.titles())==0:
      return 0.0
    n_matches=0
    for pos in range(self.length()):
        aa=self.seq_of(self.titles()[0])[pos]
        index_title=1
        this_pos_is_mismatch=False
        while index_title < len( self.titles() ) and not this_pos_is_mismatch:
          if self.seq_of(self.titles()[index_title])[pos]!=aa:
            this_pos_is_mismatch=True
          index_title+=1
        if not this_pos_is_mismatch:
          n_matches+=1
    length=self.length()
    if not count_terminal_gaps:
      n_boundaries, c_boundaries =[], [] 
      for i in self.titles():
        a=self.boundaries_of(i)
        n_boundaries.append(a[0])
        c_boundaries.append(a[1])
      nterminal_gap_length= max (n_boundaries )
      cterminal_gap_length= self.length()-1 - min ( c_boundaries )
      length-= (nterminal_gap_length+cterminal_gap_length)
    return float(n_matches)/length      
  
  def boundaries_of(self, title):
    """returns the limits, starting from 0, of a seq in the alignment, meaning the first and last positions with something different than a gap """
    met_the_start=False
    max_gap_pos=-1
    for pos in range(self.length()):
      if self.seq_of(title)[pos]=='-':
        if not met_the_start:
          max_gap_pos=pos
      else:
        last_non_gap_pos=pos
        met_the_start=True
    return (max_gap_pos+1, last_non_gap_pos)

  boundaries=boundaries_of

  def all_positions_of(self, aa, except_in_seq=[], minimum_number=1):
    """it returns a list of positions (starting from 0) corresponding to occurences of the aminoacid specified. This can occur in any of the seqs of the alignment.
     """
    out_list=[]
    for pos in range(self.length()):
      found=0
      cont_prot=0
      while found<minimum_number and cont_prot< self.nseq():
        if self.seq_of(   self.order[cont_prot] )  [pos] == aa and not (self.order[cont_prot] in except_in_seq) :
          found+=1
        cont_prot+=1
      if found>=minimum_number:        out_list.append(pos)
    return out_list

  def identity_matrix(self):
    id_matrix=[]
    for i in range(self.length()):
      i_seq=0
      all_the_same=1
      aa=self.seq_of(self.titles()[0])[i]
      while i_seq < self.nseq() and all_the_same:
        if self.seq_of(self.titles()[i_seq])[i]!=aa:
          all_the_same=0
        i_seq+=1
      id_matrix.append(all_the_same)

    return id_matrix

  def local_mismatches(self, min_length=4, min_match_length=4, invert=0):
    """parse a 2 seq alignment and return the coordinates of local mismatches. These are defined as windows of minimal length min_length,  surrounded by strecthes of 100% conservation of minimum length: min_match_length.  indexes are starting with 0 !!!! DO NOT RELY ON THIS FUNCTION, REBUILD IT IF YOU HAVE TO USE IT.
    """
    if self.nseq()!=2:
      raise Exception, 'cannot use this function (alignment.local_mismatches) on more or less than 2 sequences'
    out=[]
    id_matrix=self.identity_matrix()
    id_matrix_string=''
    for i in id_matrix:
      id_matrix_string+=str(i)

    i_string=0
    start_mismatch=''
    gaps=0
    while i_string+ min_match_length< len(id_matrix_string):
      #print      id_matrix_string[i_string:i_string+ min_match_length+1]
      #print self.best_ali[0][i_string:i_string+ min_match_length+1]

      if self.seq_of(self.order[1])[i_string]=='-':
        gaps+=1
      if (id_matrix_string[i_string:i_string+ min_match_length+1]=='1'*min_match_length+'0' and not invert) or (invert  and  id_matrix_string[i_string:i_string+ min_match_length+1]=='0'*min_match_length+'1'):
        start_mismatch=i_string+ min_match_length #cannot be 0 by definition, so I can test if start_mismatch
#        print "start found"
      if (id_matrix_string[i_string:i_string+ min_match_length+1]=='0'+('1'*min_match_length)    and not invert) or ( invert and id_matrix_string[i_string:i_string+ min_match_length+1]=='1'+('0'*min_match_length)):
  #        print "stop found"
        stop_mismatch=i_string
        if start_mismatch and stop_mismatch-start_mismatch +1 >= min_length:          
          stop_mismatch-=gaps
          start_mismatch-=gaps

          out.append( (start_mismatch,stop_mismatch) )  ###if negative frame, you can have: [(564, 507)]

      i_string+=1
    return out

  def all_mismatches(self, count_gaps=True):
    """parse a 2 seq alignment and return the coordinates of local mismatches in a list  (first position is 0)
    """
    out=[]
    if self.nseq()!=2:
      raise Exception, 'cannot use this function (alignment.local_mismatches) on more or less than 2 sequences'
    for pos in range(self.length()):
      if (self.seq_of(self.order[0])[pos]!=self.seq_of(self.order[1])[pos] ) and (count_gaps  or (not count_gaps and self.seq_of(self.order[0])[pos]!='-' and self.seq_of(self.order[1])[pos]!='-')):
        out.append(pos)
    return out


  def local_matches(self, min_match=3 ):
    return self.local_mismatches(min_match, 1, 1)

  def matches(self, min_match=3):
    """same logic as the previous functions. This is RELIABLE. Output is basically a condensation of the identity matrix, additionally with the filtering out of short matches
    output: list of [start_match, end_match] , with indexes which are included in the match and 1-based.
"""
    out=[]
    id_matrix=self.identity_matrix()

    start_match=None
    now_in_match=False
    for i in range(len(id_matrix)):
      if id_matrix[i]==1:
        if not now_in_match:
          start_match=i
          now_in_match=True
      else:
        if now_in_match:
          now_in_match=False
          if (i-start_match)>=min_match:
            out.append((start_match, i-1))

    if now_in_match:
      now_in_match=False
      if (i-start_match)>=min_match:
        out.append((start_match, i-1))

    return out


  def sequence_identity_of(self, title_1, title_2, dont_count_gaps=False, dont_count_terminal_gaps=False):
    """ Returns the sequence identity of title_1 and title_2.  The columns which are gaps in both sequences are removed before computing. 
  if dont_count_terminal_gaps , the unaligned tails are not considered; if dont_count_gaps, all columns in which one of the sequences have a gap are not considered.
    """
    seq_1=   self.seq_of(title_1)
    seq_2=   self.seq_of(title_2)
    identical_pos=0;     gaps_pos=0 #those in title1 plus those in title2 not aligned to those in title1
    #removing positions with gaps in both seqs
    pos_to_remove=[]
    for pos in range( self.length() ):
      if seq_1[pos]=='-' and seq_2[pos]=='-': pos_to_remove.append(pos)
    for i in range( len(pos_to_remove)-1, -1, -1):
      seq_1= seq_1[:pos_to_remove[i]]+seq_1[pos_to_remove[i]+1:]
      seq_2= seq_2[:pos_to_remove[i]]+seq_2[pos_to_remove[i]+1:]
    # extreme cases: sequences are just gaps are empty strings
    if any( [all([c=='-' for c in seq_1]), all([c=='-' for c in seq_2])] )  or not seq_1 or not seq_2: return 0.0 
    #computing nterminal and cterminal gap tails
    if dont_count_terminal_gaps:
      c_terminal_gaps_seq1=0
      n_terminal_gaps_seq1=0
      n_terminal_gaps_seq2=0         
      c_terminal_gaps_seq2=0
      while seq_1[n_terminal_gaps_seq1]=='-':   n_terminal_gaps_seq1+=1
      while seq_1[-1-c_terminal_gaps_seq1]=='-':  c_terminal_gaps_seq1+=1
      while seq_2[n_terminal_gaps_seq2]=='-':   n_terminal_gaps_seq2+=1
      while seq_2[-1-c_terminal_gaps_seq2]=='-':  c_terminal_gaps_seq2+=1
      n_terminal_gaps=max(n_terminal_gaps_seq1, n_terminal_gaps_seq2) #only one of the two can be different than 0, for construction
      c_terminal_gaps=max(c_terminal_gaps_seq1, c_terminal_gaps_seq2) #only one of the two can be different than 0, for construction
    for pos in range( len(seq_1) ):
      if   seq_1[pos]=='-': 			    gaps_pos+=1
      elif seq_2[pos]=='-': 			    gaps_pos+=1
      elif seq_1[pos]==seq_2[pos]:		identical_pos+=1
    if not identical_pos: return 0.0           #to avoid dividing by zero when they share no positions
    if dont_count_gaps:														    return  float(identical_pos)/( self.length() -len(pos_to_remove) -   gaps_pos)
    elif dont_count_terminal_gaps:                    return  float(identical_pos)/( self.length() -len(pos_to_remove) - n_terminal_gaps - c_terminal_gaps)
    else :							                              return  float(identical_pos)/( self.length() -len(pos_to_remove) )



  def weighted_seq_identity_of(self, title, with_coverage=True, exclude={}):
    """ Scores a sequence against the rest of the alignment, returning a float from 0.0 to 1.0 -- but the limit values depend on the alignment, it almost never can reach 1. 
    It proceeds in the following way: it compares the input title seq with any other sequence. For each sequence, a weighted score is assigned in this way:
    for each position in which none of the two sequences have a gap, a weight is considered, as the percent of conservation of the aminoacid of the other (non-input title) sequence.
    The weighted score for this pair of sequences is the sum of identity flags times the weigth for each position (identity flag: 1 if the two aminoacid are identical, 0 otherwise), divided by the total sum of weights for this sequence.
    If the input title seq has no alignemnt overlap  with any other sequence in the alignemnt, it is not counted for the final score.
    The average and std deviation of this measure is an indication of the conservation of a profile
    """
    if self.nseq()==1: return 1.0
    
    exclude_hash={title:True}
    for k in exclude: exclude_hash[k]=True 
    conservation_m = self.conservation_map(  exclude=exclude_hash )
    coverage_per_position=[]  #meaning: the proportion of seqs with something not a gap / tot n seqs 
    for element in conservation_m:
      coverage=  sum([  element[k]   for k in element if k !='-'  ])
      coverage_per_position.append(coverage)
    #start and stop: zero based
    scores=[]
    for t in self.titles():
      if not t in exclude_hash:
        total_weight=0; score_this_seq=0
        pos_in_ali=0
        for pos in range(self.length()):
          #now pos_in_ali reflects position we're in, 1 based, in the original alignment with no seq added.        
          candidate_seq_here=self.seq_of(title)[pos]
          profile_seq_here=  self.seq_of(t)[pos]
          if  profile_seq_here =='-': continue 
          if not with_coverage and candidate_seq_here=='-': continue
          
          weight=       conservation_m [pos][profile_seq_here]
          if with_coverage: weight*= coverage_per_position[pos] 
          is_identity= candidate_seq_here == profile_seq_here
          score_this_seq+=  int( is_identity  )*weight
          total_weight+=weight
         
        if total_weight:   scores.append(       score_this_seq/float(total_weight)     )
        #print title.split()[0]+' -- '+t.split()[0]+' w: ', total_weight, 's:', score_this_seq, 'ws:', score_this_seq/float(total_weight)
        #raw_input('next?')

    if not len(scores): raise Exception, "alignment-> weighted_seq_identity_of ERROR the sequence "+title.split()[0]+' has not any aligned position with any of the other sequences!'
    return sum(scores)/float(len(scores))

  def sequence_identity_hash(self, dont_count_gaps=False, dont_count_terminal_gaps=False ):
    """Returns a simmetrical hash with the sequence identities of every sequence in the alignment against each other. the simmetrical hash can be queried using .get(title1, title2)
    Nterminal gaps are ignored if dont_count_terminal_gaps==True; 
    All gaps are ignored if dont_count_gaps==True    """
    h=simmetrical_hash()
    for index_1 in range(len(self.titles() )):
      title_1= self.titles()[index_1]
      h[title_1]={}
      for index_2 in range(index_1+1, len(self.titles())):
        title_2= self.titles()[index_2]        
        h[title_1][title_2]= self.sequence_identity_of(title_1, title_2,  dont_count_gaps=dont_count_gaps, dont_count_terminal_gaps=dont_count_terminal_gaps)

    return h				


  def average_sequence_identity(self, dont_count_gaps=False, dont_count_terminal_gaps=False):
    """ """
    if self.nseq()==1: return 1.0
    values=[]
    h=self.sequence_identity_hash(dont_count_gaps=dont_count_gaps, dont_count_terminal_gaps=dont_count_terminal_gaps)
    for k1 in h:
      for k2 in h[k1]:
        values.append(h.get(k1, k2))
    return sum(values)/float(len(values))


  def average_sequence_identity_of(self, title='', dont_count_gaps=False, dont_count_terminal_gaps=False ):
    """Returns the average sequence identity of the sequence with "title" with all the other sequences in the profile. positions in which "title" has a gap are considered always as mismatch, unless dont_count_gaps is set to True. in this case those positions are excluded from computation (also when dividing for total length).
    If dont_count_terminal_gaps is set to True, the unaligned tails are excluded from computation.
     """
    if len(  self.titles() )>1:
      seq=self.seq_of(title)
      seq_identity_per_title={}
      for tit in self.titles():
        if tit !=title:
          seq_identity_per_title[tit]=self.sequence_identity_of(title, tit, dont_count_gaps=dont_count_gaps, dont_count_terminal_gaps=dont_count_terminal_gaps )
      return sum( seq_identity_per_title.values() ) / len(seq_identity_per_title.keys())
    else: return 1.0

  def conservation_score(self, title='', matrix={}):
    """This function computes a conservation score of one sequence (title) with the rest of the alignment. At each position, the average blosum score between the aa at this position and the aas of each other sequence which doesn't bring a gap here is multiplied IMP times, where IMP is the relative importance of this position of the alignment and is defined by the proportion of profile sequences (all apart from title) which don't carry a gap here. When all other sequences carry a gap, 0 is added. When the "title" sequence carry a gap, -4 (gap_score) is added,  multiplied IMP times. Lastly, the score is divided by the lenght of the alignment.
    The variable matrix provided should be a hash-like scoring matrix, a hash of hashes with as final values integers or floats for the evaluation of pairwise amino acid alignments.
    another possibilites is providing the string "identity", which simulates the identity matrix
    if no matrix is provided, blosum62 is used
    """
    if not title:
      title=self.titles()[0]

    gap_score, insertion_score =-4, 0
    if not matrix:      matrix=load_blosum()

    seq=self.seq_of(title)
    index=self.remove(title, True)
    score=0
    for pos in range(len(seq)):
      aa1=upper(seq[pos])

      score_this_pos, gapped=0, 0
      for tit in self.titles():
        aa2=upper(self.seq_of(tit)[pos])
        if aa2=='-':
          gapped+=1
        elif aa1!='-':
          if matrix == 'identity':
            score_this_pos+=bool( aa1 == aa2 )
          else:
            score_this_pos+= blosum(aa1, aa2, matrix)

      relative_importance=(  1- (float(gapped)/len(self.titles()))  ) #the relative importance of this position is defined by the fraction of sequences of the profile that are not gapped
      if len(self.titles() ) == gapped:         average_score_this_pos=insertion_score
      else:                                     average_score_this_pos=float(score_this_pos)/( len(self.titles() ) - gapped )
      if aa1!='-':                              final_position_score = average_score_this_pos*relative_importance
      else:                                     final_position_score = gap_score*relative_importance
      #print pos, gapped, relative_importance, final_position_score
      score+=final_position_score
    score=float(score)/len(seq) #normalizing per length
    self.add(title, seq, index)
    return score    

  def average_conservation_score(self, matrix={}):
    """ This functions returns the average and the std deviation of the conservation scores of the sequences in the alignment. These values represent in a way how well defined (or how loose) is a profile alignment
    """
    scores=[]
    if not matrix:
      matrix=load_blosum()
    for title in self.titles():
      scores.append( self.conservation_score(title, matrix) )

    average=0
    for s in scores:
      average+=s
    average/=len(self.titles())

    std_deviation=0
    for s in scores:
      std_deviation += pow( (s-average), 2 ) 
    std_deviation/=len(self.titles())

    return average, std_deviation
  
  def sort(self, o_function=None, inplace=True):  
    """ This function can be used to change the order of the sequences in the alignment. similarly to list.sort, it accepts a sorting function which must accept two [title, seq] arguments. If the function is not specified, the alignment is sorted alphabetically by title. If inplace is False, the self alignment is left untouched and another alignment istance is returned instead.
    """
    if not o_function:
        def o_function(x, y):
          x_title, x_seq, y_title, y_seq=x[0], x[1], y[0], y[1]
          return cmp(x_title, y_title)
    if inplace:
      return self.order.sort(o_function, key=lambda x: (x, self.seq_of(x)) )
    else:
      a=alignment(self.diz)
      a.sort(o_function, True)
      return a

  def copy(self):
    return deepcopy(self)

  def columns(self, position, length=1, inplace=False, remove_empty_seqs=False):  
    """ columns(self, position, length=1, inplace=False, remove_empty_seqs=False)
    This function extracts a subalignment of a certain portion, starting from position (first is 0) and long length.
    """
    if self.length()<position+length:
      raise Exception, "alignment->columns ERROR: position+length called too high, > than length of alignment: "+str(position)+'+'+str(length)+'>'+str(self.length())
    if not inplace:
        a=self.copy()
    else:
        a=self; 
    a.reset_derived_data()
    for title in a.titles():
        a.diz[title]=a.seq_of(title)[position:position+length]
    if remove_empty_seqs:      a.remove_empty_seqs()
    if not inplace:
        return a

  def remove_empty_seqs(self):
    """ Delete those sequences which contain only gap in the alignment, and return the list of their title."""
    out=[]
    all_titles=list(self.titles())
    for t in all_titles:
#      if 'TR_hit3' in t:
#        print '*******************', [replace_chars(self.seq_of(t), '-', '')]
      if not replace_chars(self.seq_of(t), '-', ''):
        out.append(t)
        self.remove(t, True)
    if out: self.reset_derived_data()
    return out
    
  def replace_columns(self, position, length, subalignment, inplace=True):  
    """ This function replace a portion of the alignment, starting from position (first is 0) and long length, with another alignment whose titles must be identical
    """
    if self.length()<position+length:
      raise Exception, "alignment->replace_columns ERROR: position+length called too high, > than length of alignment: "+str(position)+'+'+str(length)+'>'+str(self.length())
    if not inplace:
        a=self.copy()
    else:
        a=self; self.reset_derived_data()
    if sorted(self.diz.keys())!=sorted(subalignment.diz.keys()):  #checking titles. 
      #write_to_file( join( sorted(self.diz.keys()), '\n'), 'a')  
      #write_to_file( join(sorted(subalignment.diz.keys()), '\n'), 'b')
      raise Exception, "alignment->replace_columns ERROR: the two alignments must have the same titles! "

    for title in self.titles():
        a.diz[title]=a.diz[title][:position]+subalignment.diz[title]+a.diz[title][position+length:]
    if not inplace:
        return a

  def concatenate_with(self, other_ali, inplace=False):  
    """ This function sum the self with another alignment with identical titles: the sequences are concatenated
    """
    if self.titles()!=other_ali.titles():  #checking titles. 
      raise Exception, "alignment->concatenate_with ERROR: the two alignments must have the same titles! "
    if not inplace:
        a=self.copy()
    else:
        a=self; self.reset_derived_data()
    for title in self.titles():
        a.diz[title]=a.diz[title]+other_ali.diz[title]
    if not inplace:
        return a

  def delete_columns(self, position, length, inplace=True):  
    """ This function delete subalignment of a certain portion, starting from position (first is 0) and long length.
    """
    if self.length()<position+length:
      raise Exception, "alignment->delete_columns ERROR: position+length called too high, > than length of alignment: "+str(position)+'+'+str(length)+'>'+str(self.length())
    if not inplace:
        a=self.copy()
    else:
        a=self; self.reset_derived_data()
    for title in self.titles():
        a.diz[title]=a.diz[title][:position]+a.diz[title][position+length:]
    if not inplace:
        return a

  def realign(self, program='mafft', inplace=True, verbose=False, char_replace='*U', protein=True, mafft_options=" --auto "):
    """ This function uses mafft to realign the sequences in the alignment. The alignment is returned, unless inplace==True, in this case the self object is modified. """
    if self.titles():
      test_writeable_folder(temp_folder, 'defined temp folder')
      #checking presence of stars, which will be lost when aligning with mafft
      forbidden_dict={}
      for c in char_replace: forbidden_dict[c]=1
      stars_chars={} #title -> list of positions of any forbibbed char (char_replace), 0 based
      for title in self.titles():
        seq_no_gap=nogap(self.seq_of(title))
        for index, char in enumerate( seq_no_gap ):
          if char in forbidden_dict: 
            if not stars_chars.has_key(title): stars_chars[title]=[]           
            stars_chars[title].append( [index, char]  )
        if stars_chars.has_key(title): self.set_sequence(title, replace_chars(self.seq_of(title),  forbidden_dict, 'X') )

      if program=='mafft':
        self.display(temp_folder+'realigning_with_mafft')
        mafft_options+=' --amino ' *int(bool(protein))
        bbash('mafft '+mafft_options+' '+temp_folder+'realigning_with_mafft > '+temp_folder+'realigning_with_mafft.aligned', verbose)
#        if int(bbash('ls -s '+temp_folder+'realigning_with_mafft.aligned').split()[0])<9:
          #file empty. Maybe the input sequence is aminoacid and the program mafft didn't understand it (often happens for small alignments)
#          bbash('mafft --auto --amino '+temp_folder+'realigning_with_mafft > '+temp_folder+'realigning_with_mafft.aligned', verbose)

        b=alignment()
        seq_index=0
        for title, seq in parse_fasta(temp_folder+'realigning_with_mafft.aligned'):                               
          old_full_title= self.titles()[seq_index]
          #assert           self.fill_title(title) == old_full_title
          b.add( old_full_title, seq)          
          seq_index+=1

        if not b.titles(): raise Exception, "ERROR realign with mafft failed! output alignment is empty."

        for title in stars_chars:
          seq=b.seq_of(title)
          for pos_seq, char in stars_chars[title]:
            pos_ali=b.position_in_ali(title, pos_seq+1)
            seq= seq[:pos_ali-1]+ char  +seq[pos_ali:]
          b.set_sequence(title, seq)
        if not b.titles():          raise AliError, "ERROR mafft failed in realigning!"
      else: raise Exception, "ERROR only mafft is currently supported"      
      #elif ... #add more programs if needed
      if not inplace:      return b
      else:                self.diz=b.diz; self.order=b.order; self.reset_derived_data()
      
  def realign_columns(self, position=0, length=0, inplace=True, input_list=[], protein=True ):   
    """ This function realigns only a portion of the alignment, keeping the rest identical. The alignment is returned, unless inplace==True, in this case the self object is modified.
 The list_input keyarg can be used to provide more than one column portion at the time (provide a list of [pos, length] elements ). in this case input position and length are not taken into account
     """
    if not input_list and length:      input_list=[[position, length]]
    if not inplace:        a=self.copy()
    else:                  a=self; self.reset_derived_data()
    shrinked_length=0
    for pos, l in input_list:
      if self.length()<pos-shrinked_length+l:        raise Exception, "alignment->realign_columns ERROR: position+length called too high, > than length of alignment: "+str(pos)+'+'+str(l)+'>'+str(self.length())
      cols=a.columns(pos-shrinked_length, l)
      removed_titles= cols.remove_empty_seqs()

      if (cols.nseq()>1):
        try:         
          #printerr( str(pos)+'-'+str(l)+', ready to realign!', 1)
          cols.realign(protein=protein)    #realigning on a try statement cause some particular columns may make crash the alignment program. In this case, we keep the cols as they were before realigning
        except:      pass
        for title in removed_titles:
          cols.add(title, '-'*cols.length() )

        a.replace_columns(pos-shrinked_length, l, cols)
        shrinked_length+=l-cols.length()
    if not inplace:
      return a

  def find_desert_columns(self, max_non_gaps=1, length_clusters=2, only_titles=[], join_neighbours=0):
    """ This functions parse the alignment looking for columns in which a single sequence has something different from a gap. Contiguos desert columns are grouped, and those whose length is > than length_clusters are returned. They are returned in the format: [[pos1, length1], [pos2, lenght2], ... ]   positions are 0-based!
max_non_gaps determines how many seqs can have a non gap in the columns returned. If it is a integer, it is taken as the maximum n of seqs, while if it is a float it is intended as a proportion over the total number of seqs (rounded by excess! careful this is to put the min value to 1 instead than to 0).
If only_titles is specified, it is necessary that the columns to realign have non-gaps only for these titles.

"""
    if type(max_non_gaps)==float:      max_non_gaps=int(max_non_gaps*self.nseq())+1
   
    out=[]
    if len(only_titles)!=1:
      n_last_columns_matching=0
      for pos in range(self.length()):
        non_gaps=0
        non_gaps_titles=[]
        for title in self.titles():
          if self.seq_of(title)[pos]!='-':
            non_gaps+=1
            non_gaps_titles.append(title)
        if non_gaps<=max_non_gaps and not only_titles or all([i in only_titles  for i in non_gaps_titles]):
          n_last_columns_matching+=1
        else:
          if n_last_columns_matching >= length_clusters:
            out.append([pos-n_last_columns_matching, n_last_columns_matching])
          n_last_columns_matching=0
      if n_last_columns_matching >= length_clusters:
        out.append([pos-n_last_columns_matching+1, n_last_columns_matching])

    ## joining neighbours
    if join_neighbours:
      removed_indexes=[]
      for region_index  in range(len(out[:-1])):
        start_region, length_region=out[region_index]
        end_region=start_region+length_region-1
        next_start, next_length=out[region_index+1]
        next_end=  next_start+next_length-1
        if next_start-end_region -1  <= join_neighbours:
        #updating the next one to include the current one. marking the next one to be removed afterwards
          removed_indexes.append(region_index)
          out[region_index+1]= [  start_region, next_end-start_region+1 ]
      for i in removed_indexes[::-1]: out.pop(i)
    #printerr(out, 1)  #DEBUG debug

    return out
  
  def shrink(self, only_titles=[]):
    """ This functions detects the desert columns clusters in the alignment and realigns them"""
    desert_columns=self.find_desert_columns(1, 2, only_titles=only_titles) #, join_neighbours=2)
    self.realign_columns(input_list=desert_columns)
    self.reset_derived_data()
    
  def conserved_residues(self, mode='count'):
    """ This function compute the number of columns in the alignment that have high conservation, that is to say, all aminoacids are identical or similar to each other according to similar_aas function.
    if mode=='count', the number of the columns is reported.
    if mode=='proportion', the proportion of the columns (respect to the length of the alignment) is reported.
    if mode=='positions', the list of positions is returned ( NB 0 based! )

    """ 
    positions=[];     n_conserved=0
    for pos in range( self.length() ):
      if self.is_column_conserved(pos):
        if mode=='positions':          positions.append(pos)
        n_conserved+=1
    if mode=='count':            return n_conserved
    elif mode=='positions':      return positions
    elif mode=='proportion':     return float(n_conserved)/self.length()
      
  def is_column_conserved(self, pos, mode='similar'):
    """ Utility used by previous function, conserved residues. This function checks if a column of the alignment is conserved, this meaning either that all characters must be identical (mode=='identical') or that they all must be all similar aminoacids, computed with function similar_aas (mode=='similar', default). NB position is 0 based."""
    chars_in_this_pos={}
    for title in self.titles():
      if mode=='identical':
        for c in chars_in_this_pos:
          if self.seq_of(title)[pos]!=c:            return False
      elif mode=='similar':
        for c in chars_in_this_pos:
          aa_this_pos=self.seq_of(title)[pos]
          if not c==aa_this_pos and not  similar_aas(aa_this_pos, c):            return False
      chars_in_this_pos[self.seq_of(title)[pos]]=1
    return True

  def sort_by_completeness(self, inplace=True):
    """ This function sort the titles by completeness of the sequence. The first ones will be the more complete ones, more suitable to be queries. if inplace == False, the ordered list of titles is returned instead of set to the .order attribute of this alignment. """
    cons_map=self.conservation_map()
    temp_diz_to_sort={}
    for title in self.titles():
        score=0
        for pos in range(len(cons_map)):
          importance_of_this_pos=0
          max_cons_lett= ['', 0]
          for k in cons_map[pos]:
            if cons_map[pos][k] > max_cons_lett[1] and k!='-':
              max_cons_lett= [k, cons_map[pos][k] ]

          if self.seq_of(title)[pos]!='-':            score+=max_cons_lett[1]
          else:                                       score-=max_cons_lett[1]
        temp_diz_to_sort[title]=   float(score) / self.length()
    sorted_titles=temp_diz_to_sort.keys()
    sorted_titles.sort(key=temp_diz_to_sort.get, reverse=True)
    if inplace:      self.order= sorted_titles
    else:            return      sorted_titles

  def set_sequence(self, title, sequence):
    """ Like add, but it is meant to change existing sequences. If the alignment does not have title, an expection is raised"""
    if not self.diz.has_key(title): raise Exception, "alignment->set_sequence ERROR the alignment does not have title: "+str(title)
    self.diz[title]=sequence
    self.reset_derived_data()      

  def transfer_alignment(self, other_ali, neutral_chars='UX*', dont_shrink=False):
    """ This function computes a common alignment between two existing alignments, exploiting common sequences. If the two alignments have no common titles, an exception is raised"""
    def find_prot_no_gap(names_dict, position, ali):
      for name in names_dict:
        if ali.seq_of(name)[pos]!='-':            return name
      return False
    this_ali=self.copy(); other_ali=other_ali.copy()
    diz_common={}  #contains the names of the common prot between this and other alignment
    for title_self in this_ali.titles():
      if title_self in other_ali.titles():    diz_common[title_self]=1
    if len(diz_common)==0:     raise Exception, 'alignment->transfer_aligment ERROR: no common proteins between the two alignments'

    position_map=[-1 for i in range(this_ali.length())]   #maps positions of this_ali to the positions of other_ali
    
    prot=diz_common.keys()[0]           
    pos=0 #this keeps track of the position in this_ali  ; pos_global keeps track of the position in other ali
    while pos <this_ali.length():
      if this_ali.seq_of(prot)[pos]=='-':              prot=find_prot_no_gap(diz_common, pos, this_ali)
      if not prot:    #none of the common proteins has something different from a gap in position: pos.
        pos_global=position_map[pos-1]+1   #becomes 0 if pos is 0, for construction
        position_map[pos] = pos_global
        for p_name in other_ali.titles():        
          new_sequence= other_ali.seq_of(p_name)[:pos_global]+'-'+other_ali.seq_of(p_name)[pos_global:]
          #if 'sps_Plasmodium_knowlesi_2007-2008/1-811' in p_name: print ">before: \n"+other_ali.seq_of(p_name)+'\n>after: \n'+new_sequence
          other_ali.set_sequence(p_name, new_sequence)
        prot=diz_common.keys()[0]
      else:    
        aa=this_ali.seq_of(prot)[pos]
        pos_global=position_map[pos-1]+1     #becomes 0 if pos is 0, for construction
        while position_map[pos]==-1:
            aa_global=other_ali.seq_of(prot)[pos_global]
            if aa_global!='-':
              if lower(aa_global)==lower(aa) or (aa_global in neutral_chars) or (aa in neutral_chars):
                      position_map[pos]=pos_global
              else:                  
                printerr( other_ali.seq_of(prot)[pos_global-4:pos_global+4] +' != '+ this_ali.seq_of(prot)[pos-4:pos+4] , 1)
                raise Exception, 'alignment->transfer_alignment ERROR AminoAcids dont correspond in pos '+str(pos)+' | '+str(pos_global)+' ('+aa_global+' != '+aa+') for title: '+ prot
            else:                pos_global+=1
      pos+=1
    
    #over. let's add initial and final gaps, and add the sequences to other_ali
    for p_name in this_ali.titles():                                    # ->up 
      seq=''    
      for pos in range(len(position_map)):
        if pos==0:
          for i in range(position_map[0]):    seq+='-'    #adding initial gaps
        else:
          for i in range(position_map[pos]-position_map[pos-1]-1):                seq+='-'
        seq+=this_ali.seq_of(p_name)[pos]
      for i in range(other_ali.length()-1 -position_map[-1]):  seq+='-'    #adding final gaps
      if not diz_common.has_key(p_name):          other_ali.add(p_name, seq )
    #correcting desert columns defect:
    if not dont_shrink:      other_ali.shrink()
    return other_ali
               
               
  def remove_redundancy_with_t_coffee(self, max_pair_identity=0.5, inplace=True):
    """ This function removes sequences keeping the most representative ones. It uses the seq_reformat +trim subroutine of t_coffee, therefore it should not be used on datasets too large.
    To avoid problems with t_coffee titles, the sequences are written in a temporary fasta file using their id (index in the alignment), and read again later.
    
    """
    temp_ali_h=open(temp_folder+'temp_ali_for_removing_redundancy.fa', 'w')    
    for index, title in enumerate( self.titles() ):
      print >> temp_ali_h, ">"+str(index)+'\n'+fasta(self.seq_of(title))
    temp_ali_h.close()

    bbash('t_coffee -other_pg seq_reformat -in '+temp_folder+'temp_ali_for_removing_redundancy.fa'+' -action +trim _aln_%%'+str(int(max_pair_identity*100))+'_ -output fasta_aln > '+temp_folder+'temp_ali_for_removing_redundancy.trimmed.fa')
    
    new_ali=alignment()
    for title, seq in getfastalite(open(temp_folder+'temp_ali_for_removing_redundancy.trimmed.fa', 'r')):
      original_index=int(title)
      new_ali.add(self.titles()[original_index], seq)
      
    if inplace:      self.reset_derived_data(); self.__dict__=deepcopy(new_ali.__dict__); 
    else:            return new_ali        

  def positions_of_similar_aas(self):
    """ This function returns the positions (1-based) of the columns which contains aminoacids all identical OR similar to each other according to the similar_aas function. Thought for pairwise alignment, it works also with larger ones. Positions where gaps are present are never returned since they are not considered similar to anything """
    if len(self.titles())>=1: 
      outlist=[]
      for pos in range(self.length()): #pos is 0 based
        main_aa=self.seq_of(self.titles()[0])[pos]
        similar_aa_in_this_pos=True
        for title in self.titles()[1:]:
          aa=self.seq_of(title)[pos]
          similar_aa_in_this_pos = similar_aa_in_this_pos and (  similar_aas(main_aa, aa)  or main_aa==aa )
        if similar_aa_in_this_pos:           outlist.append(pos+1)
      return outlist

  def convert_sequences(self, sfunction, *args, **keyargs):
    """ This apply a certain function on every sequence of the alignment. Can be useful to convert the alignment in uppercase or lowercase: for example, provide the function upper from the module string as sfunction"""
    for title in self.titles():
      try:       new_seq=  sfunction(   self.seq_of(title) , *args, **keyargs)
      except:    printerr("alignment->convert_sequences ERROR can't use provided function"); raise
      self.set_sequence(title, new_seq)

  def remove_duplicated_transcripts(self, identity_threshold, naming_function=None, inplace=False, verbose=False):
    """ This function is thought to remove redundant predictions coming from a p2g prediction on a set of RNAs like a EST database.
    The function forms clusters of the alignment using the function self.clustering and the specified identity threshold (mode=0 is used --> dont_count_terminal_gaps, and no reclustering)
    If then derives a single sequence from each cluster, computing a consensus for each column to exploit all information. The id of the longest sequence in the cluster is used. 
    optionally, you can provide a naming function which must have the form:   function(cluster_alignment) -> new_title
    """

    out_ali=   self.__class__() #creating a new instance of the same class of self
    clusters= self.clustering(threshold=identity_threshold, dont_remove_empty_columns=True)
    for cluster_ali in clusters:
      if len(cluster_ali.titles())>1:
        if naming_function is None:   cluster_title=        sorted(  cluster_ali.titles(), key=lambda x: len( nogap(cluster_ali.seq_of(x) )), reverse=True   )[0]
        else:                         cluster_title=        naming_function(cluster_ali)
        cluster_seq=''
    
        cmap=cluster_ali.conservation_map(dont_save=True)
        for pos in range(cluster_ali.length()):
          sorted_aa_this_pos =  sorted( cmap[pos].keys(), key= lambda x: cmap[pos].get, reverse=True) #sorted by representation
          if len(sorted_aa_this_pos)>1 and sorted_aa_this_pos[0]=='-':    consensus_this_pos=sorted_aa_this_pos[1]
          else:                                                           consensus_this_pos=sorted_aa_this_pos[0]
          cluster_seq+=consensus_this_pos

        if verbose: 
          printerr(  join([ t.split()[0] for t in cluster_ali.titles()], ', ')+' --> '+cluster_title, 1   )

      else:
        cluster_title=cluster_ali.titles()[0]
        cluster_seq= cluster_ali.seq_of( 0 )

      out_ali.add(cluster_title, cluster_seq)    

    #sort ali?  
    
    if inplace: self.__dict__ = out_ali.__dict__
    else:       return out_ali    




  def clustering(self, threshold, mode=0, reclustering=False, reclustering_factor=2, dont_remove_empty_columns=False, outclass=None, min_overlap=10):
    """ This function analyze the self alignment and split it in clusters of similar sequences. 
    Threshold (sequence identity) is used when doing pairwise comparisons to decide whether or not to cluster the two sequence together.
    mode defines how the sequence identity is computed: 0 [DEFAULT] -> dont_count_terminal_gaps    ; 1 -> dont_count_gaps   ; 2 -> count gaps
    you can also directly provide a sequence_identity_hash-like data structure as mode, in this case you can provide any measure of identity-
    for details on howo the sequence identity is computed, see function alignment.sequence_identity_hash
    If reclustering==True, then an additional procedure is performed after clustering for the clusters formed by a single sequence: a second step of clustering is performed, with a more permissive threshold (threshold/reclustering_factor)
    by default, empty columns are removed from every output alignment, to turn this off use dont_remove_empty_columns=True
    the list of output alignments is returned. Their class is alignment by default, but you can choose it with the outclass argument (anyway, no attributes are copied)
    """

    output=[] #will contain: list of output alignments
    clusters = {}; clusters_next_id = 0
    if outclass is None: outclass=alignment

    if self.nseq()<2: 
      o=outclass()
      for title in self.titles():
        o.add( title, self.seq_of(title) )
      return [o]

    # building identity_hash {title:{other_titles_in_alignment:identity}}
    if mode == 1:      h = self.sequence_identity_hash(dont_count_gaps=True)
    elif mode == 2:    h = self.sequence_identity_hash()
    elif mode==0:      h = self.sequence_identity_hash(dont_count_terminal_gaps=True)
    else:              h=mode

    ## main program algorythm
    for index_1 in range(len(self.titles())):
      title_1 = self.titles()[index_1]
      
      if clusters.has_key(title_1):
        this_cluster_id = clusters[title_1]
      else:
        this_cluster_id = clusters_next_id
        clusters_next_id += 1
        clusters[title_1] = this_cluster_id

      for index_2 in range(index_1+1, len(self.titles())):
        title_2 = self.titles()[index_2]
        if not min_overlap or  len([ True     for p in range(self.length())   if  self.seq_of(title_1)[p]!='-' and self.seq_of(title_2)[p]!='-'   ]) >= min_overlap: # checking in how many columns we have some non gap in both sequences. if the number is very low, we can have high seq id even idf the sequences don't look alike
          if h.get(title_1, title_2) >= threshold:
            if clusters.has_key(title_2) and not clusters[title_2]==this_cluster_id : # title_2 was already in a cluster. let's update all other titles to be in this cluster
              clusters_to_update = clusters[title_2]
              
              for t in self.titles():
                if clusters.has_key(t) and clusters[t] == clusters_to_update:
                  clusters[t] = this_cluster_id
            else:
              clusters[title_2] = this_cluster_id

    c = {}#{c_id:alignment}
    for t in self.titles():
      c_id = clusters[t]
      if not c.has_key(c_id):
        c[c_id] = outclass()
      c[c_id].add(t, self.seq_of(t))

    if reclustering:
      solitaire_seqs={} #sequences left alone (a cluster of a single sequence)
      #print """The following sequences don\'t get the minimum threshold.
    #For them, a second clustering is performed reducing the threshold in one half.
    #"""
      for k in c.keys():
        if c[k].nseq() == 1:
          for title in c[k].titles():
            solitaire_seqs[title]=1
            max_title=''; max_seq_id=0.0
            for title2 in self.titles():
              if title2!=title: 
                if h.get(title, title2) > max_seq_id:
                  max_seq_id=h.get(title, title2)
                  max_title=title2
            if h.get(title, max_title) >= threshold/float(reclustering_factor) and not solitaire_seqs.has_key(max_title):
              c[clusters[max_title]].add(title, self.seq_of(title))
              #print ' --> '+title.split()[0].ljust(40)+' clustered with   '+max_title.split()[0].ljust(50)+'identity '+str(h.get(title, max_title))
              del c[clusters[title]]
      
    if not dont_remove_empty_columns:
      for k in c.keys():      c[k].remove_useless_gaps()      

    sorted_c_ids = c.keys()
    def nseq_for_key(k):      return c[k].nseq()
    sorted_c_ids.sort(key=nseq_for_key, reverse=True)    
        
    for final_cluster_id in range(len(sorted_c_ids)): # the final_cluster_id is actually not used.
      old_c_id = sorted_c_ids[final_cluster_id]
      output.append(c[old_c_id])

    return output
     
  def remove_redundancy(self, threshold, mode=0, dont_remove_empty_columns=False, outclass=None, min_overlap=10, inplace=True, silent=True, correspondance_hash={}, choose_repr=None):
    """ Modify in place or returns a copy of the alignment removing redundancy, so that in output all sequences will have sequence identity > threshold (between 0.0 and 1.0). Initially alignment is split into clusters with function clustering, then a single sequence is chosen for each cluster; normally this is done by maximizing the average sequence identity with the sequences removed in favor of this (the sequence closest to each cluster baricentre act as representative). But! you can instead provide a function to choose_repr ; this must be a function that accepts an alignment (self) and a list of titles, and returns a single title. If it returns None, the function will use the built-in representative choice procedure for that particular cluster."""
    # building identity_hash {title:{other_titles_in_alignment:identity}}
    if   mode == 1:      h = self.sequence_identity_hash(dont_count_gaps=True)
    elif mode == 2:      h = self.sequence_identity_hash()
    elif mode == 0:      h = self.sequence_identity_hash(dont_count_terminal_gaps=True)
    clusters= self.clustering(threshold, mode=h, dont_remove_empty_columns=True, min_overlap=min_overlap)  ## computing all clusters; providing identiti hash to avoid computing it twice
    out_titles=[]
    for c in clusters:
      if len(c.titles())>1:
        best_title=None
        if not choose_repr is None:
          best_title   = choose_repr (self, c.titles() )
        if best_title is None :   #either choose_repr is None, or it returned None in this case so we want to use the default function
          average_seq_identities_within_cluster={}
          for t in c.titles():
            average_this_t= sum(  [ h.get(t, t2)     for t2 in c.titles()   if t!=t2 ]) / float(len(c.titles())-1)
            average_seq_identities_within_cluster[t]=average_this_t
          these_titles_sorted= sorted( c.titles(),   key=average_seq_identities_within_cluster.get, reverse=True   )
          #for t in these_titles_sorted:          print t.split()[0], average_seq_identities_within_cluster[t]
          best_title   = these_titles_sorted[0]
        other_titles = [t for t in c.titles() if t!=best_title]
      else: 
        best_title=c.titles()[0]
        other_titles=[]
      out_titles.append([best_title, other_titles])
      #o.add(best_title,  self.seq_of(best_title) )

    for out_t, other_ts in out_titles:
      for discarded_t in other_ts: 
        correspondance_hash[discarded_t]=out_t
        if not silent:          write('REMOVED: '+ discarded_t.split()[0]+' >KEPT: '+out_t.split()[0], 1)

    if inplace:
      out_titles_h={}
      for  t, other_ts    in out_titles: out_titles_h[t]=1
      titles= self.titles()
      for t in titles:
        if not t in out_titles_h:  self.remove(t, True)
      if not dont_remove_empty_columns:    self.remove_useless_gaps()
    else:
      if outclass is None: outclass=alignment    
      o=outclass()
      for t, other in out_titles:          
        o.add(t,  self.seq_of(t) )
      if not dont_remove_empty_columns:    o.remove_useless_gaps()
      return o      
    

  def build_tree(self,    folder=None,  tree_class=None,  trimal_options=' -phylip -gt 0.1 -cons 33.33 ', phylogeny_options= '-c /home/mmariotti/software/selenoprofiles/libraries/salva_pipeline.config.for_selenoprofiles'    ):
    """ Uses Salva pipeline to build a tree. Note: use only for amino acid sequences.
    A folder dedicated to compute the phylogeny is created  (folder argument -- done in temp folder if not provided).
    it the tree_class argument is not defined (as default), the path to newick tree file computed is returned. Otherwise, one may provide a tree class (such as ete2.Tree), and this will be initialised with the mentioned tree file     """
    
    if self.nseq()<3: raise Exception, "ERROR can't build a tree with less than 3 sequences!"
    if folder is None: folder=temp_folder+'building_tree'
    folder=Folder(folder) # creating if necessary, adding "/" 
    input_of_trimal_filename =                          folder+'ali.'+'.raw_ali'
    output_of_trimal_filename =                         folder+'ali.'+'.alg.clean'
    fasta_copy_of_output_of_trimal_filename =           folder+'ali.'+'.aligned_fa'
    ungapped_fasta_copy_of_output_of_trimal_filename =  folder+'ali.'+'.fa'
    ali_temp=alignment()
    for title in self.titles():
      ali_temp.add(title, self.seq_of(title) )

    ali_temp.convert_sequences(replace_chars, '*', 'X')
    ali_temp.remove_useless_gaps()
    ali_temp.display(input_of_trimal_filename)
    bbash('trimal -in '+input_of_trimal_filename+' -out '+output_of_trimal_filename+' '+trimal_options)
    bbash('trimal -in '+output_of_trimal_filename+' -out '+fasta_copy_of_output_of_trimal_filename+' -fasta ')
    write_to_file( join([">"+t+'\n'+nogap(s)       for t, s in parse_fasta(fasta_copy_of_output_of_trimal_filename)], '\n')      , ungapped_fasta_copy_of_output_of_trimal_filename)  
    bbash('cd '+folder+' ; ReconstructingTrees.py '+phylogeny_options+' -f '+base_filename(output_of_trimal_filename)+' -d ./ ')

    ml_tree_files=bbash('find '+folder+' -name "'+'ali.tree.ml.*.nw'+'"').split('\n')
    if len(ml_tree_files)>1:
      rank_file=folder+'ali.tree.rank.ml'
      ranks={}
      for line in open(rank_file, 'r'):
        ranks[line.split()[0]]=float(line.split()[1])
      best_model=sorted(ranks.keys(), key=ranks.get)[-1] #since they are negative, we want the highest one, the closest the zero
      ml_tree_file=folder+'ali.tree.ml.'+best_model+'.nw'
    else:
      ml_tree_file=ml_tree_files[0]

    tree_link=folder+'ali.best_tree'
    bbash('ln -s '+base_filename(ml_tree_file)+' '+tree_link)   
  
    if tree_class is None: return tree_link
    return tree_class( tree_link ) 


      
def complete_word(uncomplete_word, dictionary, word_match=False):
  """given a word, returns the entry in the dictionary from where this word where cut. If it is not unique returns -1, if none is matching returns False
  ex: ('cat' {'carpet':1, 'catastrophe':1, 'bunny':1}) -->return 'catastrophe'
  if word_match==True: when the title would be called as non-unique, an additional check is performed: if a single title has a space or a tab right after the match, this is considered a good match and it is returned.   ex: ('dog' , {'doggie':1, 'dogs':1, 'dog one':1}) --> 'dog one' is returned. if no word_match option is set, this would return -1
  """
  if dictionary.has_key(uncomplete_word): return uncomplete_word
  title_length=len(uncomplete_word)
  title_that_match=''
  for t in dictionary.keys():
#    if len(t)>= title_length:
#      'nothing'
#      print t[:title_length]
    if t[:title_length] == uncomplete_word:
      if title_that_match:
        #already on!       not unique
        if word_match:
          #unless: word_match is on and a single title has a space or a tab right after the match. this match has more value than the other one.
          if    t[title_length] in '\t ' and not title_that_match[title_length] in '\t ':           title_that_match=t
          elif  t[title_length] in '\t ' and title_that_match[title_length] in '\t ':               return -1
        else:   return -1
      else:
        title_that_match=t
  if not title_that_match:    return False
  else:                       return title_that_match


def remove_items_from_list(alist, item_list, inplace=False):
  """ Utility to remove efficiently multiple items from a list. The arguments are the input list and the list of items to remove. If inplace==True, the operation is performed directly on the input list, other wise a new list is returned"""
  if not inplace: outlist=list(alist)
  else: outlist=alist
  items_hash={}
  for item in item_list: items_hash[item]=1
  for index in range(len(outlist)-1, -1, -1):
    if items_hash.has_key(outlist[index]): outlist.pop(index)
  return outlist

def remove_items_from_list_with_indexes(alist, index_list, inplace=False, index_are_ordered=False):
  """ Same as the remove_items_from_list function , but the list of indexes to remove is provided instead of the item list. the function orders the list of indexes as first step. If they are already ordered, you can save this time by setting index_are_ordered = True """
  if not inplace: outlist=list(alist)
  else: outlist=alist
  if index_list:
    if not index_are_ordered:    index_list.sort()
    if index_list[-1] >= len(outlist): raise Exception, "remove_items_from_list_with_indexes ERROR one or more indexes are > than the length of the list: "+str(index_list)+' > '+str(len(outlist))
    for item_index in range(len(outlist)-1, -1, -1):
      if item_index== index_list[-1]:     #item found
        index_list.pop(-1)
        outlist.pop(item_index)
    if index_list: raise Exception, "remove_items_from_list_with_indexes oops"
  return outlist



def getfastalite(cfile,order=1):
  line=cfile.readline()
  seq=''
  ord_list=[]    #[ [title1, seq1], [title2, seq2] ... ]
  diz={}      #{ title: seq   }
  while line:
    if line[0]=='>':
      if seq:    #not first line
#        print [title]
        if order==1:
          ord_list.append(  [title, seq]   )
        else:
          diz[title]=seq

      title=del_white(replace(line[1:], '\n', ''))
      seq=''
    else:
      seq+=replace( replace(line, '\n', ''), '#', '')

    line=cfile.readline()

  if order==1 and seq:
    ord_list.append(  [title, seq]   )

  elif seq:
    diz[title]=seq  

  if order==1:
    return ord_list
  else:
    return diz





def getfastaaln(cfile,order=1): #clustalW format
  line=cfile.readline()
  seq=''
  ord_list=[]    #[ [title1, seq1], [title2, seq2] ... ]
  diz={}      #{ title: seq   }
  while not line or line[0]=='#' or line.startswith('CLUSTAL'):
                line=cfile.readline()
            
#        while not finished:
  cont_seq=0
  while line:
    while line!='\n' and line.strip() and not line.strip()[0] in '*.:':
      title= join(line.split()[:-1], ' ')
      if  title and  line.split()[0]!='cons':
#        print [line]
        if len(ord_list)<=cont_seq:
          ord_list.append([title, ''])

        ord_list[cont_seq][1] += line.split  ()[-1]
        cont_seq+=1

      line=cfile.readline()
    line=cfile.readline()
    cont_seq=0

  if order==1:
    return ord_list
  else:
    return diz






def getfasta(text, order=1):
    """text=file.readlines(); the file should contain fasta sequences like ">title \n sequence".
    If order=1 it returns an ordered list such as [[title1, seq1], [title2, seq2], ...] ; 
    if order=0 it returns a dictionary containing {key: title, value: related sequence}.

    Works with:     fasta files (also alignments), commented fasta files, extended fasta format, modeller alignment format
                yap output (uses external function)**** not reported
    doesn't work with clustalW alignments.
    """
    if len(text)==1:
        text=text[0].split('\r')
    text.append('!...END...!')
    alignchars=AA_LETT+'-*' #allowed characters in alignments; the * is meant to be found only at the end of the seq
    cont=0      #number of proteins (number of '>title')
    pos={}      #key: protein_number (the cont associated to it); value: cline of the line with >title
    seq={}          #output
    ord_list=[]     #output ordered
    cline=0     #contatore linea
    for line in text:               #fill pos (dictionary)
        if line!='':
            if line[0]=='>':
                pos[cont]=cline
                cont=cont+1
        cline=cline+1    
    for n in pos:
        name=replace(replace(text[pos[n]][1:],'P1;',''),'\n','')         #replacing P1 is to get Modeller format alignments         
#        print 'name:'+name
        seq[name]=''
        done=0
        i=0
        ord_list.append([name,''])
        ctrl=0  #volte che trova la sequenza dopo un >
        while done!=1:
            seqtoadd=replace(text[pos[n]+i],'\n','')
#            print 'seqtoadd:'+seqtoadd
            if set(alignchars).issuperset(set(seqtoadd)) and seqtoadd!='':
                if seqtoadd[-1]=='*':
                    seqtoadd=seqtoadd[:-1]

#                print 'ordlist='+str(ord_list)
                ord_list[n][1]=ord_list[n][1]+seqtoadd
                seq[name]=seq[name]+seqtoadd
    
    
    
                ctrl=ctrl+1
            else:
                if ctrl!=0:
                    done=1
            i=i+1
    if order==1:
        return ord_list
    return seq    

reverse_complement_diz={'A':'T', 'T':'A','G':'C','C':'G',  'N':'N', 'X':'X',   'a':'t', 't':'a','g':'c','c':'g',}
def reverse_complement(seq):
  out=''
  for i in range(len(seq)):
    if reverse_complement_diz.has_key(seq[len(seq)-1-i]   ):      out+=reverse_complement_diz[seq[len(seq)-1-i]   ]
    else:                                                         out+=seq[len(seq)-1-i]  
  return out

def smith_waterman(seq1, seq2, gap_open=-4, gap_extension=-2, matrix={}, pssm=[]):
  """ This function computes the smith water alignment (local) between the two sequences in input, using the given gap_open and gap_extension parameters. The matrix used can be provided as an argument, or blosum62 called with blosum(a, b[, matrix]) is used. Alternatively, the pssm argument can be provided. This must be a list of hash with scores for each aminoacid.
  Returns:  [aligned_seq1, aligned_seq2, start_1, end_1, start_2, end_2, score]  
  notice that sequences returned are typically shorter than the input ones, as the alignment is local. indexes returned can be used to complete them.
  if no good scoring alignment is computed (all scores <0), then None is returned.
"""
  m,n =  len(seq1),len(seq2)#length of two sequences

  #generate DP table and traceback path pointer matrix
  #the DP table  
  score=    [ [ 0  for i in range(n+1) ] for j in range(m+1) ]
  pointer=  [ [ 0  for i in range(n+1) ] for j in range(m+1) ]

  P=0;
  #initial maximum score in DP table
  max_score=P;

  #calculate DP table and mark pointers
  for i in range(1,m+1):
    for j in range(1,n+1):
      #score up
      if pointer[i-1][j]==1:        penalty=gap_extension
      else:                         penalty=gap_open        #can be only 0; 2 is impossible for construnction to be selected
      score_up=score[i-1][j]+penalty;
      #score_down
      if pointer[i][j-1]==2:      penalty=gap_extension
      else:                       penalty=gap_open        #can be only 0; 1 is impossible for construnction to be selected
      score_down=score[i][j-1]+penalty;
      #score diagonal
      if not pssm:
        match_score= blosum(seq1[i-1],seq2[j-1], matrix);
      else:
        if len(pssm)<=[i-1]: raise Exception, "ERROR the pssm provided is not long enough... position requested: "+str((i-1))
        elif not pssm[i-1].has_key(seq2[j-1]): raise Exception, "ERROR key error in pssm. position: "+str((i-1))+' does not have a value for: '+seq2[j-1]
        match_score = pssm[i-1][seq2[j-1]]

      score_diagonal=score[i-1][j-1]+match_score

      score[i][j]=max(0,score_up,score_down,score_diagonal);
      if score[i][j]==0:
        pointer[i][j]=0; #0 means end of the path
      if score[i][j]==score_up:
        pointer[i][j]=1; #1 means trace up
      if score[i][j]==score_down:
        pointer[i][j]=2; #2 means trace left
      if score[i][j]==score_diagonal:
        pointer[i][j]=3; #3 means trace diagonal
      if score[i][j]>=max_score:
        max_i=i;
        max_j=j;
        max_score=score[i][j];
  #END of DP table

  align1,align2='','';#initial sequences
  i,j=max_i,max_j;#indices of path starting point

  #traceback, follow pointers
  while pointer[i][j]!=0:
    if pointer[i][j]==3:
      align1=align1+seq1[i-1];
      align2=align2+seq2[j-1];
      i=i-1;
      j=j-1;
    elif pointer[i][j]==2:
      align1=align1+'-';
      align2=align2+seq2[j-1];
      j=j-1;
    elif pointer[i][j]==1:
      align1=align1+seq1[i-1];
      align2=align2+'-';
      i=i-1;
  #END of traceback

  align1=align1[::-1];#reverse sequence 1
  align2=align2[::-1];#reverse sequence 2

  if not align1: return None
  return [align1, align2, i+1, max_i, j+1, max_j, score[max_i][max_j]]


def correct_sequence(seq1, seq2, special_chars=['U']):
  """ This function corrects a sequence in input using a second reference sequence which is identical except in positions where seq1 has special chars (typically U). 
  seq1 is returned with all such mismatch positions changed into the letter found in seq2. 
  If more mismatches are present (not only aligned to special chars) an exception is raised.
  
  The sequences can contain gaps, which are ignored (those in seq1 are returned in the same positions)
  """
  index_no_gap=0; index2=-1 

  for index1 in range(len(seq1)):
    aa1=seq1[index1]
    if aa1!='-':
      index2+=1
      aa2=seq2[index2]
      while aa2=='-':
        index2+=1
        aa2=seq2[index2]
      if aa1!=aa2:
        if aa1 in special_chars:              seq1=seq1[:index1]+aa2+seq1[index1+1:]
        else: raise Exception, "correct_sequence ERROR the sequences are different not only for special characters: "+seq1+'  !=  '+seq2 
  return seq1

def dereference(filename):
  """ This function returns the absolute path of the destination of a symbolic link (if a symbolic link is linked to another symbolic link, it goes all the way to the last file) .
  If the file provided is not a symbolic link, its absolute path is simply returned. If it does not exists, an exception is raised.
  """
  if not is_file(filename):
    raise Exception, "dereference ERROR file: "+filename+' was not found '
  filename=abspath(filename)
  out_ls=bbash('ls -l "'+filename+'"')
  while out_ls.split()[-2]=='->':
    destination=  out_ls.split()[-1]
    if destination[0]!='/':
      filename=directory_name(filename)+'/'+destination
    else:
      filename=destination
    out_ls=bbash('ls -l '+filename)    
  return filename  

def fileid_for_temp_folder(filename):
  """ This returns the absolute and dereferenced path for a filename, with / and spaces replaced with _ """ 
  return replace_chars(dereference(filename), '/ ', '_')

def fastaindex(target_genome, silent=0, force=0):
  """Run fastaindex from the exonerate package on the target_genome and return the path to the index file, or just return the path if the file is found (unless force=1)
  If the user has no writing permissions in the folder where the genome file is, the index is created in the temporary folder and it is returned.
  When the index is created, a message is printed using functio write() unless silent=1
    """
  if '.' in base_filename(target_genome):    index_file=join( target_genome.split('.')[:-1],'.'   )+'.index'
  else:                                      index_file=target_genome+'.index'
  if is_file(index_file) and not force: return index_file
  temp_index_file=temp_folder+fileid_for_temp_folder(target_genome)+'.index'
  if is_file(temp_index_file) and not force: return temp_index_file   # again this command if the file is in the temporary folder
  cmnd='fastaindex '+target_genome+' '+temp_index_file
  if not silent:    printerr( 'indexing the nucleotide database '+target_genome  , 1)
  bbash(cmnd)
  try:        
    test_writeable_folder( directory_name(index_file) )
    bash('mv '+temp_index_file+' '+index_file)
    return index_file
  except:         return temp_index_file
  
def formatdb(target_file, is_protein=False, silent=0):
  """ Utility to format for blast with formatdb"""
  if not silent:    printerr("attempting to format database: "+target_file+" (will crash if you don't have permissions) ", 1 )
  test_writeable_folder(temp_folder, 'temp_folder'); test_writeable_folder( directory_name(target_file)  )
  temp_folder_for_formatting= temp_folder+'formatting_database'
  bbash('mkdir '+temp_folder_for_formatting+' && cd '+temp_folder_for_formatting+' && ln -s '+abspath(target_file)+' .  &&  formatdb -i '+base_filename(target_file)+' -o T -p '+{False:'F', True:'T'}[bool(is_protein)]+' && mv '+base_filename(target_file)+'.* '+directory_name(abspath(target_file)))

def fastafetch(split_folder, chromosome, target_genome, verbose=0, chars_not_allowed=[':']):
  """fecthing chromosome routine. File are fetched to a file named after the fasta title. chars_not_allowed is an iterable with characters which cannot appear in the output filename  """
  target_genome=abspath(target_genome)
  species_subfolder=Folder(split_folder+replace_chars(target_genome, '/', '_')) 
  final_filename=species_subfolder+chromosome
  
  for char in chars_not_allowed:
    if char in final_filename: final_filename=  replace_chars( final_filename , char, '{Ch'+str(ord(char))+'}')
  
  temp_filename=temp_folder+'fetching_chromosome.fa'
  if is_file(final_filename):    'file already there'
  else:  
    #printerr('*** '+fastaindex(target_genome), 1)
    cmnd='fastafetch '+target_genome+' '+fastaindex(target_genome)+' "' +chromosome +'" > "'+  temp_filename+'"' 
    service( '  ...fetching chromosome: '+chromosome )
#    debug(cmnd)
    b=bash(cmnd, verbose)
    if b[0]==256 and "Could not open fasta index" in b[1]:      #db not indexed? i'll do it!
      fastaindex(target_genome)
      bbash(cmnd, verbose)  #trying again to fetch the chromosome after indexing the db
    elif b[0]==256 and "Could not find identifier" in b[1] and "|" in chromosome: #trying adding a pipe symbol in the end. This is because blast sometimes removes that.
      cmnd='fastafetch '+target_genome+' '+fastaindex(target_genome)+' "' +chromosome +'|" > "'+  temp_filename+'"' 
      b=bash(cmnd, verbose)
    elif b[0]==256 and '_' in chromosome and  is_number( chromosome.split('_')[0], True):  ##SOME BLASTOUT DONT HAVE >GI BUT ARE LIKE >N_fastafilename. THIS IS FOR THEM. 
      nseq_index=chromosome.split('_')[0]
      #chromosome=nseq_index
      service( '  ...using as index of sequence the number in blast output:  '+nseq_index +' from header: '+chromosome )
      cmnd='get_sequence_numberN.g -v NSEQ='+nseq_index+' '+target_genome+'  > '+  temp_filename
      #debug(cmnd)
      bbash(cmnd, verbose)  
    #checking if there is sequence in temp_folder+'tchrom' 
    if bbash('head -n2 "'+  temp_filename+'"')=='':      raise Exception, 'ERROR fetching chromosome; command_line: '+cmnd
    bbash('mv '+temp_filename+' "'+final_filename+'"' )
  return  final_filename ######## NB different from the function in profiles_classes


def fastasubseq(subj_file, start, clength, out_file, pipecommand='', warning=False ):  #start can be <0, in that case it becomes 0; lenght can be > than the nts at the right of start, in that case... ###NB starts with 0
  """Utility to subseq fasta sequences. pipecommand is used for expert use: if you want to pass results through a pipe before going in out_file, you can use this. Also, if you want to append results instead of writing, use pipecommand='>' (so in the final command line it will appear >>)
  """
  
  start=max(start, 0)
  cmnd='fastasubseq "'+subj_file+'" '+str(start)+ ' '+str(clength)+" "+pipecommand+"> "+out_file
  ss = bash(cmnd)
  try:
    if ss[0]!=0 and "Subsequence must end before end" in ss[1]:
      old_clength=clength
      clength= int(   ss[1].split('\n')[0].split('(')[1][:-1]     ) - start
      cmnd='fastasubseq "'+subj_file+'" '+str(start )+ ' '+str(clength)+" "+pipecommand+"> "+out_file
      ss=bash(cmnd)
      if warning: printerr('Fastasubseq: '+str(old_clength)+' bp not available, cutting the first '+str(clength), 1)
    if ss[0]!=0:
      raise Exception, "COMMAND "+cmnd+" ERROR in fastasubseq: \""+ss[1]+"\""
    return clength    #returning message for how much the sequence was actually cut
  except Exception, e:
    raise Exception, "COMMAND "+cmnd+" ERROR in fastasubseq: \""+ss[1]+"\""
    

def fasta(string, char_per_line=60):
    """ gets a sequence or a fasta file and return the same thing properly formatted
    """
    out=''
    seq=''
    for cline in string.split('\n'):
      if cline[0]=='>':
        if seq:
          for i in range(    (len(seq)-1) / char_per_line   +1):
            out+= seq[i*60:(i+1)*60] +'\n'

        out+= cline+'\n'
        seq=''
      else:
        seq+=cline
    for i in range(   (len(seq)-1) / char_per_line   +1):
      out+= seq[i*60:(i+1)*60]+'\n'

    return out[:-1]



def getfastaxff(text, order=1):     #
    """text=file.readlines(); text (type= list) contains fasta sequences like ">title \n sequence \n #secondary str".
    If order=1 it returns an ordered list such as [[title1, seq1, ss1], [title2, seq2, ss2], ...] ;
    if order=0 it returns a dictionary containing {key: title, value: [related sequence }."""
    text.append('!...END...!')
    alignchars=uppercase+'-*#0. ' #allowed characters in alignments
#    print alignchars
    cont=0      #number of proteins (number of '>title')
    pos={}      #key: protein_number (the cont associated to it); value: cline of the line with >title
    seq={}
    ord_list=[]
    cline=0     #contatore linea
    for line in text:
        if line!='':
            if line[0]=='>':
                pos[cont]=cline
                cont=cont+1
        cline=cline+1
    for n in pos:
        name=replace(replace(text[pos[n]][1:],'P1;',''),'\n','').split(' ')[0]    #for modeller file capting
#        print 'name:'+name
        seq[name]=['','']
        done=0
        i=0
        ord_list.append([name,'',''])
        ctrl=0  #volte che trova la sequenza dopo un >. se trova la str secondaria fa *-1
        while done!=1:
            seqtoadd=replace(text[pos[n]+i],'\n','')
            if set(alignchars).issuperset(set(seqtoadd)) and seqtoadd!='':
                if seqtoadd[-1]=='*':
                    seqtoadd=seqtoadd[:-1]
#                print 'ordlist='+str(ord_list)
                if seqtoadd[0]=='#':
                    seqtoadd=seqtoadd[1:]
                    ctrl=ctrl*-1
#                print 'seqtoadd:'+seqtoadd
                if ctrl<0:                    
                    ord_list[n][2]=ord_list[n][2]+seqtoadd
                    seq[name][1]=seq[name][1]+seqtoadd
                else:
                    ord_list[n][1]=ord_list[n][1]+seqtoadd
                    seq[name][0]=seq[name][0]+seqtoadd
                    ctrl=ctrl+1
            else:
                if ctrl<0:
                    done=1
            i=i+1
    if order==0:
        return seq
    return ord_list


#functions for gene class
def are_overlapping_exons(exon1, exon2, strand='+', check_phase=True):
  """the input here are exon_o: [start, stop, phase] where start is always > than stop"""
  if exon1[0] > exon2[1]:  # e1.start>e2.end
    return False
  if  exon2[0] > exon1[1]: 
    return False
  if not check_phase or not strand:
    return True

  if strand=='+' and mod3(exon1[0]-exon1[2]) != mod3(exon2[0]-exon2[2]):   ###frame checking!!!     (pos_start of exon + frame)  mod3 must be the same for both  --> equivalent: (pos_start of exon - phase)  mod3
    return False
  elif strand=='-' and mod3(-exon1[1]-exon1[2]) != mod3(-exon2[1]-exon2[2]): #just reversing coordinates. mod3 of negative numbers must work fine e.g. mod3(-1) = 2;  mod3(-8) = 1
    return False
  return True


def mod3(number):
  #to better implement the mod3 with high numbers        STILL TO FINISH
  return number%3  

def exon_length(exon):
  """exon here being just [start, stop]"""
  return ( exon[1]-exon[0]+1 )

def summary_overlap(ov, g1, g2):
  """Prints a summary of the overlap information contained in the object returned by g1.overlaps_with(g2)"""
  if not ov:    return "No overlapping"
  o=''
  for e_index in range(len(g1.exons)):
    o+='1E'+str(e_index+1)+' '+str(g1.exons[e_index][0])+'\t'+str(g1.exons[e_index][1])+' P'+str(g1.phase_of(e_index))+' -- '
    if ov[e_index]:
      index_list=ov[e_index]
      if type(index_list[0])==list:        index_list=[a[0] for a in index_list]
      for e_index2 in index_list:        o+='2E'+str(e_index2+1)+' '+str(g2.exons[e_index2][0])+'\t'+str(g2.exons[e_index2][1])+' P'+str(g2.phase_of(e_index2))+' ++ '
      o=o[:-4]
      if len(index_list)==1 and g2.exons[e_index2]==g1.exons[e_index]:        o+=' ***'
    else:
      o+='None'
    o+='\n'
  return o[:-1]

sequence_db={}
def load_sequence_db(fasta_file, title_function=None):
  """Load a fasta file in memory for being used with function fast_sequence  """
  global sequence_db
  if title_function is None: title_function=lambda x:x.split()[0]
  for title, seq in parse_fasta(fasta_file):
    title_id= title_function(title)
    sequence_db [title_id] = seq      

def get_sequence_db(title=None):
  if title is None: return sequence_db
  else:             
    try:    return sequence_db[title]
    except KeyError:  raise Exception, "MMlib - sequence_db ERROR cannot find a sequence with this title in memory {0}".format(title)
  


class gene(object):
  """This class handles any kind of genomic coordinate, with each object that can consist by more than one range. The names of attributes/methods are thought for protein coding genes but the class can be used for a variety of other purposes and concepts. Typically, the genomic coordinates refer to a target file which is multifasta.
  Important: genomic coordinates are specified 1-based, and all indexes are 1-based unless otherwise specified.

+ Fundamental attributes (more can be added for specific purposes):
  .strand       string    + or -
  .chromosome   string    name of the first word of the fasta header where this gene resides. Can include gi codes. 
  .exons        list      contains elements of the type [start, stop] for each exon. Start is always < stop, but the order of these elements is always upstream-to-downstream, so if the gene is on the negative strand,  .exons[i][0] > .exons[i+1][0]
  .id           string    identifier of the gene. Typically useful when you load a lot at once
  .phases       list      automatically filled. Make sense only for protein coding genes and is used basically only when you compute the overlap between some of them. It's a list of: the number of bases of the first codon of a exon which are on the previous exon. there's a phase element for each exon. Careful: this is not flushed if you modify a gene object.
  .target       string    path to the target multifasta file. Necessary for the .fasta_sequence() method to be called

+ Most useful methods, divided per category: (for a full list, type in python: import MMlib; import inspect; print inspect.getmembers(gene, predicate=inspect.ismethod) )
 # load / create / modify gene objects
  .__init__()     you can create instances in the standard python fashion, calling  gene(). Optionally, you can provide keyword arguments that will be interpreted as attributes of the newly created object. Example:     gene(strand='+', chromosome='chrX', target='/data/genome.fa', id='gene1'). 
  .add_exon(start,stop)     with start< stop, add a new exon to this gene. It will be placed in the right place inside the .exons list, which is determined by the .strand attribute
  .load_gff(file)     to be used when you have a gff file with only a single gene inside. Take care of providing the right tag (default is 'cds'; '*' can be used to consider all gff lines)
  .load_from_header(header) loads a header line, in the format in which is output by the .header() method
  .add_exons_from_positions_summary(p_sum)  loads the .exons attribute from a string like the one returned by the .positions_summary method
  .remove_exon(indices*)    remove one or more exons given their 0-based indices (0 is the first exon, and so on)
  .copy()         returns a full, deep copy of this gene and its attributes. 

 # access information / outputing
  .summary()           returns a string with a summary of the information in the gene
  .boundaries()        returns [start, stop] where these are the two boundaries of the gene
  .length()            returns the total length of all exons in the gene
  .span()              returns the total length of the gene (include introns)
  .header()            returns a header that can be used in a fasta, which includes all the information about this gene. See full __doc__ for details
  .gff()               returns a gff like text. Optional keyarguments can be used to control what it contains  -- see full __doc__ for details
  .positions_summary() returns a string with a summary of the positions of the exons in this gene object. Example:  123-145,265-300,454-565,666-678,988-999
  .fasta_sequence()    returns a tuple [title, seq] where seq is the nt sequence as extracted from the .target file, taking only the entry corresponding to .chromosome (and the right positions)
 
 # derive gene objects from self / or combine with other genes
  .get_exon(exon_index)   returns a copy of this gene, but with a single exon: the one specified by the (0-based!) exon_index argument
  .introns()              returns a copy of self, but with positions corresponding to the introns of self instead. 
  .subseq(start_subseq, length_subseq)    returns a copy of this gene but only with a subseq of it. This may contain less exons than the .self (start_subseq is 1-based)
  .downstream(distance, length)           returns a gene object that is downstream of self. Distance 0 means just downstream (no overlap). .upstream(d, l) exists as well
  .extend(left, right)                    returns a copy of this gene with extended boundaries
  .check_boundaries(chromosome_length)    this function fixes invalid positions that may result from modifying gene object, because pos<1 or pos>chromosome_length
  .is_downstream_of(other_g)              returns True if self is downstream of the other_g (note that .chromosome is not checked). .is_upstream_of(g) exists as well
  .intersection_with(other_g)        returns a gene object with only the intersection with an other_g; the two genes should overlap or a empty gene will be returned
  .union_with(other_g)               returns a gene object with the union with an other_g; the two genes should be on the same chromosome and strand or a empty gene will be returned.
  .subtracted_of(other_g)           returns a copy of self removing all positions in common with other_g
  .overlaps_with(other_g, full=False, phase=True, summary=False, strand=True)    check if this gene overlaps the other_g object. See the full .__doc__ of this method for details
  .restore_absolute_coordinates(parent_gene)   useful when you have coordinates relative to a subseq of a chromosome. This will restore the coordinates to absolute, given the original subseq range
  .relative_coordinates(...)   This is like the reverse of the last function. You have a certain (self) gene on a chromosome. You want to cut the chromosome in a subseq and know the new coordinates of the gene in the obtained subseq.

## other MMlib functions related to the gene class (not methods of gene though):
load_all_genes(gff_file)     returns a list of gene objects corresponding to the entries in a gff files. Define the wanted tag or any other .load_gff argument. This function determines a uniq id for each gff line which decides which exons go together. This defaults to taking the first word of the last gff field but the function can be defined with get_id
merge_genes(gene_list)       merge a list of genes to remove redundancy. It has lots of options, see its __doc__ for details
  """

  def __init__(self, load_file='', species=None, tag=None, strand=None, chromosome=None, target=None, seed=None,  **other_features):
    """ standard python fashion to create gene intances. Optionally, you can provide keyword arguments that will be interpreted as attributes of the newly created object. Example:     gene(strand='+', chromosome='chrX', target='/data/genome.fa', id='gene1'). These can be also non-standard attributes that will be created. 
    The load_file argument can be specified to load the gene from a gff file (using load_gff). The tag argument will be passed to this function. A tag is a gff element class, for example "cds", "exon" and so on... "cds" are loaded by default. "*" can be provided to load all lines. 
    The seed argument can be provided to load another gene object into this instance. It is pretty much like copying, but changing the class to the one of self.
    """
    self.exons=[] #first-coding exon first (if negative strand, position-based reverse order); each exon is [start, stop] where always tart<stop 
    self.species=species
    self.target= target
    self.chromosome=chromosome
    self.strand=strand
    self.phases=[]
    self.id=None
    self.sec_pos=[]
    if load_file and tag:      self.load(load_file, tag)
    elif load_file:            self.load(load_file)  
    for k in other_features:   self[k]=other_features[k]
    if not self.id: self.id= str(uniq_id(self))
    if not seed is None: 
      g=seed.copy(newclass=self.__class__)
      self.__dict__= g.__dict__
      return 

  def __getitem__(self, index):
    if index in self.__dict__:       return self.__dict__[index]
    if type(index)==int:             return self.exons[index]
    else:                            return None
  def __setitem__(self, key, value):     self.__dict__[key]=value
  def __str__(self):                     return self.summary()
  def __nonzero__(self):                 return bool(self.exons)
  def get(self, attribute):
    """ shortcut to get an attribute.. if not present, returns None (does not raise exception).
    """
    return self.__getitem__(attribute)
  def __iter__(self):                   return self.exons
  def __delitem__(self, attribute):
    """Utility meant to delete non std attributes"""
    del self.__dict__[attribute]
  def copy(self, newclass=None):
    """All object features will be copied. if no new_class is provided, the one of the parent (self) will be used """
    if newclass is None:    a=self.__class__()
    else:           a=newclass()
    for i in self.__dict__:
      if i !='exons':      a.__dict__[i]=deepcopy(self.__dict__[i])
    for i in self.exons:      a.exons.append(list(i)) #hard copying the list of exons
    return a

  def summary(self, description='GENE OBJECT', other_fields=[], print_exons=True):
    """ returns a summary string for this object. A description="" argument can be provided, it will be displayed as a title. The other_fields argument can be provided with a list of strings corresponding to a list of attributes that you want to printed as well.  The argument print_exons specifies if the list of exons has to be printed for the object (by default, it is)  """
    o=''
    o+='### '+description+' ###'+'\n'
    o+='# ID:\t'+str(self.id)+'\n'
    if self.species: o+='# SPECIES:\t'+str(self.species)+'\n'
    if self.target:  o+='# TARGET:\t'+str(self.target)+'\n'
    o+='# CHROMOSOME:\t'+str(self.chromosome)+'\n'
    o+='# STRAND:\t'+str(self.strand)+'\n'
    o+='# BOUNDARIES:\t'+str(self.boundaries())+'\n'
    for f in other_fields:         o+="# "+f.upper()+':\t'+ str(self[f])+'\n'
    if self.sec_pos:
      o+='# SEC_POS:\t'
      for i in self.sec_pos:        o+=str(i)+' '
      o+='\n'
    for e_index in range(len(self.exons)):      o+='EXON'+str(e_index+1)+'\t'+str(self.exons[e_index][0])+'\t'+str(self.exons[e_index][1])+'\tPHASE:'+str(self.phase_of(e_index))+'\n'
    return o[:-1]

  def boundaries(self):
    """ Returns [start, stop] where these are the two boundaries of the gene (basically taking first and last exons)"""
    if self.exons:      return [ min(self.exons[0][0], self.exons[-1][0]), max(self.exons[0][1], self.exons[-1][1])] #works both with + and - strand
    else:               return None

  def boundaries_gene(self, minimal=False):
    """ Return a copy of the self gene with a single exon that span its boundaries"""
    if minimal: 
      g=gene(chromosome=self.chromosome, strand=self.strand)
      g.exons=list(self.exons)
    else:       g=self.copy()
    bounds=g.boundaries()
    while g.exons:      g.remove_exon(0)
    g.add_exon(bounds[0], bounds[1])
    return g    

  def length(self):
    """ Returns the total nucleotide length of this gene"""
    l=0
    for st, end in self.exons:       l+=(end-st)+1
    return l
  def __len__(self): return self.length()  
  def span(self):
    """ Returns the total length in nucleotide that this gene spans. It differs from .length() because .span() counts the intron lengths as well"""
    if not self.exons:       return 0
    b=self.boundaries()  
    return b[1]-b[0]+1
    
  def phase_of(self, exon_index):
    """the number of bases of the first codon of a exon which are on the previous exon"""
    ######################### how phases changes from exon to exon
    #phase1#len%3  #phase2#
    #########################
    #  0  #  0  #  0  #
    #  0  #  1  #  1  #
    #  0  #  2  #  2  #
    #########################
    #  2  #  0  #  2  #
    #  2  #  1  #  0  #
    #  2  #  2  #  1  #
    #########################
    #  1  #  0  #  1  #
    #  1  #  1  #  2  #
    #  1  #  2  #  0  #
    #########################
    if len(self.phases)!=len(self.exons):      self.phases=[]
    if not self.phases:
      for e_index in range(len(self.exons)):
        if e_index==0:          self.phases.append(0)
        else:                   self.phases.append(    mod3( self.phases[e_index-1]+ mod3(exon_length(self.exons[e_index-1])) )    )
    return self.phases[exon_index]
    
  
  def load_gff(self, gff_file, tag='cds', check_sec_pos=True, keep_program=False, parse_keywords=False, keep_tag=False, process_line=None):
    """can accept as input direclty text, or path to a file, or a file handler; tag=* means any tag
    If check_sec_pos==True (default), then "Sec_position:" (as present in selenoprofiles gffs) are read from the gff input and this information is stored into the .sec_pos attribute.
    If keep_program==True (not default), a .program attribute is used to keep the program field of the gff being loaded (2nd field)
    If parse_keywords==True (not default), the 9th field of the gff is searched for expressions like:   SOMEKEY:VALUE ; for each one found, a .SOMEKEY attribute is created and assigned to VALUE, taking care of converting VALUE to the appropriate type (string, number, float). A .keywords attribute is created to keep the list of attribute names added this way. You may specify another separator instead of ":" as argument of parse_keywords, as long as just one will be found per text block.
    For maximum flexibility, you have the process_line keyword. You can provide a function that accepts two arguments: a list of the fields (the result of line.split('\t' on the current line)  and  the self gene object being filled. In this way you can keep any other attribute in the line you're interested into.
    """
    if parse_keywords==True: parse_keywords=':'
    if type(gff_file)==str and os.path.isfile(gff_file):
      gff_lines=open(gff_file, 'r').readlines()
      self.id=gff_file
    elif type(gff_file)==file:
      gff_lines=gff_file.readlines()
      self.id=gff_file.name
    else: #string
      gff_lines=gff_file.split('\n')
    try:
      line=gff_lines[0]
      self.chromosome=line.split('\t')[0]
      self.strand=line.split('\t')[6]
      if not self.species:
        path_splt=os.path.abspath(gff_file).split('/')  
        if 'output' in path_splt and len(gff_file.split('/')[-1].split('.'))>=4:
          #selenoprofiles gff
          self.species=path_splt[path_splt.index('output')-1]
      for line in gff_lines:
        if line:
          splt=line.split('\t')
          if tag=='*' or  lower(splt[2]) == lower(tag):
            self.add_exon( int(splt[3]), int(splt[4]) )
            if check_sec_pos and len(splt)>8 and 'Sec_position:' in splt[8]:
              for i in splt[8].split('Sec_position:')[1:]:
                self.sec_pos.append(int(i.split()[0]))
            if parse_keywords and len(splt)>8:
              commentsplt=splt[8].split()
              for word_index, word in enumerate(commentsplt):
                if word.count( parse_keywords  )==1:
                  keyword=word.split(parse_keywords)[0]
                  if word.endswith(parse_keywords):  value=option_value(commentsplt[word_index+1])
                  else:                   value=option_value(word.split(parse_keywords)[1])
                  self[keyword]=value
                  if not hasattr(self, 'keywords'): self['keywords']=[]
                  self['keywords'].append(keyword)

            if keep_program:                 self['program']=splt[1]
            if keep_tag:                     self['tag']=splt[2]
            if not process_line is None:     process_line(splt, self)

    except Exception, e:                     raise Exception, "ERROR loading gff: "+str(gff_file)+' '+str(e)
  load=load_gff
  
  def add_exon(self, start, stop):
    """ usage: gene_obj.add_exon(start, stop) where start<stop. exons can be added in any order, they will be placed in order.
    """
    if stop<start:
      raise Exception, "gene->add_exon ERROR stop must be > than start. called with arguments: "+str(start)+' '+str(stop)
    index_to_append=0
    if self.strand=='+':
      while index_to_append<len(self.exons) and self.exons[index_to_append][0]<start:
        index_to_append+=1
    if self.strand=='-':
      while index_to_append<len(self.exons) and self.exons[index_to_append][0]>start:
        index_to_append+=1
    self.exons.insert(index_to_append, [start, stop])

  def remove_exon(self, *indices):
    """ remove one or more an exon given the index(es). After that, indices will change. can take as arguments both indices or list of indices. NOTE: indices are 0-based here!
    """
    index_list=[]
    for i in indices:
      if i or i==0:
        if type(i)==list:
          index_list.extend(i)
        else:
          index_list.append(i)
    index_list.sort(reverse=True)
    for i in index_list: #reversing so the change of index does not affect this function
      self.exons.pop(i)

  def get_exon(self, exon_index, minimal=False):
    """ Return a copy of the self gene, with a single exon given by the exon index (0 is for the first one); id put by default is the id of the self gene plus the flag "_exonN" where N is the exon index  """
    if minimal:     g=gene(chromosome=self.chromosome, strand=self.strand)
    else:           g=self.copy()
    g.id=self.id+"_exon"+str(exon_index)
    g.exons=[]
    g.add_exon( self.exons[exon_index][0], self.exons[exon_index][1]   )
    return g


  def overlaps_with(self, other_g, full=False, phase=True, summary=False, strand=True):
    """Main function. Computes the overlap between two gene objects. By default phase is checked as well, use kwarg phase=False to disable it (when not dealing with cds). Same for strand.
       It returns False if the chromosome or strand are different, or one of the two objects is empty, or if the gene boundaries are not overlapping, or if they are but no exon is overlapping and same phase of another one in the second gene.
       If there is overlap, it returns: (if Full=False, by default it is)
           a list of lists, one for each exon of the self gene. Each of these contains the indexes of the exons of the second gene overlapping even partially with this exon of the self gene.
           e.g.     [[0], [2], [], [5]] --> G1exon1 overlaps with  G2exon1, G1exon2 overlaps with G2exon3, G1exon3 is not overlapping with anything, G1exon4 is overlapping with G2exon6

      If Full=True,
          instead of just an index, each of these contains a three element list for each exon of the second gene overlapping even partially with this exon of the self gene:  [index, start, stop], where start and stop are the positions of the overlapping intersection.

      if summary=True, instead of normal (or full) output, a tuple (output, summary) is returned, summary being computed by the function summary_overlap()

    """
    if self.chromosome != other_g.chromosome or (strand and self.strand != other_g.strand) or not self.exons or not other_g.exons: #or self.species != other_g.species   ###assuming nobody is so stupid to look for overlap in  genes in different species
      o= False
    #producing false exons (with phase 0) to check boundaries
    else:
      if strand:
        strand=self.strand
      if are_overlapping_exons(self.boundaries(), other_g.boundaries(), strand, False):
          #they are sharing some genomic space. I assume that now testing brutally each exon of gene1 against each exon of gene2 would not be too expensive:
          overlap_list=[]
          for e1_index in range(len(self.exons)):    
            exon_g1=list(self.exons[e1_index])
            exon_g1.append(self.phase_of(e1_index))

            overlapping_this_exon=[]
            for e2_index in range(len(other_g.exons)):
              exon_g2=list(other_g.exons[e2_index])
              exon_g2.append(other_g.phase_of(e2_index))

              if are_overlapping_exons( exon_g1, exon_g2, self.strand, phase   ):
                if full:
                  a=[e2_index]
                  a.extend(intersection_of(exon_g1, exon_g2))
                  overlapping_this_exon.append( a   )
                else:
                  overlapping_this_exon.append(e2_index)

            overlap_list.append(overlapping_this_exon)
          if not any(overlap_list):
            o= False
          else:
            o= overlap_list
      else:
        o= False
    if summary:
      return o, summary_overlap(o, self, other_g)
    else:
      return o

  def subseq(self, start_subseq, length_subseq=None, minimal=False, newclass=None):
    """ This function returns a gene object obtained cutting a subseq of the self object. If the a single exon do not contain the whole length, the returned object will have more than one exon.  
        NB start_subseq is 1 based.    """
    if not minimal: out_gene=self.copy(newclass=newclass); out_gene.exons=[]
    else:           out_gene=gene(chromosome=self.chromosome, strand=self.strand, target=self.target)
    if length_subseq is None: length_subseq = self.length() - start_subseq + 1 

    first_exon_index, length_of_previous_exons = 0, 0
    while first_exon_index<len(self.exons) and self[first_exon_index][1] - self[first_exon_index][0]+1 +length_of_previous_exons   <=  start_subseq-1 : #exiting this loop, the start of secis is included in this exon
      length_of_previous_exons+=  self[first_exon_index][1]-self[first_exon_index][0]+1
      first_exon_index+=1
    if first_exon_index>=len(self.exons):      raise Exception, "ERROR subseq function called with start too high for this gene object: "+gene.summary(self)
    
    #print "first_exon_index:"+str(first_exon_index)
    if self.strand=='+':  
      position_start_of_subseq = self[first_exon_index][0] + (start_subseq-1-length_of_previous_exons)    
      out_gene.add_exon( position_start_of_subseq ,           min(    position_start_of_subseq + length_subseq -1 , self[first_exon_index][1]   )           )
      current_exon = first_exon_index
      remaining_length = length_subseq - (out_gene[0][1]-out_gene[0][0]+1) 
      while remaining_length :
        #print "remaining_length"+str(remaining_length)
        current_exon+=1
        if current_exon>=len(self.exons):          raise Exception, "ERROR subseq function called with length too high for this gene object."
        out_gene.add_exon(  self[current_exon][0]  ,        min(  self[current_exon][0]+remaining_length-1 ,  self[current_exon][1]  )        )
        remaining_length -= (out_gene[-1][1]-out_gene[-1][0]+1)

    elif self.strand=='-':
      position_start_of_subseq = self[first_exon_index][1] - (start_subseq-1-length_of_previous_exons)    
      out_gene.add_exon(  max(    position_start_of_subseq - (length_subseq -1) , self[first_exon_index][0]   )    , position_start_of_subseq  )
      current_exon = first_exon_index
      remaining_length = length_subseq - (out_gene[0][1]-out_gene[0][0]+1) 
      while remaining_length :
        #print "remaining_length"+str(remaining_length)
        current_exon+=1
        if current_exon>=len(self.exons):
          raise Exception, "ERROR subseq function called with length too high for this gene object."
        out_gene.add_exon(    max(  self[current_exon][1]-(remaining_length-1) ,  self[current_exon][0]  )    , self[current_exon][1]     )
        remaining_length -= (out_gene[-1][1]-out_gene[-1][0]+1)      

    return out_gene

  def reverse_subseq(self, start_subseq, length_subseq):
    """This function is equivalent to subseq, except for the fact that the index refers to the end to the gene. so, a start_subseq of 1 and length_subseq of 3 would return a gene object with only the last codon """
    a=self.copy()
    if    self.strand=='+':    a.strand='-'
    elif  self.strand=='-':    a.strand='+'
    a.reorder_exons()
    b=a.subseq(start_subseq, length_subseq)
    b.strand=self.strand 
    b.reorder_exons()
    return b

  def reorder_exons(self):
    """ This function put the exons uin the correct order, that is to say, upstream to downstream, which depends on the .strand attribute."""
    if    self.strand=='+':   self.exons.sort(key=lambda x:x[0])
    elif  self.strand=='-':   self.exons.sort(key=lambda x:x[0], reverse=True)
    

  def restore_absolute_coordinates(self, parent_gene, inplace=True):
    """ This function is useful when the coordinates on the self gene are relative to some portion of the chromosome, they are not absolute. 
    By providing a gene object representing the portion to which they were relative, coordinates are set back to absolute.
    NB: the strand of the result is computed as in maths by the product of the strand of the self and of the parent. 
    If inplace==True, the coordinates are modified in the self object, otherwise a new object with absolute coordinates is returned
    """
    
    g=gene(chromosome=parent_gene.chromosome)
    if self.strand!=parent_gene.strand:
      g.strand='-'
    else:
      g.strand='+'
    for start, stop in self.exons:
      s=parent_gene.subseq(start, (stop-start+1))
      for s_start, s_stop in s.exons:
        g.add_exon(s_start, s_stop)
        
    if inplace:
      self.exons=g.exons
      self.chromosome=g.chromosome
      self.strand=g.strand
    else:
      return g

  def relative_coordinates(self, subseq_start, subseq_end, reverse=False, chromosome='subseq_sequence'):
    """ This function can be seen as the reverse of restore_absolute_coordinates.
    You have a certain (self) gene on a chromosome. You want to cut the chromosome in a subseq and know the new coordinates of the gene in the obtained subseq.
    If you also reverse complemented the subseq, used reverse=True.  To define the chromosome attribute of the returned gene, use the chromosome=... keyarg    
    All coordinates are 1 based
    """
    new_g=gene(chromosome=chromosome)
    if not reverse: 
      #subseq was cut normally, no rev comp   
      new_g.strand = self.strand
      for start, end in self.exons:
        new_g.add_exon(    start - subseq_start +1 ,  end-subseq_start+1      )
    else:
      #subseq was cut and then rev comp   
      if self.strand=='-':   new_g.strand='+'
      else:                  new_g.strand='-'      
      for start, end in self.exons:
        new_g.add_exon(    subseq_end - end +1,   subseq_end - start +1    )
    return new_g      

  def downstream(self, distance, region_length, **keyargs):
    """return a gene object with the region downstream this gene. vars: distance, region_length. distance=0 means the region is right downstream - starting from the very next base. Using a negative distance value, you will obtain a partial overlap with this gene (e.g. distance=-3 and region length=6 will return a gene object with the last codon of this gene plus the next codon)
    You can provide as keyargs any other attribute that you want to set in the resulting gene object. If you don't set the id, it will be set to its uniq_id in the python environment.
    """
    g=gene(strand=self.strand, chromosome=self.chromosome, target=self.target)
    if self.exons:
      if self.strand=='+':          g.add_exon(self.exons[-1][1]+1+distance, self.exons[-1][1]+1+distance+region_length-1)  
      elif self.strand=='-':        g.add_exon( self.exons[-1][0]-1-distance-region_length+1, self.exons[-1][0]-1-distance)   
    for k in keyargs:      g[k]=keyargs[k]
    if not g.id:      g.id=str(uniq_id(g))
    return g

  def upstream(self, distance, region_length, **keyargs):
    """see downstream function.
    """
    g=gene(strand=self.strand, chromosome=self.chromosome, target=self.target)
    if self.exons:
      if self.strand=='+':
        g.add_exon(self.exons[0][0]-1-distance-region_length+1 ,self.exons[0][0]-1-distance)        
      elif self.strand=='-':
        g.add_exon(self.exons[0][1]+1+distance, self.exons[0][1]+1+distance+region_length-1)        
    for k in keyargs:
      g[k]=keyargs[k]
    if not g.id:
      g.id=str(uniq_id(g))
    return g

  def is_upstream_of(self, other_g, max_overlap=0):
    """checks if gene self is upstream (not a single base overlap, unless max_overlap!=0) of other_g ; if the genes are not on the same strand, an exception is raised.
    It the output is "true", actually te distance between the two genes is reported. If the two genes are exactly adjacent, the distance reported is 1 (to avoid a 0 in a bool evaluation). If max_overlap is not 0, in case this rescues an evaluation , the n of bases in overlap are returned as a negative value. 
    NOTE: the chromosome attribute is not checked
     """
    if self.strand!=other_g.strand:       raise Exception, "ERROR is_upstream_of function: genes are not on the same strand"
    if self.strand=='+':                  dist= other_g.exons[0][0]-self.exons[-1][1]
    elif self.strand=='-':                dist= self.exons[-1][0]-other_g.exons[0][1]

    if dist<=0:      
        if dist+max_overlap>0:   return dist-1
        else:                     return False
    else:            return dist
      
  def is_downstream_of(self, other_g, max_overlap=0):
    """ see is upstream of"""
    return other_g.is_upstream_of(self, max_overlap=max_overlap)

  def extend(self, left=0, right=0, inplace=False, minimal=False, down=None, up=None):
    """ This function returns a copy of the self gene object with extended boundaries. The new coordinates could be out of the allowed range... it is suggested to run check_boundaries right after. If inplace==True, the gene is modified inplace and not returned. If down or up are specified, the direction of the extension depends on the strand
"""
    if not down is None or not up is None:
      if down is None: down=0
      if up   is None: up=0
      if   self.strand=='+':  return self.extend(left=up, right=down, inplace=inplace, minimal=minimal)
      elif self.strand=='-':  return self.extend(left=down, right=up, inplace=inplace, minimal=minimal)

    if not self.exons: raise Exception, "gene->extend ERROR can't extend an empty gene! no exons were found for "+str(self.id)
    if inplace:     g=self 
    elif  minimal:  
      g=gene(chromosome=self.chromosome, strand=self.strand, target=self.target)
      g.exons=list(self.exons)
    else:           g=self.copy()
    if g.strand=='+':
      g.exons[0] = [g.exons[0][0]-left,  g.exons[0][1]]
      g.exons[-1]= [g.exons[-1][0],      g.exons[-1][1]+right]
    elif g.strand=='-':
      g.exons[0]=   [g.exons[0][0],        g.exons[0][1]+right]
      g.exons[-1] = [g.exons[-1][0]-left,  g.exons[-1][1]]
    if not inplace:         return g

  def check_boundaries(self, chromosome_length=False):
    """ This function checks if the gene object has some illegal exons, that is to say, exons that spans regions which don't exist, either with positions <1 or > than the length of the chromosome.
    This must be provided as an argument, otherwise this fact is not checked.
    Illegal exons are removed when they are totally illegal, or set to legal positions when at least a small region of it is legal. A string message with the removed exon is returned if any was removed.
    Returns 0 if nothing has changed, 1 if some positions were changed, a string with details if some entire exons were removed
    """    
    modified=0
    to_remove={} #keep the indices of the exons which are out of bound
    for exon_index in range(len(self.exons)):
      start, stop= self.exons[exon_index]
      if stop<1:        to_remove[exon_index]=1
      if start<1:
        self.exons[exon_index][0]=1
        modified=1
      if chromosome_length:
        if start>chromosome_length:
          to_remove[exon_index]=1          
        if stop>chromosome_length:
          self.exons[exon_index][1]= chromosome_length
          modified=1

    if to_remove:
      msg='removed exons: '
      for k in sorted(to_remove.keys()):
        msg+=str(k)+' ('+str(self.exons[k][0])+','+str(self.exons[k][1])+') '
      self.remove_exon( to_remove.keys() ) 
      return msg
    return modified

  
  def intersection_with(self, geneobj, check_phase=False):
    """This method returns another gene object, with only the intersection of the two gene object in input. The phase is considered only if check_phase is True"""
    out_gene=self.copy(); out_gene.exons=[]; out_gene.chromosome=None; out_gene.strand=None
    if not self.chromosome==geneobj.chromosome:
      return out_gene
    out_gene.chromosome=self.chromosome

    if not self.strand==geneobj.strand:
      return out_gene
    out_gene.strand=self.strand

    ov=self.overlaps_with(geneobj, full=True, phase=check_phase)
    if ov:
      for selfindex in range(len(ov)):
        if ov[selfindex]:
          for geneindex, start, stop in ov[selfindex]:
            out_gene.add_exon(start, stop)
    return out_gene  


  def union_with(self, geneobj, check_phase=False, id='SUM', **keyargs):
    """This method returns another gene object, with the union of the two gene objects in input. The phase is considered to check overlap, only if check_phase is True.
The variable id defines what the id of the result will be like. If it is "SUM", the id will be the results of joining the two original ids with "_+_".
If id is "LONGEST", the resulting gene object will have the id of the longest of the joined gene object only (or the first one if length is the same). 
Lastly, if you defined an id which is not SUM or LONGEST, it becomes the id of resulting gene object.

  Also, this can be checked for all other key arguments given. (So you can define also program='SUM' or program='some_program')
    
    """
    #print self.summary('BASE1')
    #print geneobj.summary('BASE2')
    
    out_gene=self.copy(); out_gene.exons=[]
    if not self.chromosome==geneobj.chromosome:
      return out_gene
    out_gene.chromosome=self.chromosome
    if not self.strand==geneobj.strand:
      return out_gene
    out_gene.strand=self.strand
    
    for start, end in self.exons + geneobj.exons:
      #print "cycling exon", start, end
      t_gene=gene(strand=self.strand, chromosome=self.chromosome)
      t_gene.add_exon(start, end)
      o=t_gene.overlaps_with(out_gene, phase=check_phase)
      if o: #the exon which I am about to add overlaps with something that I already added. Therefore, I must instead replace everything that overlaps with this with a new exon that includes them all
        exons_to_replace = o[0] #list of indices
        start, end = min( [out_gene[i][0] for i in o[0]]+[t_gene[0][0]] ) , max ([out_gene[i][1] for i in o[0]]+[t_gene[0][1]] ) 
        #print "removed existing: ", o[0], '; adding: ', start, end
        out_gene.remove_exon(o[0]) #removing all overlapping exons in out_gene
      out_gene.add_exon(start, end)

    if not keyargs.has_key('id'):   keyargs['id']=id
    for k in keyargs:
      if keyargs[k]=='SUM':
       out_gene[k]= str(self[k])+'_+_'+str(geneobj[k])
      elif keyargs[k]=='LONGEST':
        if self.length()>=geneobj.length():
          out_gene[k]= str(self[k])
        else:
          out_gene[k]= str(geneobj[k])
      else:
          out_gene[k]= keyargs[k]

    #print out_gene.summary('after union')
    return out_gene  

  
  def subtracted_of(self, geneobj):
    """This returns a gene object, result of the self object subtracted of the second gene object. All overlapping regions are removed. """
    out_gene=self.copy(); out_gene.exons=[]
    if not self.chromosome==geneobj.chromosome:       return self
    out_gene.chromosome=self.chromosome

    if not self.strand==geneobj.strand:
      return self
    out_gene.strand=self.strand

    ov=self.overlaps_with(geneobj, full=True, phase=False)
    if not ov:      return self
    else:
      for selfindex in range(len(ov)):
        if ov[selfindex]:
          for geneindex, start, stop in ov[selfindex]:
    
            if not self[selfindex][0]==start:
              out_gene.add_exon(self[selfindex][0], start-1)

            if not self[selfindex][1]==stop:
              out_gene.add_exon(stop+1, self[selfindex][1])

        else:
          out_gene.add_exon( self[selfindex][0],  self[selfindex][1]   )  
      
    return out_gene  

  def introns(self, minimal=False, skip_null_introns=False):
    """ Returns a gene object containing the introns of the self object. The self is copied so that all features (apart the genomic coordinates of exons) will be present in the returned gene.
    The id assigned will be equal to the id of self gene plus the flag "_introns". NOTE: if the self object has exons which are adjacent, an exception is raise; instead, it is tolerated (and these introns will just not be returned)  if skip_null_introns=True
     """
    if minimal:      g=gene(chromosome=self.chromosome, strand=self.strand)
    else:            g=self.copy()
    g.id=self.id+"_introns"
    g.exons=[]
    for exon_index in range(1, len(self.exons)):
      if self.strand=='+':
        intron_start =      self.exons[exon_index-1][1]+1
        intron_end   =      self.exons[exon_index]  [0]-1
      elif self.strand=='-':
        intron_start =      self.exons[exon_index]  [1]+1
        intron_end   =      self.exons[exon_index-1][0]-1
      if intron_start == intron_end+1: 
        ## adjacent exons! 
        if skip_null_introns: continue
        raise Exception, "gene->introns  ERROR there are adjacent exons! cannot return a null intron!! Use skip_null_introns=True to tolerate this and skip these introns"
      g.add_exon(intron_start, intron_end)
    return g          
       
  def bed(self, exon_index=None, show_id=True, strand=True, score=False, other_fields=[], sep='\t'):
    """Generic function to produce a bed from a gene object. E.g. default: (spacers are tab):
chr4  9236903   9236953   BE_0001  0   +
chr4  49110668  49110718  BE_0001  0   +
By default (exon_index==None) this includes one line for each (start,end) in self.exons. If you provide a 0-based exon_index, it gives  a single line for that certain exon.
By default (show_id==True) the name in the bed out (fourth field) is printed and it is self.id; provide a string argument to show_id to print this as name instead. Provide show_id=False to show_id to show minimal bed (e.g. "chr4  9236903   9236953"). In this case, all following args are ignored.
  If strand is False, then only the first 4 fields in bed are printed, (minimal bed + name/id).    
  If strand is True (default), at least the first 6 fields in bed are printed (up to strand including score), with score set to 0 unless provided as argument score.
    If strand==True, then other_fields can be provided, so each item in this list is converted to string and added at every line.
By default output is tab-separated; use sep==' ' to produced space separated
"""
    out=''
    if exon_index is None: considered_exons=self.exons
    else:                  
      try:   considered_exons=self.exons[exon_index]
      except IndexError: raise Exception,"gene->bed ERROR exon index is invalid! {0} Gene: {1}".format(exon_index, self.header())
    for start, end in considered_exons: 
      out+='\n{1}{0}{2}{0}{3}'.format(sep, self.chromosome, start, end)
      if not show_id:        break
      elif show_id==True:    out+='{0}{1}'.format(sep, self.id)
      else:                  out+='{0}{1}'.format(sep, show_id)        
      if strand and not score: score=True
      if not score:         break
      else:                 
        if score==True: score=0
        out+='{0}{1}'.format(sep, score)      
      if not strand:         break
      else:                  out+='{0}{1}'.format(sep, self.strand)

      for the_field in other_fields: out+='{0}{1}'.format(sep, the_field)
    return out[bool(out):]
    

  def gff(self, tag='', is_gtf=False, comment='', program='', sec_pos=False, id=None, score=None, position_features=[], last_field=None):
    """Generic function to produce a gff from a gene object. The options control the open fields. sec_pos was used for selenoprofiles < 3.1; position_features is required for >= 3.1. it is thought to add (one or few) single position features as comment in the same line of the exon they belong to. format: [[position, what_to_write], ... ]  
Last field overrides the id and the comment arguments """
    if id is None: its_id=str(self.id)
    else:          its_id=id
    if not tag:
      if not hasattr(self, 'tag'): tag='cds'
      else:                        tag=getattr(self, 'tag')
    if not program and self['program']:      program=self['program']
    elif not program:                        program='generic_program'
    if score is None: score_txt='.'
    elif type(score) in (str, float, int):  score_txt=str(score)

    out=''
    for start, end in self.exons:
      add_sec_pos=''
      #################
      ### deprecated
      if sec_pos:   ##### note: this was added for selenoprofiles, but it is not used anymore
        sec_positions_list=[] 
        if self.sec_pos:        sec_positions_list=self.sec_pos    # self.sec_pos is not used anymore. All the necessary information is already in the aminoacid sequence in the alignment attribute. So, I built a method available for p2ghit gene class that returns a position list as the self.sec_pos list was.
        else:                   
          try: sec_positions_list=self.sec_positions_list()
          except: pass
        for sec_position in sec_positions_list:   ### deprecated
          if sec_position>= start and sec_position<= end:
            add_sec_pos+='Sec_position:'+str(sec_position)+' '
      #################      
      if position_features:
        for pos, what_to_write in position_features:
          if pos>= start and pos<= end:
            add_sec_pos+=what_to_write+' '

      if last_field is None:
        if is_gtf:          last_field_txt='gene_id "'+its_id+ '"; transcript_id "'+its_id+'"; '+add_sec_pos+comment
        else:               last_field_txt=its_id+' '+add_sec_pos+comment
      else:                 last_field_txt=last_field
      out+=self.chromosome+'\t'+program+'\t'+tag+'\t'+str(start)+'\t'+str(end)+'\t'+score_txt+'\t'+self.strand+'\t.\t'+last_field_txt+'\n'
    return out[:-1]

  def geneid_gff(self, program='', id=None, score=None, no_frame=False):
    """ Output a gff to be used with  geneid option -R or -O
    score: None --> . is set ; which means: element is compulsory (like infinite score)
           float --> this value is set for all exons
           list ---> list of values, of exons, in order, upstream to downstream
    no_frame: False --> assume this gene is protein coding, so we derive the information of the frame and we feed it to geneid,
              True  --> if this is not a protein coding gene structure, set frame column to ".", which is equivalent to set it to unknown 
    """
    if type(score)==list and len(score)!= len(self.exons): raise Exception, "gene->geneid_gff ERROR list of scores provided does not have exactly a score for each exon"
    if id is None: its_id=str(self.id)
    else:          its_id=id
    if not program and self['program']:      program=self['program']
    elif not program:                        program='generic_program'
    out=''
    for index, element in enumerate(self.exons):
      start, end = element
      last_field=its_id

      if score is None: score_txt='.'
      elif type(score) in [float, int]:   score_txt=str(round(score, 2))
      elif type(score)==list:             score_txt=str(round(score[index], 2))
      else:          raise Exception, "gene->geneid_gff ERROR score not accepted: "+str(score)
      
      if no_frame:   frame_txt= '.'
      else:          
        phase=  self.phase_of(index)
        frame=  (3-phase)%3
        frame_txt=str(frame)
        
      if len(self.exons)==1: tag='Single'
      ######## finish!
        

      out+=self.chromosome+'\t'+program+'\t'+tag+'\t'+str(start)+'\t'+str(end)+'\t'+score_txt+'\t'+self.strand+'\t'+frame_txt+'\t'+last_field+'\n'
    return out[:-1]


  def extend_orf(self, chromosome_length=None, stops={'TGA':1, 'TAG':1, 'TAA':1}, starts={'ATG':1}, up=True, down=True, get_seq=lambda x:replace(upper(x.fasta_sequence()[1]),'U','T'),  extension_parameter=1000, keep_seq=False):
    """ Extend a gene both upstream and downstream so that the first codon will be a start (the most upstream one without including a stop), and the last codon will be a stop. It returns a tuple like   (nt_extended_upstream, nt_extended_downstream), how much was extended in the two directions.
You must provide chromosome_length, which is the length of the sequence to which self.chromosome points to; this is used to check if during the extension, this boundary is passed. You can avoid providing this only in the case in which the variable chromosome_lengths is defined in MMlib and has self.chromosome as key. 
The allowed Starts and Stops can be provided as argument, in the form of any iterable (checked with "in"; use hashes for max speed). If you want to extend upstream to the next stop regardless of start codon, use starts={}
The arguments up and down (True by default) control the direction in which the extension is attempted (upstream and/or downstream).
The argument get_seq tells how the sequence of the gene object should be retrieved. There are at least two possibilities: 
- lambda x:x.fasta_sequence()[1]     #default; the .target attribute of the gene object must be defined, and the split_folder variable must be available in MMlib (see set_local_folders)
- lambda x:x.fast_sequence()         #faster but more memory intensive; to call this, the full target sequence database must be loaded in memory with function load_sequence_db
The argument extension_parameter tells how big are the chunks of sequence retrieved at once; this parameter is applied to both sides of     
The argument keep_seq, if True, sets the sequence of the new object as its .seq attribute before returning it.
   """
    if chromosome_length is None: 
      try:    chromosome_length=chromosome_lengths[self.chromosome]
      except: raise Exception, "extend_orf ERROR you must provide chromosome_length as argument; this may be skipped if chromosome_lengths is defined in MMlib and has self.chromosome as key. But this chromosome was not found: "  +str(self.chromosome)
    if not hasattr(self, 'original_bounds'): self.original_bounds=self.boundaries()
    allowed_letters='ATGC'
      ### extending and getting sequence. Doing it at once to avoid calling get_seq multiple times.
    big_extended_g= self.extend(right=extension_parameter, left=extension_parameter, inplace=False)
    big_extended_g.check_boundaries(chromosome_length)
    down= down and big_extended_g.downstream(0, 1).boundaries()[0] != self.downstream(0, 1).boundaries()[0] 
    up=   up   and big_extended_g.upstream(0, 1).boundaries()[0]   != self.upstream(0, 1).boundaries()[0] 
    big_extended_seq= upper( get_seq(big_extended_g) )
    if self.strand=='+':        
      offset= self.boundaries()[0]-big_extended_g.boundaries()[0]
      extended_out_downstream=  ( big_extended_g.boundaries()[1] - self.boundaries()[1] ) != extension_parameter
    elif self.strand=='-':      
      offset= big_extended_g.boundaries()[1]-self.boundaries()[1]
      extended_out_downstream=  (  self.boundaries()[0] - big_extended_g.boundaries()[0] ) != extension_parameter
    extended_out_upstream=    offset!= extension_parameter
      #seq_g= big_extended_seq [offset:offset+  self.length() ]

    if up or down:
      last_codon=big_extended_seq [offset+self.length()-3:offset+self.length() ]
      while down and not last_codon in stops:  #  last codon is not stop codon: keep extending. out of boundaries or weird sequence cause it to break
        ## going downstream
        codon_downstream= big_extended_seq [offset+  self.length(): offset+self.length()+3]
        ### forcing to stop if letters are not standard (e.g. Ns)
        if not all ( [c in allowed_letters for c in codon_downstream] ):     down=False; break
        if len(codon_downstream)!=3:                                         break
        self.extend(down=3, inplace=True)
        last_codon=codon_downstream

      ## first we extend up to the stop, then we subseq
      first_codon = big_extended_seq [offset:offset+3]
      while up and not first_codon in stops:
        codon_upstream= big_extended_seq [offset-3:offset]
        if not all ( [c in allowed_letters for c in codon_upstream] ):     up=False; break
        if len(codon_upstream)!=3:                                         break
        self.extend(up=3, inplace=True)
        offset-=3
        first_codon=codon_upstream
      
      rerun_up=False;       rerun_down=False
      if down and not last_codon in stops and not extended_out_downstream  \
            and  abs(big_extended_g.downstream(0, 1).boundaries()[0] - self.downstream(0, 1).boundaries()[0]) <3:         rerun_down=True    ## we extended to the limit
      if up and not first_codon in stops and  not extended_out_upstream    \
            and abs(big_extended_g.upstream(0, 1).boundaries()[0] - self.upstream(0, 1).boundaries()[0]) <3:              
         rerun_up=True    ## we extended to the limit
      if rerun_up or rerun_down:
        #print self
#        print "rerun", rerun_up, rerun_down
        return self.extend_orf(chromosome_length, stops=stops, starts=starts, get_seq=get_seq, up=rerun_up, down=rerun_down, extension_parameter=extension_parameter, keep_seq=keep_seq)

      ### cut to first start
    if starts:
        first_start=None
        for codon_index in range(self.length()/3): 
          codon= big_extended_seq [ offset+codon_index*3:offset+codon_index*3+3 ]
          if codon in starts: first_start=codon_index; break
        if not first_start is None: 
          if   self.strand=='+':       actually_reducing=  self.boundaries()[0]+first_start*3 > self.original_bounds[0] 
          elif self.strand=='-':       actually_reducing=  self.boundaries()[1]-first_start*3 < self.original_bounds[1] 
          if not actually_reducing:    ## cutting here
            new_coords_g=self.subseq( first_start*3 + 1 ) 
            offset+=first_start*3
            self.exons=list(new_coords_g.exons)
          elif up: ## we went up to the closest stop but found no Methionine in this upstream extension. Let's go back to the original 5' boundary 
            if   self.strand=='+':       offset+=self.original_bounds[0]-self.exons[0][0]; self.exons[0][0]=self.original_bounds[0]
            elif self.strand=='-':       offset+=self.exons[0][1]-self.original_bounds[1]; self.exons[0][1]=self.original_bounds[1]
    if keep_seq: self.seq= big_extended_seq [offset:offset+  self.length() ]
    if   self.strand=='+':      extended_up = self.original_bounds[0] - self.boundaries()[0];  extended_down = self.boundaries()[1] - self.original_bounds[1]  
    elif self.strand=='-':      extended_up = self.boundaries()[1] - self.original_bounds[1];  extended_down = self.original_bounds[0] - self.boundaries()[0]
    del self.original_bounds
    return extended_up, extended_down

  def header(self, no_id=False, no_species=False, no_target=False, no_chromosome=False, no_strand=False, compress=False, max_exons_uncompressed=5, **other_fields):
    """ This returns a header with all information about the gene... it can be used as a one line description. 
      As keyed arguments, you can provide attributes of the gene that you want to be included in the header. Example, if the gene X has a .program attribute with value "blast", you can call the function X.header(program=1) and the string " program:blast " will be included in the header, at the end.
      You can also specify function that must be called instead of attributes. In this case just specify a field called like function_xxx , and the function xxx will be called and the string xxx:output will be included in the header. No arguments can be specified for such functions
    """
    species_string=''
    if bool(self.species and not no_species):        
      species_string=str(self.species)
      if ' ' in species_string: species_string='"'+species_string+'"'
      species_string=' species:'+species_string
    out=(str(self.id)+" ")*int(not no_id) + ("chromosome:"+self.chromosome+' ')*int(not no_chromosome)  +("strand:" +self.strand+' ')*int(not no_strand)+'positions:'+self.positions_summary(compress, max_exons_uncompressed) +species_string+ (' target:'+str(self.target))*int(bool(self.target and not no_target))
    
    for k in other_fields:
      if k.startswith('function_'): 
        f=k.split('function_')[1]
        try:      out+=' '+f+':'+str(getattr(self, f)())
        except:   
          printerr("gene->header ERROR searching/executing function: "+str(f)+' in gene with id '+str(self.id), 1)
          raise
      else:
        try:      out+=' '+k+':'+str(self[k])
        except:   raise Exception, "gene->header ERROR can't find field: "+str(k)+' in gene with id '+str(self.id)
    
    return out
    
  def load_from_header(self, header_line):
    """ Reverse of last function... to load a gene back from there."""
    while header_line[0] in ['#', '>']: header_line=header_line[1:] #removing trailing characters so it can be used on fasta headers or comment lines.
    self.__init__()
    if 'chromosome:' in header_line:       self.chromosome=header_line.split('chromosome:')[1].split()[0]
    if 'strand:'     in header_line:       self.strand=header_line.split('strand:')[1].split()[0]
    self.add_exons_from_positions_summary( header_line.split('positions:')[1].split()[0] )
    if 'species:' in  header_line:
      species_string=header_line.split('species:')[1].split()[0]
      if  species_string[0]=='"': species_string=header_line.split('species:"')[1].split('"')[0]
      if  species_string!='None': self.species=species_string
    if 'target:' in  header_line:
      target_string=header_line.split('target:')[1].split()[0]
      if  target_string!='None': self.target=target_string
    first_word= header_line.split()[0]
    if not (first_word.startswith('chromosome:') or first_word.startswith('strand:') or first_word.startswith('positions:') or first_word.startswith('species:') or first_word.startswith('target:')) :       self.id=header_line.split()[0 ]

  def positions_summary(self, compress=False, max_exons_uncompressed=5):
    """ Returns a string with a summary of the positions of the exons in this gene object. 
    Example:  123-145,265-300,454-565,666-678,988-999
    If compress == True, this summary is compressed and limited to certain exons, whose number is decided by: max_exons_uncompressed
    Example compress==True, max_exons_uncompressed=4 --> 123-145,265-300,454-565,...,988-999
    """
    out=''
    if self.exons:
      if compress and len(self.exons) > max_exons_uncompressed:
        for start, stop in self.exons[:max_exons_uncompressed-1]:        out+=str(start)+'-'+str(stop)+','
        out+='...,'; out+=str(self.exons[-1][0])+'-'+str(self.exons[-1][1])+','
      else:
        for start, stop in self.exons:        out+=str(start)+'-'+str(stop)+','
      out=out[:-1]    #removing last ,
    return out

  def add_exons_from_positions_summary(self, pos_summ):
    """ Read an uncompressed  positions summary once produced by the above function and add to the self object all exons in it. If the summary is compressed, an exception is raised """
    for piece in pos_summ.split(','):
      if piece=='...': raise Exception , "gene->add_exons_from_positions_summary ERROR can't load exons from a compressed positions summary!"
      start= int(piece.split('-')[0]) ; end=int(piece.split('-')[1])
      self.add_exon(start, end)

  def fasta_title(self):
    """ Analog to header, but it is used to generate a single word with all the information to load back the gene after, and be able to understand which portion of which chromosome was exactly written in a fasta file.
     This is useful since some programs keep only the first word of the fasta header. """
    return self.chromosome+'[positions:'+self.positions_summary() +'][strand:'+self.strand+']'
    
  def load_from_fasta_title(self, title):
    """ Reverse of last function... to load a gene back from there."""
    self.chromosome  = title.split('[positions:')[0]
    self.strand      = title.split('[strand:')[1].split(']')[0]
    self.add_exons_from_positions_summary(title.split('[positions:')[1].split(']')[0])    

  def fast_sequence(self):
    """ Can be used only if the target file for this gene is loaded in memory in the sequence_db object (see function load_sequence_db). provides a much faster way to get sequences than method fasta_sequence. returns string with the sequence """    
    if not sequence_db.has_key(self.chromosome): raise Exception, "ERROR fast_sequence() cannot find chromosome identifier: "+str(self.chromosome)
    seq_out=''
    if self.strand=='-':
      for start, end in reversed(self.exons):      seq_out+= sequence_db[self.chromosome] [ start-1:end ]      
      seq_out= reverse_complement(seq_out)
    else:
      for start, end in self.exons:      seq_out+= sequence_db[self.chromosome] [ start-1:end ]
    if len(seq_out)!=self.length(): raise Exception, "ERROR fast_sequence() wrong sequence length in memory! aborting "
    return seq_out

  def fasta_sequence(self, to_file='', target='', chromosome_file='', split_exons=False, title=''):
    """ Cut the sequence corresponding to this gene object. If to_file is defined, nothing is returned, and the subseqed sequence is written to "to_file" argument. ; if it is not defined, a single (title, seq) object is returned.
  target can be defined as key arg to override the target attribute of the gene object. It should point to the genome (multifasta or not) file.
  if split_exons is True, then multiple exons will be present in the output as fata entries (so, if split_exons and to_file is not defined, a list of (title, seq) is returned instead of a single entry.
  if you don't define title, this will be computed with the function header. 

  If title is 'fasta_title', the method fasta_title is used to determine the title of the output fasta.
  If split_exons is True, then the word "_EXON" plus the index is added to the FIRST WORD of the title (defined as argument or not).

    """
      
    if not chromosome_file and not target and not self.target :      
      if 'reference_genome_filename' in globals() and is_file(reference_genome_filename):        target=reference_genome_filename
      else:       raise Exception, "gene-> fasta_sequence ERROR the target is not defined "
    elif not chromosome_file and not target:      target=self.target
    if chromosome_file:
      if not is_file(chromosome_file):       raise Exception, "gene-> fasta_sequence ERROR the chromosome_file "+chromosome_file+" was not found"
      chrom_file=chromosome_file
    else:      chrom_file=fastafetch(split_folder, self.chromosome, target)
    if not to_file:      file_out=temp_folder+'gene_fasta_seq'
    else:                file_out=to_file

    bbash('>'+file_out)      

    if title=='fasta_title':             title=self.fasta_title()
    elif not title:                      title=self.header()

    if not split_exons:
      write_to_file(">"+title, file_out)
      for start, stop in self.exons:
        if self.strand=='+':          fastasubseq(chrom_file, start-1, stop-start+1, file_out, pipecommand=" | gawk 'NR>1' >") #appending
        elif self.strand=='-':
          fastasubseq(chrom_file, start-1, stop-start+1, temp_folder+'gene_fasta_seq_to_revcomp') 
          bbash("fastarevcomp "+temp_folder+"gene_fasta_seq_to_revcomp  "+"| gawk 'NR>1' "+" >> "+file_out) 
          
    else:
      for exon_index in range(len(self.exons)):

        start, stop = self.exons[exon_index]
        this_title=title.split()[0]+'_EXON'+str(exon_index+1)  +' '*int(  len(title.split())>1 )+   join(title.split()[1:], ' ')

        bbash('echo ">'+this_title+'" >> '+file_out)    
        if self.strand=='+':          fastasubseq(chrom_file, start-1, stop-start+1, file_out, pipecommand=" | gawk 'NR>1' >") #appending
        elif self.strand=='-':
          fastasubseq(chrom_file, start-1, stop-start+1, temp_folder+'gene_fasta_seq_to_revcomp') 
          bbash("fastarevcomp "+temp_folder+"gene_fasta_seq_to_revcomp   "+"| gawk 'NR>1' "+" >> "+file_out) 

    if not to_file:
      if not split_exons:        return parse_fasta(file_out).all()[0]
      else:                      return parse_fasta(file_out).all()

def intersection_of(range1, range2):
  """assuming they are overlapping, and that they're like [start, stop]  with start < stop"""
  return [ max(range1[0], range2[0]), min(range1[1], range2[1]) ]

uniq_idfunctions_hash={'gtf': lambda s:s.split('transcript_id "')[1].split('"')[0], 'selenoprofiles1_gff': lambda s:s.split('\t')[1], 'selenoprofiles2_gff':lambda s:s.split('\t')[8].split()[0], 'gff':lambda s:s.strip().split('\t')[-1].split()[0]}

def get_gff_format(gff_file):
  gff_format=''
  check_file_presence(gff_file)
  sample_lines=bbash('head -15 '+gff_file)
  if 'transcript_id "' in sample_lines:
      gff_format='gtf'
  elif '\tSP.' in sample_lines:
      gff_format='selenoprofiles1_gff'
  elif 'selenoprofiles' in sample_lines: 
      gff_format='selenoprofiles2_gff'
  elif gff_file.endswith('.gff'):   gff_format='gff'
  if gff_format:    return gff_format
  else:             raise Exception, "ERROR unknown gff format for file: "+gff_file


def load_all_genes(gff_file, tag='cds', get_id='', add=None, is_sorted=False, **load_gff_args):
  """ load and returns all genes from a gff file, determining which line belong to which gene using the function get_id, given as input.
  This function must take a line as input and return its id, which is the same for lines describing the same gene object to be loaded. If not provided, it uses defaults function which depend on the extension of the file loaded
  Argument add can be provided: a function that takes the string of gff lines and the gene object, and may manipulate the gene object reading information from the lines.
  Normally the function makes no assumptions on the distribution of the lines belonging to the same entry (gene) in the file. If the file is sorted (the lines of the same gene entry are consecutive) you can specify is_sorted=True to make the function faster and less memory expensive.
  You can specify additional keyword arguments to control the behaviour of the load_gff function used to load the gff lines. available arguments are (see gene.load_gff documentation):
  check_sec_pos=True, keep_program=False, parse_keywords=False
  """
  #I changed the name of the keyarg in this function. This is to work with old programs
  if 'uniqid_function' in load_gff_args:    
    uniqid_function=load_gff_args['uniqid_function']
    del load_gff_args['uniqid_function']
  else:                                     uniqid_function=get_id  ### default

  if not uniqid_function:    uniqid_function=uniq_idfunctions_hash[get_gff_format(gff_file)]
  genes=[]
  cfile=open(gff_file, 'r')

  if not is_sorted: 
    #more memory expensive, but handles better any gff
    id2lines_list={}
    for line in cfile:        
      line=line.strip()      #;print [line]
      if line and not line[0]=="#" and (tag=='*' or lower(line.split('\t')[2]) == lower(tag) ):
        the_id=uniqid_function(line)
        if not id2lines_list.has_key(the_id): id2lines_list[the_id]=''
        id2lines_list[the_id]+=line +'\n'
    for the_id in sorted( id2lines_list.keys() ):
      gff_lines= id2lines_list[the_id]
      g=gene()
      g.load_gff( gff_lines, tag, **load_gff_args) 
      g.id=the_id
      if not add is None:  add(gff_lines, g)
      genes.append(g)            

  else: 
    ### old code. not proud. and:   it can't work if the lines referring to the same entry are not consecutive
    current_id=''
    new_id=None
    gff_lines=''
    cline=cfile.readline()
    while cline:
      try:
        if not cline[0]=="#":
          if tag=='*' or lower(cline.split('\t')[2]) == lower(tag):
            new_id=uniqid_function(cline)
            if current_id!= new_id:
              if current_id and gff_lines: #if it is not the first entry, append to the output list the gene we just finished parsing
                g=gene()
                g.load_gff(gff_lines, tag, **load_gff_args) 
                g.id=current_id
                if not add is None:  add(gff_lines, g)
                genes.append(g)
                gff_lines=''
              current_id=new_id
            gff_lines+=cline
      except Exception, e:
        print "ERROR loading gff  line: "+cline+" ### "+str(e)
        raise
      cline=cfile.readline()
    #last entry
    if current_id and gff_lines: #if it is not the first entry, append to the output list the gene we just finished parsing
      g=gene()
      g.load_gff(gff_lines, tag, **load_gff_args) 
      g.id=current_id
      if not add is None:  add(gff_lines, g)
      genes.append(g)

  cfile.close()
  return genes
  

def order_genes_for_strand_chr_pos(x, y):
  """Order fucntion for genes: by strand, chromosome, positions start
  """
  if x.strand!=y.strand:            return (x.strand+' '+y.strand).index('+')-1
  if x.chromosome!=y.chromosome:    return cmp(x.chromosome, y.chromosome)  
  else:                             return cmp(x.boundaries()[0], y.boundaries()[0])
order_genes_to_merge=order_genes_for_strand_chr_pos

def order_genes_for_chr_strand_pos(x, y):
  if x.chromosome!=y.chromosome:    return cmp(x.chromosome, y.chromosome)  
  if x.strand!=y.strand:            return (x.strand+' '+y.strand).index('+')-1
  else:                             return cmp(x.boundaries()[0], y.boundaries()[0])

def order_genes_for_chr_pos(x, y):
  if x.chromosome!=y.chromosome:    return cmp(x.chromosome, y.chromosome)  
  else:                             return cmp(x.boundaries()[0], y.boundaries()[0])
  
try:   
  from pygraph.classes.graph import graph
  from pygraph.algorithms.accessibility import connected_components

  class gene_overlap_graph(graph):
    """ Class to store overlaps between gene classes. Can't be initialized manually, it's only returned by function genes_overlap """
    def are_overlapping(self, g1, g2):
      """ Tell if gene1 and gene2 (arguments) are overlapping, meaning, are connected in the graph"""
      return g2 in self.node_neighbors[g1]
      
    def all_overlapping_with(self, g1):
      """ Return the list of genes overlapping with argument g1"""
      return self.node_neighbors[g1]
      
    def overlap_clusters(self, min_size=1, sort=True):
      """ returns a list of lists of genes overlapping each other.   
the genes are not required to be all overlapping to each other (e.g.   g1 overlaps g2 ;   g2 overlaps g3;    [ g3 do not overlap g1] -->  [g1, g2, g3] will be returned
Cluster minimal size is two. output list is sorted to have biggest clusters first, unless sort==False
      """
      cc= connected_components(self)  ## e.g .{a: 1, b: 1, c: 1, d: 1, e: 2, f: 2}
      if min_size>=2:
      ## removing clusters of size 1
        for g in self.nodes:
          if not self.all_overlapping_with(g): del cc[g]
      ## here cc has already at least two elements with same label, meaning: gene not overlapping with anything are not here
      ##### transforming label hash such as cc in lists of elements to output  
      connected_comp_lists_hash={}
      for k in cc: 
        if not connected_comp_lists_hash.has_key(   cc[k]   ):   connected_comp_lists_hash[ cc[k] ]= []
        connected_comp_lists_hash[ cc[k] ].append( k )  
        
      out=connected_comp_lists_hash.values()
      if min_size>=3:         out = [v for v in out if len(v)>=min_size]
      if sort:                out.sort(key=len)
      return out
 
except: pass

def bedtools_intersect(gene_list, strand=True, options={}): ## implicit:{'s':True} if strand==True
  """ Open a temp file preparing the input to bedtool intersect, runs it and returns a filehandler on the result. The file is provided to bedtools both as input A and input B, to compute all against all overlaps. Basic commadnline executed: #bedtool intersect -a ALL_GENES.fa -b ALL_GENES.fa  -wa -wb
  Options are passed as key of the dict options. When the value to a key is boolean True, that option has no argument. In any other case, the value is converted to string.
  E.g. bedtools_intersect(gene_list,  strand=True, options={'sorted':True, 'f':0.5 }  --> bedtool intersect -a ALL_GENES.fa -b ALL_GENES.fa  -wa -wb   -f 0.5 -sorted -s
  in output each gene is identified by its index in the list, e.g. chr4  9236903   9236953   11  0   +   chr4  9236903   9236944   12  0   + ### 11 and 12 are the ids
"""
  global temp_folder;   temp_overlap_file= temp_folder+'gene_overlap_file.bed'
  test_writeable_folder(temp_folder, 'temp folder ! Not defined maybe? [Use set_local_folders(temp_folder)]' )
  temp_overlap_file_h= open(temp_overlap_file, 'w')
  for g_index, g in  enumerate(gene_list):     print >> temp_overlap_file_h, g.bed( show_id = str(g_index), strand=strand )
  temp_overlap_file_h.close()  
  bedtools_intersect_command= "bedtools intersect {1} -wa -wb -a {0} -b {0} ".format(temp_overlap_file, {True:'-s', False:''}[bool(strand)])
  for k in options.keys():     
    bedtools_intersect_command+='-'+str(k)+' '
    if not options[k]==True:   bedtools_intersect_command+=str(options[k])+' '
#  print bedtools_intersect_command
  try:   bp=   bash_pipe(bedtools_intersect_command, return_popen=True)
  except IOError: 
    #raise Exception, 
    printerr("bedtools_intersect ERROR bedtools not installed??", 1)
    raise
  return bp #bbash(bedtools_intersect_command)
#  return bp    


def gene_clusters(gene_list, strand=True):
  """ Use bedtools intersect (temp_folder from MMib is used. See function set_local_folders) to compute overlaps between the genes.
  Returns two dictionaries: gene2cluster, cluster2genes
  where cluster is a numeric index (not consecutives)
"""
  geneid2cluster_index={}; cluster_index=1; cluster_index2geneids={}
  geneid2gene={}; 
  for g_index, g in enumerate(gene_list):    geneid2gene[str(g_index)]=g
  bp=bedtools_intersect(gene_list, strand=strand)
  for line in bp.stdout:
#    print line
    splt=line.rstrip().split('\t')
    if strand:    id_left = splt[3]; id_right= splt[9]
    else:         id_left = splt[3]; id_right= splt[7]
    if   id_left != id_right:
      if   (  not id_left in geneid2cluster_index )  and  ( not id_right in geneid2cluster_index ):  #new cluster
        cluster_index2geneids [cluster_index] = [id_left, id_right]
        geneid2cluster_index[id_left]=cluster_index;              geneid2cluster_index[id_right]=cluster_index
        cluster_index+=1
#        print '1 creating with ids: ',id_left, id_right, ' the cluster ', cluster_index-1
      elif (  id_left in geneid2cluster_index     )  and  ( not id_right in geneid2cluster_index ):
        cluster_index2geneids [ geneid2cluster_index[id_left] ].append( id_right )
        geneid2cluster_index[id_right]=   geneid2cluster_index[id_left]
#        print '2 moving id: '+id_right+' to cluster ', cluster_index
      elif ( not id_left in geneid2cluster_index  )  and  ( id_right in geneid2cluster_index ):
        cluster_index2geneids [ geneid2cluster_index[id_right] ].append( id_left )
        geneid2cluster_index[id_left]=    geneid2cluster_index[id_right]
#        print '3 moving id: '+id_left+' to cluster ', cluster_index
      elif geneid2cluster_index[id_left] != geneid2cluster_index[id_right]    : #not id_left in geneid2cluster_index  and  not id_right in geneid2cluster_index   is implicit
        #putting those in cluster of id_right into the cluster in id_left, unless the reverse is more efficient
        if len( cluster_index2geneids [ geneid2cluster_index[id_right] ] )  > len( cluster_index2geneids [ geneid2cluster_index[id_left] ] ): id_left, id_right = id_right, id_left
#        print '4 LEFT:', id_left,  'c:', geneid2cluster_index[id_left],  'RIGHT:',  id_right, 'c:', geneid2cluster_index[id_right]

        cluster_index_to_remove = geneid2cluster_index[id_right]
        for gid in cluster_index2geneids [ geneid2cluster_index[id_right] ]:   
#          print '4 moving id: '+gid, 'from cluster', cluster_index_to_remove, ' to cluster ', geneid2cluster_index[id_left]
          geneid2cluster_index[ gid ] =  geneid2cluster_index[id_left]
        cluster_index2geneids[  geneid2cluster_index[id_left]  ].extend(   cluster_index2geneids[  cluster_index_to_remove  ]   )
#        print '4 removing cluster', cluster_index_to_remove
        del cluster_index2geneids[  cluster_index_to_remove  ]

  # producing dictionaries that can be used in output
  gene2cluster={}; cluster2genes={}
  for gid in geneid2gene:
    g=geneid2gene[gid]
    if gid in geneid2cluster_index:    gene2cluster[ g ] = geneid2cluster_index[gid]
    else:                              gene2cluster[ g ] = cluster_index; cluster_index+=1
    if not gene2cluster[ g ] in cluster2genes: cluster2genes[ gene2cluster[ g ] ]= []
    cluster2genes[   gene2cluster[ g ]   ].append(g)
  return gene2cluster, cluster2genes


def remove_overlapping_gene_clusters(gene_list,  scoring=len,  cmp_fn=None, phase=False, strand=True, out_removed_genes=[], remember_overlaps=False, verbose=False, fix_ties=True,  overlap_fn=None  ):
  """ Returns a reduced version of the list in input, removing genes that overlaps.
  When two genes overlap, a score is assigned to each gene to decide which to keep -- similarly to the key argument to sort, you can provide a function as argument of scoring. by default, it's the gene lenght. Alternatively, you can use cmp_fn in a similar fashion to cmp in sort. It must accept two gene arguments, and return -1 if you want to keep the first argument, +1 if you want to keep the last one.  
  Important: when you use cmp_fn, take care of ties! don't let python decide, or your results may not be reproducible. when scoring is chosen, this is avoided with a trick here, which adds very small quantities (max 0.001) to the scores of each gene object which depend on their ids. 
  So, when using scoring, use integer scores or anyway make that in no case two gene objects must differ for more than 0.001!
  When you have complex overlap structures, this is what happens: 
  -clusters of overlapping genes are built. In this cluster though, you may have pairs of genes not overlapping, but overlapping to something that overlaps the other (or even more far fetched than that) --> e.g g1 overlaps g2 ;   g2 overlaps g3;    [ g3 do not overlap g1] --> a cluster is [g1, g2, g3]
  -for each cluster, initially the best scoring gene is taken (this will be output), and all things really overlapping with it are thrown away. The procedure is repeated until no gene is left in this cluster. 
  This will ensure that you have no overlaps in the output genes, and neither that you will take off genes without having something really overlapping with it in output 
  The list out_removed_genes may be used to collect the genes removed from the input list. To use it, initialize an empty list variable, then pass it as this argument: 
    # e.g.
    a_list=[]
    remove_overlapping_genes(gene_list, out_removed_genes=a_list)
    # now a_list contains the gene removed.
  When you use the out_removed_genes argument, you may want to know the correspondance between the genes removed and the ones kept, without recomputing overlaps. If you use remember_overlaps=True, the attribute .overlapping will be added to the removed genes; this is a link to the gene kept (which is present in the output, returned list)
   Normally the overlaps between any two genes is checked through two steps; first, the bedintersect tool, which can take into account the strand or not (depending on the argument of strand); second, the gene.overlaps_with function, which can take into account also the phase (frame). You can replace this second check with any given function providing it as argument of overlap_fn; this must take two gene arguments, and return True or False. If overlap_fn is provided, then the phase argument is ignored.
   """
  outlist=[]
  #overlaps_graph= genes_overlap(gene_list, phase=phase, strand=strand)
  #clusters= overlaps_graph.overlap_clusters(min_size=1)
  gene2cluster, cluster2genes = gene_clusters( gene_list, strand=strand )

  for cluster_id  in cluster2genes.keys():
    cluster= cluster2genes[cluster_id]
    #cluster is a gene list 
    ## we do the following:  we take the best scoring gene, we put this in the outlist, we throw away everything that overlaps with it, and we repeat until we finished the cluster
    while cluster:
      if cmp_fn:     
        cluster.sort(cmp=cmp_fn)
      else:          
        #I'm not giving directly scoring to sort because this may lead to different results in different runs, because of ties ( only if fix_ties is true           )
        hash_obj_to_score={}
        for obj in cluster:
          score=scoring(obj)
          if fix_ties:   score+= string_hashed_to_number( str(obj.id),  0.001 )
          hash_obj_to_score[obj]= score
          #print obj,  score
        cluster.sort(key= lambda x:hash_obj_to_score[x], reverse=True)           

      best_gene=  cluster[0]
      genes_not_overlapping_best_gene=[]
      outlist.append(best_gene)
      for g in cluster: 
        if g!=best_gene:
          they_overlap = best_gene.overlaps_with(g, phase=phase, strand=strand) if overlap_fn is None else  overlap_fn(best_gene, g)
          if not they_overlap: 
            ## gene that we're keeping for next cycle
            genes_not_overlapping_best_gene.append(g)
          else:  
            ## genes that we're throwing away
            if remember_overlaps:              g.overlapping= best_gene
            if verbose: printerr('removing: '+g.id+'  --> keeping: '+best_gene.id, 1)
            out_removed_genes.append(g)
      cluster=     genes_not_overlapping_best_gene
  return outlist


#gs=load_all_genes('/users/rg/mmariotti/Archaea/ivan_dotu/aSeblastian/janna_all_methods.knownsp.output.all_secis.gff', tag='secis')
#get_gene_overlaps(gs)
  
def genes_overlap(gene_list, phase=False, strand=True):
  try:    overlaps_graph= gene_overlap_graph()
  except NameError:  raise ImportError, "ERROR pygraph modules not installed! can't initialize subclass gene_overlap_graph "
  list_no_empty=[g for g in gene_list if g]
  overlaps_graph.add_nodes(list_no_empty)
  #ordered_gene_list= sorted(   list_no_empty, cmp=order_genes_for_chr_strand_pos )     ### old: doesn't work if strand==False
  ordered_gene_list= sorted(   list_no_empty, cmp=order_genes_for_chr_pos )  
  for index1, g1 in enumerate(ordered_gene_list):
    right_boundary_g1=g1.boundaries()[1]
    index2=index1+1
    while  index2<len(ordered_gene_list) and   ordered_gene_list[index2].boundaries()[0] <= right_boundary_g1:
      if g1.overlaps_with(  ordered_gene_list[index2],  phase=phase, strand=strand ):        overlaps_graph.add_edge( (g1, ordered_gene_list[index2]) )     
      index2+=1  
  return overlaps_graph

def remove_overlapping_genes(gene_list,  scoring=len,  cmp_fn=None, phase=False, strand=True, out_removed_genes=[], remember_overlaps=False, verbose=False, fix_ties=True):
  """ Returns a reduced version of the list in input, removing genes that overlaps.
  When two genes overlap, a score is assigned to each gene to decide which to keep -- similarly to the key argument to sort, you can provide a function as argument of scoring. by default, it's the gene lenght. Alternatively, you can use cmp_fn in a similar fashion to cmp in sort. It must accept two gene arguments, and return -1 if you want to keep the first argument, +1 if you want to keep the last one.  
  Important: when you use cmp_fn, take care of ties! don't let python decide, or your results may not be reproducible. when scoring is chosen, this is avoided with a trick here, which adds very small quantities (max 0.001) to the scores of each gene object which depend on their ids. 
  So, when using scoring, use integer scores or anyway make that in no case two gene objects must differ for more than 0.001!
  When you have complex overlap structures, this is what happens: 
  -clusters of overlapping genes are built. In this cluster though, you may have pairs of genes not overlapping, but overlapping to something that overlaps the other (or even more far fetched than that) --> e.g g1 overlaps g2 ;   g2 overlaps g3;    [ g3 do not overlap g1] --> a cluster is [g1, g2, g3]
  -for each cluster, initially the best scoring gene is taken (this will be output), and all things really overlapping with it are thrown away. The procedure is repeated until no gene is left in this cluster. 
  This will ensure that you have no overlaps in the output genes, and neither that you will take off genes without having something really overlapping with it in output 
  The list out_removed_genes may be used to collect the genes removed from the input list. To use it, initialize an empty list variable, then pass it as this argument: 
    # e.g.
    a_list=[]
    remove_overlapping_genes(gene_list, out_removed_genes=a_list)
    # now a_list contains the gene removed.
  When you use the out_removed_genes argument, you may want to know the correspondance between the genes removed and the ones kept, without recomputing overlaps. If you use remember_overlaps=True, the attribute .overlapping will be added to the removed genes; this is a link to the gene kept (which is present in the output, returned list)   """
  outlist=[]
  overlaps_graph= genes_overlap(gene_list, phase=phase, strand=strand)
  clusters= overlaps_graph.overlap_clusters(min_size=1)
  for cluster_index, cluster in enumerate(clusters):
    #cluster is a gene list 
    ## we do the following:  we take the best scoring gene, we put this in the outlist, we throw away everything that overlaps with it, and we repeat until we finished the cluster
    while cluster:
      if cmp_fn:     
        cluster.sort(cmp=cmp_fn)
      else:          
        #I'm not giving directly scoring to sort because this may lead to different results in different runs, because of ties ( only if fix_ties is true           )
        hash_obj_to_score={}
        for obj in cluster:
          score=scoring(obj)
          if fix_ties:   score+= string_hashed_to_number( str(obj.id),  0.001 )
          hash_obj_to_score[obj]= score
          #print obj,  score
        cluster.sort(key= lambda x:hash_obj_to_score[x], reverse=True)           

      best_gene=  cluster[0]
      genes_not_overlapping_best_gene=[]
      outlist.append(best_gene)
      for g in cluster: 
        if g!=best_gene:
          if not overlaps_graph.are_overlapping(best_gene, g):          
            ## gene that we're keeping for next cycle
            genes_not_overlapping_best_gene.append(g)
          else:  
            ## genes that we're throwing away
            if remember_overlaps:              g.overlapping= best_gene
            if verbose: printerr('removing: '+g.id+'  --> keeping: '+best_gene.id, 1)
            out_removed_genes.append(g)
      cluster=     genes_not_overlapping_best_gene
  return outlist

from hashlib import md5 

def string_hashed_to_number(astring, n_max=1.0):
  """ convert a string to a hashed number ranging from 0 to n_max (def:1) """
  m = md5(astring)
  number= int(   m.hexdigest(),   16  )  #getting the number corresponding to the hex code returned by the md5 hashing functino
  higher_possible_number= int( 'f'*32  ,   16  ) #highest possible number is all 'f' chars, and the string returned by hex digest is 32 chars long
  out= number*float(n_max) / higher_possible_number
  return out

def merge_genes(gene_list, phase=False, inplace=False, mode='merge', id_mode='SUM', program_mode='LONGEST', removed_genes=[], remember_overlaps=False, strand=True, strict_overlap=False):
  """ This functions accepts a list of gene objects, and it merges them (in place, if inplace==True): overlapping genes (checked with the overlaps_with function, in case with phase) are merged. If mode=="merge", when two genes overlap, they're replaced by their union. The gene.id gets the name of the two ids joined by "_+_" (this behavior can be changed with id_mode... see union_with function in gene class. If mode=="longest", only the longest prediction among the overlapping ones is kept, and it is not modified. As mode, one can also provide a function which has to accept two gene object arguments and return a new one. In this case, id_mode and program_mode are ignored
  removed genes accepts a list which is modified inplace to contain to list of removed genes
  If the variable  remember_overlaps is set to True, the genes returned in the removed_genes list contain an new attribute, .overlapping, which points to the overlapping gene which was kept
  When two genes are merged, the one which is kept will keep a list of the genes that were removed in its favor, and the overlaps will be computed also on this list of genes which it represents. For this reason, this may happen:   g1 overlaps with g2 (and has priority over it).  g2 overlaps with g3 but not with g1.       -->    g1 will be kept and g2 and g3 are removed -- even if g1 was not overlapping with g3.   This method ensure consistency of results. 
  anyway, if you want to turn this off, use strict_overlap=True. in this case, the overlaps are not computed with represented genes
  """
  if not strand: raise Exception, "ERROR this function cannot be used to merge genes on different strands! Please convert all of them to + and rerun."
  if not mode in ['merge', 'longest'] and type(mode)!=type(lambda x:x): #the mode can also be a function telling which gene to keep when two overlap
    raise Exception, "merge_genes function, mode not recognized : "+str(mode)
  exon_list=[] #### I build an exon list which I can order by position start, so then I can exploit the fact the only adiacent exons can overlap. 
  for g_index, g in enumerate(gene_list):
    for st, end in g.exons:
      e=gene(chromosome=g.chromosome, strand=g.strand)
      e.add_exon(st, end)
      e.id=g.id
      e.index=g_index           #In every exon, I add a ".index" attribute to keep track from which gene it comes from
      exon_list.append(e)
    gene_list[g_index].index=g_index   #then, I add the same attribute to the gene, to know its position in the gene list without navigating the whole list.
    gene_list[g_index].representing_genes=set()   #index collection. Every time two genes g1 and g2 are found overlapping and g1 is kept but g2 is not, the g1.representing_genes will contain the index of g2, remembering that now g1 represents also g2
  exon_list.sort(order_genes_to_merge)


  #now I cycle through each exon. At each cycle, there's one ruler, and the next one is checked. Anyway, since these are gene objects, I check the overlap of the full genes, instead of the single exons. NB: If gene1 and gene2 are overlapping, then one of the two is added a ".merged_in" attribute, which basically means, ignore this result and see gene n. (merged_in propriety) instead.
  for index_ruler in range(len(exon_list)-1):
    current_exon=exon_list[index_ruler]; next_exon=exon_list[index_ruler+1]
    gene1= gene_list[current_exon.index]; gene2=gene_list[next_exon.index]
    #dereferencing: if the exon belong to a gene which is been merged, consider the linked gene instead
    while not gene1['merged_in'] is None:        
      #print "delink1: "+gene1.id +' --> '+str(gene1['merged_in'])+' '+gene_list[gene1.merged_in].id
      gene1=gene_list[gene1.merged_in]    
    while not gene2['merged_in'] is None:        
      #print "delink2: "+gene2.id +' --> '+str(gene2['merged_in']  )+' '+ gene_list[gene2.merged_in].id
      gene2=gene_list[gene2.merged_in]    
    if gene1.index!=gene2.index:
      genes_overlaps=False
      if mode=='merge':         genes_overlaps=gene1.overlaps_with(  gene2, phase=phase, strand=strand  )
      else: # if mode != merge, then we must check if (any gene represented by)  gene1 overlaps with any gene represented by gene2. Break statements are to optimize
        for g1_ind in set( [ gene1.index ] ) | gene1.representing_genes :
          if genes_overlaps: break
          for g2_ind in  set( [ gene2.index ] ) | gene2.representing_genes:
            if genes_overlaps: break
            if gene_list[g1_ind].overlaps_with(gene_list[g2_ind], phase=phase, strand=strand): genes_overlaps=True
#      if gene1.overlaps_with(  gene2, phase=phase  ):                                        #################### replacing one of the two genes

      if genes_overlaps:
        #print 'merging '+gene1.id+ ' '+gene2.id
        #for frid in gene1.id.split('_+_'):
        #  if '_+_'+frid+'_+_' in gene2.id: 
        #    raise Exception, "time to die! "+frid
        if mode=='merge':
          gene_list[gene1.index] =  gene1.union_with( gene2, id=id_mode, program=program_mode, index=gene1.index ) #with the union
          gene_list[gene2.index].merged_in =  gene1.index
        elif mode=='longest':                   #with the longest of the two genes
          if gene1.length()>=gene2.length():                       #print "longest is "+gene1.id+' (against '+gene2.id+') ; linking '+gene2.id +' to '+str(gene1.id)
            gene_list[gene2.index].merged_in = gene1.index 
          else:                                                   #print "longest is "+gene2.id+' (against '+gene1.id+') ; linking '+gene1.id +' to '+str(gene2.id)     
            gene_list[gene1.index].merged_in = gene2.index            
        elif type(mode)==type(lambda x:x):
          #print "merging g1 ("+gene1.id+','+str(gene1.index)+') with g2 ('+gene2.id+','+str(gene2.index)+')',
          g = mode( gene1, gene2 )         ## evaluating which of the two overlapping genes should be kept
          #print ' --> '+g.id
#################################################
          #print 'merged_in', gene_list[gene1.index]['merged_in'],  gene_list[gene2.index]['merged_in'], , gene2==g
          # may return gene1, gene2 or a new gene. if it is gene1 or gene2, I put this control to be sure that the .merged_in attribute is never pointing to itself.
          if gene2==g:   #keeping gene2
            g=g.copy();    g.index=gene2.index
            for gene_index in gene_list[gene1.index].representing_genes | set([gene1.index]):
              gene_list[gene_index].merged_in=g.index
              if not strict_overlap:               g.representing_genes.add(gene_index)
            gene_list[gene2.index] =  g            
            gene_list[gene1.index].representing_genes=set()
          else:
            g=g.copy();    g.index=gene1.index
            for gene_index in gene_list[gene2.index].representing_genes | set([gene2.index]):
              gene_list[gene_index].merged_in=g.index
              if not strict_overlap:               g.representing_genes.add(gene_index)
            gene_list[gene1.index] =  g            
            gene_list[gene2.index].representing_genes=set()
          
  list_out=[]; removed_gene_indexes=[]
#  unlinked_list=[] #will contain all the genes to reconsider, since : gene1 was merged into a gene2 that was  merged into gene3 but now gene1 does not overlap with gene3. it is possible only if mode != merge
  for g_index in range(len(gene_list)):  
    del gene_list[g_index].index  #deleting non.std propriety
    if gene_list[g_index]['merged_in']!=None :  #outputing all those which are not merged_in something else.
#      if mode!='merge':
#        #checking that the gene in which this was merged is really overlapping this gene. If the mode is not merge, it may not. In this case, we have to make a subsequent analysis
#        final_gene=gene_list[ gene_list[g_index]['merged_in'] ]
#        while final_gene['merged_in']!=None:
#          final_gene=gene_list[final_gene['merged_in']]
#        if not  gene_list[g_index].overlaps_with(final_gene, phase=phase):
#          gene_list[g_index]['merged_in']=None
#          unlinked_list.append(gene_list[g_index])
      final_gene=gene_list[ gene_list[g_index].merged_in ]
      if remember_overlaps:      gene_list[g_index].overlapping= final_gene

      if bool(gene_list[g_index].representing_genes): raise Exception, "merge don't work"
      removed_genes.append(gene_list[g_index])
      removed_gene_indexes.append(g_index)
      del gene_list[g_index]['merged_in']
    else:
      if not inplace:        list_out.insert(0, gene_list[g_index])
      if inplace:        gene_list.pop(g_index)         #print "removing ", g_index
      
    del gene_list[g_index]['representing_genes']

#  if unlinked_list:    unlinked_list=merge_genes( unlinked_list, phase=phase, inplace=False, mode=mode, id_mode=id_mode, program_mode=program_mode ) 
  if inplace:    
    remove_items_from_list_with_indexes(gene_list, removed_gene_indexes, inplace=True)
    return gene_list
  else:          return list_out
    
def get_species_from_library(species_name, library=False):
  """ Parse the file species_library and returns a list: [taxid, scientific name] for the species name provided. This can differ from the name in output is the title provided is listed under the synonyms of the actual species name. If the species name is not found, returns None.
  It accepts also species names with masked characters as {chXX}, and also detects those names in which underscore means white space.
  For this method to work, the variable species_library must be accessible, and contain a string pointing to a names.dmp file as downloaded from ncbi taxonomy.
  Alternatively this file can be provided with the keyword argument library.    """
  if '_' in species_name  and not ' ' in species_name:    species_name=replace_chars(species_name, '_', ' ')
  species_name=unmask_characters(species_name)
  if not library: library=species_library
  cmnd='gawk  -v input="'+species_name+'"'+""" -F"\\t" '{ if ($7=="scientific name"){ sc_name=$3; sc_name_taxid=$1};  if ($3 == input){ taxid=$1 }; if (taxid && taxid==sc_name_taxid){print taxid "\\t" sc_name; exit}   }'    """+library
  b=bash(cmnd)
  if b[0]: raise notracebackException, "ERROR searching species. something wrong in  command: "+cmnd
  if not b[1]: return None
  else:        return b[1].split('\t')



def get_taxids_from_ncbi_db(   scientific_name_list,    ncbi_db='/users/rg/didac/NCBI/Taxonomy/names.dmp',   temp_dir=None, silent=False, full=False    ):
  """ Utility to get the taxids for a (high) number of scientific names.  Input is:  list of scientific names (not synonyms!), path to names.dmp file of ncbi_taxonomy, temporary folder (or temp_folder of MMlib will be used)  and silent flag.   If ambygous entries are present or some name is not found, warning are printed to screen, unless silent==True.
  Normally, an hash with {scientific_name -> taxid} is returned. If full is set to True, a tuple is returned like: out_hash, ambigous_hash, not_found_hash
  out_hash is like normal output. ambigous_hash collects the taxids of ambigous entries, like   {scientific_name -> [taxid1, taxid2, .. ]}.   not_found_hash collects the scientific names that were not found, like: {scientific_name -> 1 }
  NOTE: taxids returned are of type  int    """

  if temp_dir is None: temp_dir = temp_folder
  t_file=    temp_dir+'/sc_names_w_tab.txt' ; t_file_h=  open(t_file, 'w')
  for sc_name in scientific_name_list:    print >> t_file_h, "|\t"+sc_name+"\t|"
  t_file_h.close()
  grep_out= temp_dir+'/grep_out'
  bash( 'grep -Ff '+t_file+' ' +ncbi_db+' > '+grep_out)
  ambygous= temp_dir+'/ambygous_out'
  unambygous= temp_dir+'/unambygous_out'
  bash('rm '+ambygous)
  bash("gawk -F'|' '{if ($3!=\"\t\t\"){print $0 > \""+ambygous+"\"} else{ print } }' "+grep_out + " > "+unambygous)
  if is_file(ambygous):
    n_ambygous_entries   = int(  bash('wc -l '+ambygous)[1].split()[0]    )
    n_sc_names_ambygous  = int(  bash("  gawk -F'|' '{print $2}' "+ambygous+" | sort | uniq  | wc -l "  )[1].split()[0]  )
  else:   n_ambygous_entries,  n_sc_names_ambygous= 0, 0
  n_unambygous_entries = int(  bash('wc -l '+unambygous)[1].split()[0]  )
  if not silent and n_ambygous_entries:     
    printerr("get_taxids_from_ncbi_db WARNING "+str(n_sc_names_ambygous)+ " scientific names in input are ambygous: a total of "+str(n_ambygous_entries)+ " entries were found for these. Excluding them from output...", 1 )
  n_sc_names_found= n_sc_names_ambygous+ n_unambygous_entries
  if not silent and n_sc_names_found != len(scientific_name_list):
    printerr("get_taxids_from_ncbi_db WARNING "+str(  len(scientific_name_list) -  n_sc_names_found  )+ " scientific names were not found!", 1 ) 
  out_hash={}
  for line in open(unambygous):
    splt = line.split('\t')  
    out_hash[  splt[2]   ] =   int(splt[0])
  if full: 
    ambigous_hash={}
    if n_ambygous_entries:
      for line in open(ambygous):
        splt = line.split('\t')  
        if not  ambigous_hash.has_key(splt[2]):  ambigous_hash[  splt[2]   ] =   []
        ambigous_hash[splt[2]].append(       int(splt[0])    )
    not_found_hash={}
    for n in scientific_name_list:
      if not ( n in out_hash or  n in ambigous_hash): not_found_hash[n]=1
    return out_hash, ambigous_hash, not_found_hash
  return out_hash
    
class parser(object):
  """ This class handles reading from a text file, typically a sequence file. IT is meant to be a parent class for parsers. 
  Usage: (example with parse_fasta)
  p=parse_fasta(file) #--> this opens the file, prepares a filehandler ; NOTE: a filehandler can be also provided instead of a file
  for object in p:    #--> the next() method is used until a StopIteration is raised. Inside it, the method parse_next() is run. This is different for each parser metaclass, and returns object
    print object      # [title, sequence]
  
  Example of parse_next (of parse_fasta parser):

  def parse_next(self):
    title=self.last_line[1:-1]
    seq=''
    self.last_line=self.file.readline()                    # ---> read a new line, stores it in self.last_line
    while self.last_line and  self.last_line[0]!='>':      # ---> when self.last_line is false (end of the file), next() will raise StopIteration instead of calling parse_next()
      seq+=replace_chars(self.last_line, ' \n\r\t', '')
      self.last_line=self.file.readline()
    return title, seq    
  
  To build a new parser, define a child class of this superclass and define a parse_next method which parse self.last_line or, in case, the next lines and returns a desired object.
  """
  def __getitem__(self, name):
    if name in  dir(self):       return self.__dict__[name]
    else:                        return None
  def __setitem__(self, key, value):
    self.__dict__[key]=value
  def __init__(self, filename='',  **keyargs):
    for key in keyargs:      self[key]=keyargs[key]
    if filename:             self.load(filename) 
  def load(self, filename=''):
    if not filename and self.file and type(self.file) == file: raise Exception, "parser ERROR: trying to load a filehandler which has already finished. Please instanciate another parser"
    if type(filename)==file:
      self.load_filehandler(filename)
    else: #string
      self.load_filename(filename) 
  def load_filename(self, filename=''):
    if not filename:      filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file= open(filename, 'r')
    self.last_line=self.file.readline()
  def load_filehandler(self, fileh=''):
    self.file=fileh 
    self.last_line=self.file.readline()
  def __iter__(self):    return self      
  def all(self):
    outlist=[]
    for i in self:      outlist.append(i)
    if outlist==[None]: return []
    return outlist
  def next(self, skip_comments=True):
    if self.file.closed:      self.load()
    if not self.last_line:    self.stop()
    if skip_comments and not   ( self.__dict__.has_key('skip_comments') and not self['skip_comments'] )  :
      while self.last_line and (self.last_line[0]=='#'):                self.last_line=self.file.readline()  
    try:      return self.parse_next()
    except StopIteration: raise
    except Exception, e:
      print "ERROR parsing file: "+self.file.name+' ; line: '+self.last_line
      raise          
  def stop(self):
    self.file.close()
    raise StopIteration

  def parse_next(self):
    """This method is the key of the parser class and must be implemented for each parser. The method should read self.last_line and next line (in case it is necessary). Before returning the desired object, it should move the cursor self.last_line to the next line
    """
    raise Exception, "ERROR the generic parser class has no parse_next method: you must define it in the metaclass"

class parse_fasta(parser):
  remove_chars=set(['\n', '\r', '\t', ' '])
  def parse_next(self):
    title=del_white(self.last_line[1:-1])
    seq=''
    self.last_line=self.file.readline()
    while self.last_line and  self.last_line[0]!='>':
      seq+=replace_chars(self.last_line, self.remove_chars, '')
      self.last_line=self.file.readline()
    return title, seq    

class parse_sam(parser):
  """ Returns a gene object for each line of the sam input. nornmally just positions, strand, chromosome, id and sequence are kept. if you want other attributes, define them as attributes of this parser object. You can have these attributes: 'flag', 'mapq', 'cigar', 'rnext', 'pnext', 'qual'] ( see SAM1 manual)
  Example, to have the qualities:
  p=parse_sam( somefile )
  p.qual=True
  for g in p:
    print g.id, g.qual # --> the qualities are stored in the .qual attribute
  """
  def parse_next(self):
    #write(self.last_line, 1, how='blue')
    splt=    self.last_line.split('\t')
    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = splt[:11]
    flag, pos, mapq, pnext, tlen= int(flag), int(pos), int(mapq), int(pnext), int(tlen)
    g=gene(chromosome=rname, strand='+')
    g.add_exon(  pos, pos+len(seq)-1 )
    g.seq=seq 
    g.id=qname
    for i in ['flag', 'mapq', 'cigar', 'rnext', 'pnext', 'qual']:
      if hasattr(self, i) and getattr(self, i):        exec('g.'+i+' = '+i)
    self.last_line=self.file.readline()
    return g


#class parse_gff(parser):
#  def parse_next(self):


class rnazhit(object):
  """ This class read input from RNAz output as channeled by the program rnazWindows.pl run on clustalw format nucleotide alignments """

  features_names=["Sequences", "Columns", "Reading direction", "Mean pairwise identity", "Shannon entropy", "G+C content", "Mean single sequence MFE", "Consensus MFE", "Energy contribution", "Covariance contribution", "Combinations/Pair", "Mean z-score", "Structure conservation index", "Background model", "Decision model", "SVM decision value", "SVM RNA-class probability", "Prediction"]

  def __init__(self, rnaz_string=''):
    self.data={} #stores every numeric feature in the header of rnaz output, which is like:   Mean z-score:  -0.09   
    # we index by the key: "Mean z-score"  -- value is converted to approprioate type e.g. -0.09 -> float
    self.ali=alignment(); self.title2ss={}; self.title2zscore={}; self.title2mfe={};   self.title2code={};
    if rnaz_string:   self.load(rnaz_string)

  def load(self, rnaz_string):
    """ parse a string coming from parsing a rnaz output file"""
    self.__init__()  #resetting data
    line_index=0; lines= rnaz_string.split('\n')
    ##parsing "key: value" file header
    while line_index<len(lines) and not lines[line_index].strip().startswith("Sequences"):      line_index+=1
    for line in lines[line_index:]:   
      splt=line.split(':')
      if len(splt)<2: break
      self.data[splt[0].strip()]= option_value( splt[1].strip() )
    ##checking for missing values
    err_msg=''
    for key in self.features_names:
      if not self.data.has_key(key): err_msg+=key+', '
    if err_msg: raise Exception, "rnazhit load ERROR feature"+int(err_msg.count(',')>1)*"s"+" not found: "+err_msg[:-2]
    ### parsing alignment
    while line_index<len(lines) and not lines[line_index].strip().startswith(">"):      line_index+=1
    title, seq, ss= None, None, None
    for line in lines[line_index:]:
      if   line.startswith(">"):       title=line.strip()[1:]
      elif title and not seq:          seq  =line.strip()
      elif title and seq and not ss:   
        ### line is now: #   ..(((--))..)..ETCETERA..((.-))...  ( -447.90, z-score =  -0.85, S)
        ss   =line.split()[0]      
        self.ali.add(title, seq)
        self.title2ss[title]= ss
        if title!='consensus': 
          self.title2zscore[title]=float(line.split("z-score = ")[1].split(",")[0])
          self.title2mfe[title]=float(line.split(",")[0].split()[-1]) 
          self.title2code[title]=line.rstrip().split()[-1][0] # R or S
        title, seq, ss= None, None, None
      if line.startswith("#"): break
    #checking evertyhing is there
    n=self.ali.nseq()-1
    if n!= len(self.title2zscore) or n!= len(self.title2mfe) or n!=len(self.title2code) or n!=self.data['Sequences']: 
      print self.summary()
      raise Exception, "rnazhit load ERROR some data was not found for all "+str(self.data['Sequences'])

  def positions(self):
    """ if RNAz was run channeled through rnazWindows (if not, return None) , returns two indexes which are found in the sequence titles; indexes are python style, 0 for first and end not included. e.g. (0, 120) """
    try:          return  map(int,  self.titles()[0].split('/')[-1].split('-'))
    except:       return None


  def __nonzero__(self):                 return bool(self.data)
  def probability(self):
    if not self: return None
    return self.data["SVM RNA-class probability"]
  def MFE(self):
    if not self: return None
    return self.data["Consensus MFE"]
  def zscore(self):
    if not self: return None
    return self.data["Mean z-score"]
  def titles(self):    return self.ali.titles()
  def seq_of(self, t): return self.ali.seq_of(t)
  def summary(self): 
    """returns the (almost) same string parsed to generate this hit"""
    if not self: return None
    o='############################  RNAz ?.?  ##############################\n\n'
    for k in self.features_names:      o+=' '+k+': '+str(self.data[k])+'\n'
    o+='\n######################################################################\n\n'
    for t in self.titles(): o+='>'+t+'\n'+self.seq_of(t)+'\n'+self.title2ss[t]+'\n'
    return o
  __repr__= summary

class parse_rnaz(parser):
  """ Parse a rnazhit object at each next() call"""
  def load(self, filename=''):
    if not filename:
      filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file=open(filename, 'r')
    self.last_line=self.file.readline()
    while self.last_line and not (self.last_line.startswith("#") and "RNAz" in self.last_line):
      self.last_line=self.file.readline()
    #set first line to ############################  RNAz 2.1  ##############################
      
  def parse_next(self):
    rnaz_string=''
    self.last_line=self.file.readline() # skipping line ####(..)####  RNAz   so we can parse until the next one. if there's none, self.last_line will be set on false and this method won't be called anymore
    while self.last_line and not (self.last_line.startswith("#") and "RNAz" in self.last_line):
      rnaz_string+=self.last_line
      self.last_line=self.file.readline()
    return rnazhit(rnaz_string)
      
      
class estexoneratehit(gene):
  """ puppet class to manage est2genome exonerate predictions"""
  def gff(self):    return gene.gff(self, tag='exon', program='exonerate')

class parse_exonerate_est(parser):
  """ Parse an exonerate file and returns a exoneratehit object for each prediction inside. 
  The target names obtained by fastasubseq are recognized and set to absolute coordinates. Also the target names obtained through the method from gene class "fasta_sequence", with title set to "fasta_title", are recognized (they are used in selenoprofiles).
  """
  def load(self, filename=''):
    if not filename:  filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file = open(filename, 'r')
    self.last_line=self.file.readline()

  allowed_chars={'A':1, 'G':1, 'T':1, 'C':1, 'N':1, '-':1}    
  def parse_next(self):
    line=self.last_line;     cfile=self.file #necessary to recycle code: this below is the old exonerate parser
    # reading inputfile
    cont_b=0
    while line and line !='-- completed exonerate analysis\n':
      full_query_block=''; full_target_block=''; 
      ali_target_seq='';   ali_query_seq=''; 
      while line and line !='-- completed exonerate analysis\n' and line.split(':')[0]!="  Target range":
        if line!='\n' and line.split()[0]=='Target:':   full_target_name=del_white(line[:-1].split('Target:')[1])
        line=cfile.readline()
      if line and line !='-- completed exonerate analysis\n':
        skip_intron=False
        passed_the_N_of_intron=False
        line=cfile.readline()
        line=cfile.readline()
        while line and line.split(':')[0]!='vulgar':
          align_block=[]        # reading 4 lines-alignment blocks
          for i in range(4):
            align_block.append(line)
            line=cfile.readline()
          ali_block_start= align_block[0].find(':')
          ali_block_end  = align_block[0].rfind(':') 
          query_block=   align_block[0] [ ali_block_start+2:ali_block_end-1  ]
          target_block=  align_block[2] [ ali_block_start+2:ali_block_end-1  ]

          full_query_block  += query_block
          full_target_block += target_block

        if len(full_query_block)!=len(full_target_block): raise Exception, "ERROR parsing exonerate_dna ! different lengths of full_query_block and full_target_block: \n"+full_query_block+'\n'+full_target_block
        for index, q_char in enumerate( full_query_block ):
          if q_char in self.allowed_chars and full_target_block[index] != '.':
            ali_query_seq += q_char
            ali_target_seq  += full_target_block[index]
            
        vulgar_line=line
        qname=vulgar_line.split()[1] ;      qstart=int(vulgar_line.split()[2])    + 1   ;     qend=int(vulgar_line.split()[3])   
        tname=vulgar_line.split()[5];       real_tstart=int(vulgar_line.split()[6])    +1 ;    real_tend=int(vulgar_line.split()[7])
        if real_tstart> real_tend:    #negative frame
          real_tstart = real_tstart-1
          real_tend = real_tend+1
        raw_score=vulgar_line.split()[9]
        ali=join(vulgar_line.split()[9:], ' ')
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

        e=estexoneratehit()
        e.load_gff(gff_lines, tag='exon')
        e.alignment=alignment()
        e.alignment.add('q', ali_query_seq );       e.alignment.add('t', ali_target_seq )
        e.query=gene(chromosome=qname, strand='+')
        e.query.add_exon(qstart, qend)
        e.id='est_exonerate:'+str(uniq_id(e))
        e.query.id=str(uniq_id(e))+'_query'
        e.score=int(raw_score)
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

        if e.length()!= len(nogap(e.alignment.seq_of('t'))):  
          printerr('WARNING error in the exonerate file: '+(self.file.name)+' ; the gff does not have a perfect correspondence with the alignment displayed. To avoid crash, returning None like if it were an empty exonerate prediction', 1)
          return None
        return e
        
        
  
class infernalhit(gene):
  """This class manages the predictions by infernal rna search program. The output accepted is the one by cmsearch  """   
  def __init__(self, infile=''):
    self.ss=''
    self.query=gene()
    self.score=None
    self.evalue=None
    self.alignment=alignment()
    self.cmali=None
    gene.__init__(self)
    if infile: self.load(infile)

  def load(self, infile):
    """ accepts a file with a single hit """
    self.__dict__ = parse_infernal(infile).all()[0]

  def summary(self):
    """ """
    out='CM: '+str(self.query.chromosome)+'\n'
    out+='>'+self.chromosome+'\n\n'
    out+=' Strand = '+self.strand+'\n'
    out+=" Query = "+str(self.query.boundaries()[0])+' - '+str(self.query.boundaries()[1])+', Target = '
    if self.strand=='+': out+=str(self.boundaries()[0])+' - '+str(self.boundaries()[1])+'\n'
    else:                out+=str(self.boundaries()[1])+' - '+str(self.boundaries()[0])+'\n'
    out+=" Score = "+str(self.score)
    if not self.evalue is None: out+=", E = "+str(self.evalue)
    out+='\n\n'
    out+=' '*12+self.ss+'\n'
    out+=' '*12+self.alignment.seq_of('q')+'\n'
    out+='\n' #consensus line... too lazy to do it
    out+=' '*12+self.alignment.seq_of('t')+'\n'
    out+='\n'
    return out

  def abstract_line(self):
    """One line with all the information for the infernal """
    out=self.query.chromosome+' '+self.chromosome+' '+self.strand+' TPOS:'+self.positions_summary()+' QPOS:'+self.query.positions_summary()+' S:'+str(self.score)+' E:'+str(self.evalue)+' '+self.alignment.seq_of('q')+' '+self.alignment.seq_of('t')+' '+self.ss
    return out

  def load_abstract_line(self, line):
    self.__init__()
    splt=line.split()
    self.query.chromosome=splt[0];  self.chromosome=splt[1]; self.strand=splt[2]; self.add_exons_from_positions_summary( splt[3].split('TPOS:')[1] )
    self.query.add_exons_from_positions_summary( splt[4].split('QPOS:')[1] ); self.score= float(splt[5].split('S:')[1]); 
    e=splt[6].split('E:')[1]
    if e!='None': self.evalue = e_v(e)
    self.alignment.add('q', splt[7] );     self.alignment.add('t', splt[8] );    self.ss = splt[9]     

  def remove_Xs(self, gene_seq=None):
    """This function is for the infernal hits which contain insertions, kept as Xs in the virtual infernalhit object. When run, it interrogates the target file specified in .target and recovers the missing sequence. """
    if not self.target: raise Exception, "infernalhit -> remove_Xs ERROR the .target attribute is not defined"

    if not gene_seq is None: gene_seq= replace( lower( gene_seq  ), 't', 'u')
    if "x" in self.alignment.seq_of('t'):
      gaps_target=0
      for pos in range(self.alignment.length()):
        nt=self.alignment.seq_of('t')[pos]
        if nt=='-': gaps_target+=1
        elif nt=='x': 
          p=pos
          while p<len(self.alignment.seq_of('t')) and self.alignment.seq_of('t')[p]=='x': p+=1
          x_range_start = pos
          x_range_end   = p-1
          break
      x_range = self.subseq( x_range_start-gaps_target+1, x_range_end-x_range_start+1  )
      if gene_seq is None:        
        subseq_from_target= replace( lower( x_range.fasta_sequence()[1]), 't', 'u')
      else:                       subseq_from_target=  gene_seq [x_range_start-gaps_target:x_range_end+1-gaps_target] 
      if len(subseq_from_target)!=x_range.length(): 
        raise Exception, "infernalhit->remove_Xs ERROR the sequence fetched for this hit has wrong length: subseq from target: {0}  range analyzed: {1} ".format(len(subseq_from_target), x_range.length())

      new_seq_target = self.alignment.seq_of('t')[:x_range_start]+subseq_from_target+self.alignment.seq_of('t')[x_range_end+1:]
      self.alignment.set_sequence('t',        new_seq_target       )
      return self.remove_Xs(gene_seq=gene_seq)
    if "x" in self.alignment.seq_of('q') and not self.cmali is None:
      gaps_query=0
      for pos in range(self.alignment.length()):
        nt=self.alignment.seq_of('q')[pos]
        if nt=='-': gaps_query+=1
        elif nt=='x':
          p=pos
          while p<len(self.alignment.seq_of('q')) and self.alignment.seq_of('q')[p]=='x': p+=1
          x_range_start = pos
          x_range_end   = p-1
      rf_no_insertions=replace(self.cmali.rf, '.', '')
      pos_in_rf=self.query.boundaries()[0]+ len( replace(self.alignment.seq_of('q')[:x_range_start], '-', '')  )-1   #this var is 0 based 
      rf_subseq= rf_no_insertions[pos_in_rf:pos_in_rf+x_range_end-x_range_start+1]
      new_seq_query=self.alignment.seq_of('q')[:x_range_start]+rf_subseq+self.alignment.seq_of('q')[x_range_end+1:]
      self.alignment.set_sequence('q',        new_seq_query       )
      return self.remove_Xs(gene_seq=gene_seq)

  def get_pairs(self, unaligned=False):
    """returns the list of pairs in the target, as 0 based positions in the alignment (or in the target sequence if unaligned=True) """
    pairs_in_model = ss_pairs( self.ss )
    out=[ [first, second]   for first, second in pairs_in_model  if not ( self.alignment.seq_of('t') [first] == '-' or self.alignment.seq_of('t') [second] == '-' ) ]
    if unaligned: out=[  [self.alignment.position_in_seq('t', first+1)-1, self.alignment.position_in_seq('t', second+1)-1]       for first, second in out ]
    return out
      
    
  def sequence(self):
    return upper(nogap(self.alignment.seq_of('t')))

  def RNApackage_ss(self):
    """ return a fasta like string with sequence in target and secondary structure"""
    a=alignment()
    press=self.ss
    ss=''
    for char in press:
      if    char in '<([{': ss+='('
      elif  char in '>)]}': ss+=')'
      else: ss+='.'
    seq=upper(self.alignment.seq_of('t'))
    find_index=seq.find('-')
    while find_index!=-1:
      # we're removing one character at a time: find_index
      seq=seq[:find_index]+seq[find_index+1:]
      ss=press[:find_index]+ss[find_index+1:]
      find_index=seq.find('-')
    return ">"+self.header()+'\n'+seq+'\n'+ss

def ss_pairs(ss_string):
    """ Given a string with secondary structure  -- in which pairs are represented by any parenthesis ([{<  -- returns a list of tuples with 0 based positions of each pair, from 5' to 3'"""
    height = 0;     stem = {};     pairs = []; #    pos_rf=0;     pos_t=0
    for pos, ss_char in enumerate(ss_string):
        if ss_char in '([{<':
            stem[ height ] = pos
            height += 1
        elif ss_char in ')]}>' and height > 0:
            height -= 1
            if height < 0: break
            paired_pos = stem[height]
            pairs.append( (paired_pos, pos) )
    return sorted(pairs, key = lambda x: x[0])

class parse_infernal(parser):
  """ Parse a cmsearch output (infernal rna search package) and return infernalhit objects at each next() call.
  It tries to identify the infernal version by parsing the first commented lines and look for something like: # INFERNAL 1.1.1 (July 2014)
  If your output does not have this, you can force version with this:
  the_parser = parse_infernal(  your_file )
  the_parser.infernal_version= '1.1' #for example
  for hit in the_parser:   #now this should work if version is correct
    #do stuff with hit
"""
  def load(self, filename=''):
    if not filename:      filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file=open(filename, 'r')
    self.current_cm='unknown'
    self.current_chromosome='unknown'
    self.infernal_version=None
    self.last_line=self.file.readline()
    while self.last_line and self.last_line.startswith('# '):
      if self.last_line.startswith('# INFERNAL '):        self.infernal_version=self.last_line.split()[2][:3]
      self.last_line=self.file.readline()
      
  def parse_next(self):
    if    self.infernal_version == '1.0': return self.parse_next_ver1_0()
    elif  self.infernal_version == '1.1': return self.parse_next_ver1_1()
    else: raise Exception, "ERROR infernal version not recognized or supported: "+str(self.infernal_version)

  def parse_next_ver1_1(self):
    while self.last_line and not self.last_line.startswith('>>'):
      self.last_line=self.file.readline()
    if self.last_line.rstrip() and self.last_line.startswith('>>'):
      self.current_chromosome=self.last_line.rstrip().split()[1]
      self.last_line=self.file.readline()
      hit_lines = []
      while self.last_line and not self.last_line.startswith('>>'):
        hit_lines.append( self.last_line.rstrip() )
        self.last_line=self.file.readline()
    elif not self.last_line.rstrip(): self.stop()

    g=infernalhit()
    g.chromosome=self.current_chromosome
    scores = hit_lines[2].rstrip().split()
    g.evalue =     e_v(scores[2])
    g.score =    float(scores[3])
    query_start =  int(scores[6])
    query_end =    int(scores[7])
    target_start = int(scores[9])
    target_end =   int(scores[10])
    g.strand =         scores[11]
    hand_rf = ''
    if '-' in g.strand:
      target_start, target_end = target_end, target_start
    g.add_exon(target_start, target_end)
    query_seq=''; target_seq=''
    while hit_lines:
      current_line = hit_lines[0]
      if current_line.endswith('CS'):
        g.ss += hit_lines[0].split()[0].rstrip('CS')
        self.current_cm = hit_lines[1].split()[0]
        query_seq += ' '.join( hit_lines[1].split()[2:-1] )
        target_seq += ' '.join( hit_lines[3].split()[2:-1] )
        ##  new, thanks to didac:
        del hit_lines[0:4]
        for next_line in hit_lines:
          if next_line.endswith('RF'):
            hand_rf += ' '.join( next_line.split()[0:-1] )
            break
        del hit_lines[0:]
        ##
      if hit_lines: del hit_lines[0]

    g.query.chromosome = self.current_cm
    g.query.strand='+'
    g.query.add_exon(query_start, query_end)

    while "*" in query_seq:
      pos_insert=query_seq.find('*')
      insert_length_target=  int( target_seq.split('[')[1].split(']')[0].strip()  )
      insert_length_query= int( query_seq.split('[')[1].split(']')[0].strip() )
      if scores[8].startswith('~'):
        scores[8]=''
        query_seq = insert_length_query*'x' + query_seq[pos_insert+1:]
        target_seq = insert_length_query*'-' + target_seq[pos_insert+1:]
        g.ss = '~'*max([insert_length_query, insert_length_target]) + g.ss.lstrip('~')
        if hand_rf:
          hand_rf = '~'*max([insert_length_query, insert_length_target]) + hand_rf
      else:
        query_seq = query_seq.split('*')[0]+ insert_length_query*'x'+"-"* (insert_length_target-insert_length_query)   + join( query_seq.split('*')[2:], '*'  )
        target_seq=target_seq.split('*')[0]+insert_length_target*'x'+"-"*(insert_length_query-insert_length_target)   + join( target_seq.split('*')[2:], '*'  )
        l_ss_gap=0
        while len(g.ss)>pos_insert+l_ss_gap and g.ss[pos_insert+l_ss_gap]=='~': l_ss_gap+=1
        g.ss = g.ss[:pos_insert]+'~'*max([insert_length_query, insert_length_target])+g.ss[pos_insert+l_ss_gap:]
        if hand_rf:
          hand_rf = hand_rf[:pos_insert]+'~'*max([insert_length_query, insert_length_target])+hand_rf[pos_insert+l_ss_gap:]

    g.alignment.add('q', replace(query_seq, '.', '-'))
    g.alignment.add('t', replace(target_seq, '.', '-'))
    g.hand_rf = hand_rf
    return g

  def parse_next_ver1_0(self):
    while self.last_line and not self.last_line.startswith('>') and not self.last_line.startswith('CM:') and not self.last_line.strip().startswith('Query') :      self.last_line=self.file.readline()
    if self.last_line.startswith('CM'): self.current_cm= self.last_line.rstrip().split('CM: ')[1]

    while self.last_line and not self.last_line.startswith('>') and not self.last_line.strip().startswith('Query'):      self.last_line=self.file.readline()
    if self.last_line.startswith('>'):    self.current_chromosome=self.last_line.rstrip()[1:]

    while self.last_line and not self.last_line.strip().startswith('Query') :      self.last_line=self.file.readline()
    if not self.last_line:       self.stop()    
    #now line is like:  Query = 1 - 86, Target = 9 - 89
    query_start = int(self.last_line.split('Query =')[1].split('-')[0].strip())
    query_end = int(self.last_line.split('Query =')[1].split(',')[0].split('-')[1].strip())
    target_start= int(self.last_line.split('Target =')[1].split('-')[0].strip())
    target_end = int(self.last_line.split('Target =')[1].split('-')[1].strip())

    g=infernalhit()
    g.chromosome=self.current_chromosome
    if target_start<target_end: g.strand='+'
    else:                       
      g.strand='-'
      target_start, target_end = target_end, target_start
    g.add_exon(target_start, target_end)

    g.query.chromosome = self.current_cm
    g.query.strand='+'
    g.query.add_exon(query_start, query_end)

    self.last_line=self.file.readline()
    #now :  Score = 22.84, GC =  51         OR :          Score = 35.67, E = 6.512e-06, P = 5.57e-12, GC = 40
    g.score= float(   self.last_line.split(',')[0].split()[-1]  )
    if " E = " in self.last_line:         g.evalue=e_v(     self.last_line.split(',')[1].split()[-1]   )
  
    self.last_line=self.file.readline()          
    self.last_line=self.file.readline()

    query_seq=''; target_seq=''
    while self.last_line.strip(): #two empty lines will make terminate this cycle. This is the signal for: end of this hit.
      #now in ss
      g.ss+=self.last_line.strip()
      self.last_line=self.file.readline()
      query_seq+= join(self.last_line.split()[1:-1], '')
      self.last_line=self.file.readline()
      self.last_line=self.file.readline()
      target_seq+=join(self.last_line.split()[1:-1], '')
      self.last_line=self.file.readline()
      self.last_line=self.file.readline()      

    while "*" in query_seq:
      pos_insert=query_seq.find('*')
      insert_length_target=  int( target_seq.split('*[')[1].split(']')[0].strip()  )
      insert_length_query= int( query_seq.split('*[')[1].split(']')[0].strip() )
      query_seq= query_seq.split('*')[0]+ insert_length_query*'x'+"-"* (insert_length_target-insert_length_query)   + join( query_seq.split('*')[2:], '*'  )
      target_seq=target_seq.split('*')[0]+insert_length_target*'x'+"-"*(insert_length_query-insert_length_target)   + join( target_seq.split('*')[2:], '*'  )
      l_ss_gap=0
      while len(g.ss)>pos_insert+l_ss_gap and    g.ss[pos_insert+l_ss_gap]=='~': l_ss_gap+=1
      g.ss = g.ss[:pos_insert]+'~'*max([insert_length_query, insert_length_target])+g.ss[pos_insert+l_ss_gap:]

    g.alignment.add('q', replace(query_seq, '.', '-'))
    g.alignment.add('t', replace(target_seq, '.', '-'))

    return g

class covelshit(gene):
  """Gene class to manage predictions by covels """
  def __init__(self, **kargs):
    self.score=None
    self.cm_file=None
    self.sequence_data=None
    gene.__init__(self, **kargs)

  def sequence(self):
    """ Returns the nucleotide sequence of this hit, obtained with lazy computing principle"""
    if self:
      if not self.sequence_data:        self.sequence_data = upper(self.fasta_sequence()[1])
      return self.sequence_data

  def summary(self):
    """ """
    if not self: return 'Empty covelshit'
    out='cm_file: '+str(self.cm_file)+'\n'
    out+='>'+self.chromosome+'\n\n'
    out+=' Strand = '+self.strand+'\n'
    out+=" Score = "+str(self.score)+'\n\n'
    
    if self.strand=='+': out+=str(self.boundaries()[0]).ljust(15)+' '
    else:                out+=str(self.boundaries()[1]).ljust(15)+' '
    out+=self.sequence()
    if self.strand=='+': out+=str(self.boundaries()[1])
    else:                out+=str(self.boundaries()[0])
    return out

class parse_covels(parser):
  """ Parse a covels-SE output covels objects at each next() call"""
  def load(self, filename=''):
    if not filename:
      filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file=open(filename, 'r')
    self.target=None
    self.cm_file=None
    self.last_line=self.file.readline()
    while self.last_line and not ' : ' in self.last_line:
      if self.last_line.startswith('Database to search/score'):   self.target= self.last_line.split()[-1]
      if self.last_line.startswith('Model'):                      self.cm_file=self.last_line.split()[-1]
      self.last_line=self.file.readline()

  def parse_next(self):
    while self.last_line and not ' : ' in self.last_line:                   self.last_line=self.file.readline()
    if self.last_line:
      g=covelshit(target=self.target, cm_file=self.cm_file)
      g.chromosome = self.last_line.split()[-1]
      g.score= float( self.last_line.split()[0] )
      start = int(  self.last_line.split()[1] )
      end   = int(  self.last_line.split()[2] )
      if start > end: 
        g.strand='-'
        start, end= end, start
      else:       g.strand='+'
      g.add_exon(start, end)
      self.last_line=self.file.readline()
      return g

class erpinhit(gene):
  """Gene class to manage predictions by covels """
  def __init__(self, **kargs):
    self.score=None
    self.evalue=None
    self.epn_file=None
    self.seq=None
    gene.__init__(self, **kargs)

  def summary(self):
    o= '  ERPIN hit -- ID: '+str(self.id)+'\n'
    o+='-Epn_file: '+str(self.epn_file)+'\n'
    o+='-Target: '+str(self.target)+'\n'
    o+='-Chromosome: '+str(self.chromosome)+' -Strand: '+str(self.strand)+' -Positions: '+str(self.boundaries()).strip('[').strip(']')+'\n'
    o+='-Score: '+str(self.score)+' -Evalue: '+str(self.evalue)+'\n'
    o+='-Seq: '+str(self.seq)
    return o
    

class parse_erpin(parser):
  def clean_from_service_msg(self, line):
    if '\r' in line and 'Kb:' in line:
      st=  line.find('Kb:')
      end= line.rfind('\r')
      line=line[:st]+line[end+1:]
    return line

  def load(self, filename=''):
    if not filename:
      filename=self.file.name
    check_file_presence(filename, 'filename')
    self.file=open(filename, 'r')
    self.target=None
    self.epn_file=None
    self.current_chromosome=None
    self.last_line=self.file.readline()
    while self.last_line and (len(self.last_line)<2 or ( not  (self.last_line[:2] in ['FW', 'RC'] or self.last_line.startswith('>')) )):      
      if self.last_line.startswith('Training set:'): self.epn_file=self.last_line.split('"')[1]
      elif self.last_line.startswith('Database:'): self.target=self.last_line.split('"')[1]
      self.last_line=self.clean_from_service_msg( self.file.readline() )

  def parse_next(self):
    while self.last_line and (len(self.last_line)<2 or ( not  (self.last_line[:2] in ['FW', 'RC'] or self.last_line.startswith('>')) )):       
      self.last_line=self.clean_from_service_msg( self.file.readline() )
    if self.last_line:
      if self.last_line.startswith('>'): 
        self.current_chromosome = self.last_line.split()[0][1:] #taking only the first word for consistency with other parsers. here we would have the complete fasta title, anyway
        self.last_line=self.clean_from_service_msg(  self.file.readline() )

      g=erpinhit( chromosome= self.current_chromosome, strand={'FW':'+', 'RC':'-'}[self.last_line.split()[0]], epn_file=self.epn_file, target=self.target )
      g.id = self.last_line.split()[1]
      g.score = float(self.last_line.split()[3])
      g.evalue= float(self.last_line.split()[4].strip())
      start_str, end_str= self.last_line.split()[2].split('..')
      g.add_exon(int(start_str), int(end_str))
      self.last_line=self.clean_from_service_msg( self.file.readline() )
      g.seq=self.last_line.strip()
      self.last_line=self.clean_from_service_msg( self.file.readline() )
      return g
    else: self.stop()
      
def secis_alignment(secis_list):
  """Input is a list of secis (gene) objects, which must have the primary sequence. The sequences must be split by white spaces in the ss components, as the output of Patscan does.
  The function aligns initially all these bits with each other, then realign the portions that may be improved 
  """
  pieces_lengths=[10, 8, 13, 5, 13, 2, 25, 14, 5, 10, 9, 10]  
  alignment_pieces=[alignment() for i in range(12)]
  for secis in secis_list:
    #print secis.seq
    splt=secis.seq.split()
    if len(splt)!=12:
      raise Exception, "secis_alignment ERROR sequence provided in the secis must be split in the 12 components by white spaces"
    for i in range(12):
      if pieces_lengths[i]-len(splt[i])<0: raise Exception, "secis_alignment ERROR pieces length is too short for piece: "+str(i)+' '+str(len(splt[i]))+'>'+str(pieces_lengths[i])
      alignment_pieces[i].add( secis.id, '-'*(pieces_lengths[i]-len(splt[i])) + splt[i] ) #completing seq to desired lenght using gaps
  final_ali=alignment()
  for secis in secis_list:
    seq=''
    for i in range(12):
      seq+=alignment_pieces[i].seq_of(secis.id) #+'Z'
    final_ali.add(secis.id, seq)

  pos, cols_list = 0, []
  for piece_l in pieces_lengths:
    cols_list.append([pos, piece_l])
    pos+=piece_l
  final_ali.realign_columns(input_list=cols_list)
  final_ali.remove_useless_gaps()
  final_ali.convert_sequences(upper)
  return final_ali
    
global sorted_blast_hits_temp_list
sorted_blast_hits_temp_list=[]  
class blasthit_list(list):
  """ Puppet class to handle the sorting of blast hits. This class can be initiated using another list as argument.
      The main and only method is sort. 
  """
  
  def __init__(self, inputlist=[]):
    self.features={}
    del self[:]
    for i in inputlist: self.append(i)      

  def sort(self, sorting_guide=['chromosome', 'strand', 'position' ],  inplace=True, toplevel=True): 
    """ sorting_guide contains the field to sort the list, ordered by priority. It can also be a function to apply to the self to return the key to be evaluted with a < statement.
        By default, sorting_guide is ['chromosome', 'strand', 'position' ], so the output is like this:  (genes on same chromosome clusters together, then inside those, genes on the same strand are clustered together, inside those cluster they are ordered by start position. This way of ordering facilitates merging.
        
          1	generic_program	cds	958795	959028	.	+	.	28215760 
          1	generic_program	cds	3678820	3679188	.	+	.	28215376 
          10	generic_program	cds	26931271	26931345	.	-	.	28214608 
          10	generic_program	cds	26939920	26939970	.	-	.	28214736 
          10	generic_program	cds	26996481	26996582	.	-	.	28214480 
          10	generic_program	cds	27008546	27008674	.	-	.	28214992 
          11	generic_program	cds	440442	440957	.	+	.	28216912 
          11	generic_program	cds	110088239	110088574	.	+	.	28217040 
          13	generic_program	cds	105985582	105986034	.	-	.	28217168 
          16	generic_program	cds	56119934	56120302	.	+	.	28971408 
          16	generic_program	cds	66120225	66120647	.	+	.	28971152 
          16	generic_program	cds	2527808	2528107	.	-	.	28971280 
    
     """
    global sorted_blast_hits_temp_list
    if sorting_guide:
      if toplevel:
        sorted_blast_hits_temp_list=[]
        
      key=str(sorting_guide[0])
      if len(sorting_guide)>1 :
        self.features[key]={}                        #self.features.chromosome={} ->contains possible values of chromosome, linked to the list of those. each list is then sorted

        for bh in self:
          try:
            value=bh.__dict__[key]             #chr1=bh.chromosome
          except:
            printerr('blast_list.sort ERROR can\'t obtain the field '+str(key)+' for blasthit '+str(bh), 1)
            raise

          if not self.features[key].has_key(value):    self.features[key][value]=blasthit_list()   #self.features.chromosome['chr1']= blastlist[]
          self.features[key][value].append(bh)       # --> so we have a list of every blast hit in the chromsome 1, a list for chromsome 2....        

        for value in sorted(self.features[key].keys()):
          self.features[key][value].sort(sorting_guide[1:], toplevel=False)  

      elif len(sorting_guide)==1 :      
        if type(key)==str and key=='position':
          list.sort(self, key=lambda bh : bh[0][0]) # NB blast hits must have filled (at least one exon)        
        elif type(key)==str:
          list.sort(self, key=lambda bh : bh.__dict__[key])
        elif type(key)==type(lambda a:a):
          list.sort(self, key=key)

        sorted_blast_hits_temp_list.extend(self)

      if toplevel:
        self.features.clear()
        if inplace:
          self.__init__(sorted_blast_hits_temp_list)
          return self
        return blasthit_list( sorted_blast_hits_temp_list )

class species():
  """ Will control connection with ncbi and tree information when I finish it """
  def __init__(self, name=''):
    self.name=''
    self.taxid=0
    if name: self.name=name
  def __str__(self): return self.name
  def fill_taxid(self, name='', silent=False):
    """ search the self.name or the name provided in ncbi taxonomy and assign a taxid. if none is found, raise a NameError exception. if more than one found, sets the one with shortest Scientific name as taxid, then raise a NameError exception with "WARNING" in the message. If silent == True, no exception is raised in any case.
     """
    search_string = name if name else self.name 
    search_results=ncbi_taxonomy.main({'silent':1, 'S':search_string, 'print_opt':1})
    def search_result_to_id(line): return int( line.split()[0] )
    def search_result_to_scientific_name(line): return  del_white(line[:-1].split('->')[1].split('#')[0])
    if not len(search_results.keys()): 
      if not silent: raise NameError, "ncbi_taxonomy ERROR can't find species: "+search_string  
    elif len(search_results.keys())> 1: 
      taxid=int(min(search_results.keys(), key=lambda x : search_result_to_scientific_name(search_results[x])))
      if not silent: raise NameError, "ncbi_taxonomy WARNING searching species '"+search_string+"' more than one species found: "+join([ search_result_to_scientific_name(s) for s in search_results ] , ',')
    self.taxid=search_result_to_id(search_results[search_results.keys()[0]])
    return True

  def __nonzero__(self):
    return bool(self.taxid)
    
  def load_from_db(self, db_object):
    """ Utility useful just in combination with the sqlite database of selenoprofiles. It loads  species object directly out of there. """
    self.taxid=db_object[0]
    self.name=db_object[1]
   
  def name_with_underscores(self, chars_not_allowed=':/;#@%[]()<>&_='):
    """ Returns the name of the species formatted in order to be used as a folder or file name. Forbidden characters are replaced with {chX}, where X is the numeric ASCII code for the char."""
    out=mask_characters(self.name, chars_not_allowed)
    return replace_chars(out, ' ', '_')

def center_str(s, n_char, fill=' '):
    o=s
    while len(o) < n_char:      o=fill+o+fill
    while len(o) > n_char:
      if len(o) % 2: o=o[:-1]
      else:          o=o[1:]
    return o

def ete_tree_correct_forbidden_chars(tree, forbidden_chars='():,', inplace=True):
  """ utility for ete to avoid the problem of newick conversion. In fact, in this format there are forbidden chars such as ( ) :   which cannot be present in node names. 
      This function detect the node with names including those characters and replace them to {chXX}, where XX is the ASCII integer identifier for the char.
      
      If inplace==True, the tree is modified in place, if it is not, a new tree with changed names is returned.
  """ 
  if not inplace: tree=deepcopy(tree)
  for node in tree.traverse():    node.name=mask_characters(node.name, forbidden_chars)
  if not inplace: return tree
  
  
def ete_tree_restore_forbidden_chars(tree, inplace=True):
  """ Reverse of the last function. Restore the names of a newly loaded tree changing the {chXX} occurences with the characters they mean """
  if not inplace: tree=deepcopy(tree)
  for node in tree.traverse():    node.name=unmask_characters(node.name)
  if not inplace: return tree


chars_to_replace_in_filenames='():,*./\'_=#[]'

def mask_characters(a_string, chars='DEFAULT'):
  """ This function is used to mask some "forbidden characters", provided as any iterable (list). These characters are replaced with {chXX}, where XX is the ASCII integer identifier for the char.  """
  if chars=='DEFAULT':  chars=chars_to_replace_in_filenames
  for char in chars: 
    if char in a_string:        a_string = replace_chars(a_string, char, '{ch'+str(ord(char))+'}')
  return a_string
  
def unmask_characters(a_string, replace_underscores=False):
  """ This is the inverse of the last function. """
  if replace_underscores: a_string=replace(a_string, '_', ' ')
  while '{ch' in a_string:
    try: 
      char_n=int(a_string.split('{ch')[1].split('}')[0])
      a_string=join(a_string.split('{ch'+str(char_n)+'}'), chr(char_n))
    except: pass
  return a_string

  
class pfamhit(gene):
  """ Class to manage pfamhits ... NOT FINISHED!! """
  def load_data( pfam_family=None, pfam_start=None, pfam_end=None, query_name=None, query_start=None, query_end=None, ali=None ):
    """ Utility to load all the data providing it direclty as arguments. ali must be an alignment of two sequences, the first is the query and the last is the target. the titles are not taken in to account, they are saved into self.alignment with titles 'q' and 't' respectively. """
    self.chromosome=pfam_family
    if pfam_start and pfam_end:      self.add_exon(pfam_start, pfam_end)
    self.strand='+'

    self.query=gene()
    self.query.chromosome=query_name
    
    if query_start and query_end:    self.query.add_exon(pfam_start, pfam_end)

    self.alignment=alignment()
    if ali and ali.nseq()!=2: raise Exception, "pfamhit->load ERROR alignment provided can have only two sequences... it has: "+str(ali.nseq())
    self.alignment.add('q', ali.titles()[0])
    self.alignment.add('t', ali.titles()[1])
    
    

def fasta_next_seq(filehandler, cline=''):
  """This function is thought to parse a fasta file and return a sequence at the time, to process each one without laoding the whole file.
  The input is a file handler, and a string for the line which the program may need to parse more than once...
  The use is:

  f=fasta_next_seq(  seqs_file_h, '')
  while f:
    title, seq, cline= f
    ## perform operations on title, seq
    f=fasta_next_seq(  seqs_file_h, cline)

    !!!! obsolete! use parse_fasta() class instead
   """

  title, seq='', ''
  if not cline:
    cline=filehandler.readline()
  while cline:
    if cline[0]=='>':
      if title:
        return title, seq, cline
      title=cline[1:-1]
    else:
      seq+=cline[:-1]
    cline=filehandler.readline()
  if title:
    return title, seq, cline
  else:
    return None
def check_file_presence(input_file, descriptor='input_file', exception_raised=Exception):
  if not input_file or not is_file(input_file):
    raise exception_raised, "ERROR "+descriptor+ ": "+str(input_file) +" not defined or not found. Run with option -h for help."
def check_directory_presence(input_file, descriptor='folder', exception_raised=Exception):
  if not input_file or not is_directory(input_file):
    raise exception_raised, "ERROR "+descriptor+ ": "+str(input_file) +" not defined or not found. Run with option -h for help."

def test_writeable_folder(folder, descriptor=''):
  rnumber= random.randint(1, 99999)
  folder=Folder(folder)
  filename=folder+'WrItE_TeSt.'+str(rnumber)
  if bash('echo "x" > '+filename+' && rm '+filename)[0]:
    raise Exception, "ERROR "+descriptor+ ": cannot write in "+folder

is_file=os.path.isfile
is_directory=os.path.isdir
abspath=os.path.abspath
base_filename=os.path.basename
def directory_name(*args, **kargs):
  out=os.path.dirname(*args, **kargs)
  if out=='': return '.'
  return out
file_size=os.path.getsize

def list_of_lines(inputfile):
  """ Return the list of lines in file: inputfile, removing the newline characterss \\n """
  check_file_presence(inputfile)
  out=[line.strip() for line in open(inputfile, 'r')]
  if out==['']: return []
  return out
  
def sankoff(tree, node2seqFn=None, matrix=None):
  """ Sankoff algorithm for the reconstrunction of ancestral states. tree is a ete2.Tree instance, node2seqFn is a function (leafnode)-> its seq, matrix is a hash with substitution costs (should have zeros in diagonal).
  NOTE: matrix.keys() is used to get all possible letter in the alphabet considered
  letters not present in the alphabet are considered equally (im)possibile. This have the effect of excluding those nodes for the calculation of ancestral states at that position. This makes possible its use for gapped alignments.
  A node2ancestral_sequence hash is returned
  """
  if matrix is None:
    matrix={'A':{'A':0, 'C':5, 'G':2, 'T':5},            'C':{'A':5, 'C':0, 'G':5, 'T':2},            'G':{'A':2, 'C':5, 'G':0, 'T':5},            'T':{'A':5, 'C':2, 'G':5, 'T':0}    }
  if node2seqFn is None:     node2seqFn= lambda x:x.sequence
  seq_length= len(node2seqFn( tree.get_leaves()[0] )) #testing node2seqFn 
  alphabet=matrix.keys()
  node2ancestral={}
  for sequence_index in range(seq_length):
    ### assigning costs, leaf to root
    node2cost={};   node2costPerNt={}; 
    for node in  tree.traverse(strategy='postorder'):
      if node.is_leaf():
        lett_this_node=node2seqFn(node)[sequence_index]
        node2cost[node]=0
        node2costPerNt[node]={}
        for lett in alphabet:
          if lett==lett_this_node:          node2costPerNt[node][lett]= 0
          else:                             node2costPerNt[node][lett]= sys.maxint-1000
      else: 
        #initializing ancestral seq on first round (index=0)
        if not node2ancestral.has_key(node): node2ancestral[node]=''
        node2costPerNt[node]={}
        for lett in alphabet:
          c=0
          for child in node.get_children():
            #best nt:
            best_lett2=None
            for lett2 in alphabet:
              if node2costPerNt[child][lett2]==node2cost[child]:
                best_lett2=lett2
                break
            if best_lett2 is None: cost_change=sys.maxint-1000          #special case
            else:                  cost_change= matrix[lett][best_lett2] 
            c+=min( [node2cost[child]+   cost_change   , node2costPerNt[child][lett] ]  )  # cost change, cost unchange
          node2costPerNt[node][lett]=c
        node2cost[node] = min([ node2costPerNt[node][lett]  for lett in alphabet  ])

    ## now backtracking, assigning ancestral states (root to leaves)
    for node in tree.traverse(strategy='preorder'):
      if node==tree: #root 
        for lett in alphabet: 
          if node2costPerNt[node][lett]==node2cost[node]:
            node2ancestral[node]+=lett
            break
      elif not node.is_leaf():
########################## !!! ########################## ---> triple check this ONE: isn't it matrix[something] ? 
        if node2cost[node]+1 > node2costPerNt[node][  node2ancestral[node.up][sequence_index] ]:        node2ancestral[node]+=node2ancestral[node.up] [sequence_index]
        else: 
          for lett in alphabet: 
            if node2costPerNt[node][lett]==node2cost[node]:
              node2ancestral[node]+=lett
              break
  return node2ancestral


codon2sitecount={}  #{ codon: [sites]  }  #sites as with split_nonsense=True
def count_sites(cds, silent=False, split_nonsense=False):
  """Counts the number of possible  Syn and nonsyn sites for a input nucletoide coding sequence. 
  As a single site can be partly non-syn and partly syn, the numbers returned are float (always mutiple of one third)
  The functino computes the number also separately for CpG sites.   To obtain the number of nonCpG changes, subtract the number of CpG sites from the total number
  Codons with any character different from ACTG (for example, N) are skipped and a message is printed to stderr
  Returns [nonSyn, Syn,  CpG_nonSyn, CpG_syn]
  if nonsense==True, nonsense (stop) mutations are differentiated from nonsyn mutations. the function instead returns [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ] 
  """   
  global codon2sitecount
  cds=replace(upper(nogap(cds)), 'U', 'T')
  syn=0    ;  nonsyn=0   # these will result to be three times as much the actual values: I dive them as the very last step!
  cpg_syn=0 ; cpg_nonsyn=0
  nonsense=0;   cpg_nonsense=0
  #noncpg_syn=0 ; noncpg_nonsyn=0

  if len(cds)%3!=0: raise Exception, "count_sites ERROR the sequence must be composed of codons (length multiple of 3)"
  for i_codon in range(len(cds)/3):
    ## cycling codons
    codon=cds[i_codon*3:i_codon*3+3]
    if not all([lett in 'ACTG' for lett in codon] ):       
      if not silent:        printerr('count sites WARNING skipping codon n.'+str(i_codon+1)+' : '+codon, 1)
      continue

    if codon in codon2sitecount:      n,s,x,cn,cs,cx=codon2sitecount[codon]
    else:
      n,s,x,cn,cs,cx=0,0,0,0,0,0

      for i_within_codon in range(3):
        ## cycling each position
        nt= codon[i_within_codon]
        i_cds=i_codon*3+i_within_codon
        syn_changes_this_pos=0; nonsense_this_pos=0
        for alt_nt in 'ACTG':
          if alt_nt==nt: continue
          alt_codon=   codon[:i_within_codon]+alt_nt+codon[i_within_codon+1:]
          if transl(alt_codon)==transl(codon):          
            syn_changes_this_pos+=1
          elif "*" in transl(alt_codon) +transl(codon): nonsense_this_pos+=1

        is_cpg= (   nt == 'G' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='C' )  or (i_cds!=0 and cds[i_cds-1] =='C' )    )     ) or \
              (   nt == 'C' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='G' )  or (i_cds!=0 and cds[i_cds-1] =='G' )    )     )        # # G and the next or previous is C OR #C and the

        s+=syn_changes_this_pos      
        n+=(3-syn_changes_this_pos-nonsense_this_pos)
        x+=nonsense_this_pos
        #syn+=    s
        #nonsyn+= n 
        #nonsense+=nonsense_this_pos
        if is_cpg:    
          cs+=  syn_changes_this_pos
          cn+= (3-syn_changes_this_pos-nonsense_this_pos)
          cx+= nonsense_this_pos
          # cpg_syn+=  syn_changes_this_pos
          # cpg_nonsyn+= (3-syn_changes_this_pos-nonsense_this_pos)
          # cpg_nonsense+= nonsense_this_pos
      codon2sitecount[codon]=n,s,x,cn,cs,cx
    syn+=s
    nonsyn+=n
    nonsense+=x
    cpg_syn+=cs
    cpg_nonsyn+=cn
    cpg_nonsense+=cx

      #else:
      #  noncpg_syn+=  syn_changes_this_pos
      #  noncpg_nonsyn+= (3-syn_changes_this_pos)        
  #write('Length: '+str(len(cds)), 1)
  #write('Syn: '+str(syn/float(3))+' Nonsyn: '+str(nonsyn/float(3))+' ', 1   )  
  #write('CpGSyn: '+str(cpg_syn/float(3))+' CpGNonsyn: '+str(cpg_nonsyn/float(3))+' ', 1   )  
  #write('NonCpGSyn: '+str(noncpg_syn/float(3))+' NonCpGNonsyn: '+str(noncpg_nonsyn/float(3))+' ', 1   )  
  if split_nonsense:    return [ (nonsyn)/float(3), syn/float(3), nonsense/float(3), (cpg_nonsyn)/float(3),  cpg_syn/float(3), (cpg_nonsense)/float(3)]
  else:                 return [ (nonsyn+nonsense)/float(3), syn/float(3), (cpg_nonsyn+cpg_nonsense)/float(3),  cpg_syn/float(3)]

def count_changes(cds, cds2, silent=True, split_nonsense=False):
  """ count the number of Syn and non.syn between two sequences.
  The function computes the number also separately for CpG sites.   To obtain the number of nonCpG changes, subtract the number of CpG sites from the total number.
  NOTE: for CpG call, the function looks only at the first sequence provided. So for this, results are not completely symmetrical 
  positions with a gap in any of the two  sequences are skipped
  returns [non_syn,  syn, cpg_non_syn, cpg_syn ] 
 if nonsense==True, nonsense (stop) mutations are differentiated from nonsyn mutations. the function instead returns [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ] 

  """
  cds=replace(upper(cds), 'U', 'T')
  cds2=replace(upper(cds2), 'U', 'T')

  syn=0    ;  nonsyn=0   ;  nonsense=0
  cpg_syn=0 ; cpg_nonsyn=0; cpg_nonsense=0
  #noncpg_syn=0 ; noncpg_nonsyn=0
  
  if len(cds)%3!=0 or len(cds2)%3!=0: raise Exception, "count_changes ERROR the sequences must be composed of codons (length multiple of 3)"
  if len(cds)!=len(cds2):             raise Exception, "count_changes ERROR the sequences do not have the same length"
  for i_codon in range(len(cds)/3):
    ## cycling codons
    codon=cds[i_codon*3:i_codon*3+3]
    if not all([lett in 'ACTG' for lett in codon] ):       
      if not silent:        printerr('count_changes WARNING skipping codon n.'+str(i_codon+1)+' : '+codon, 1)
      continue
    codon2=cds2[i_codon*3:i_codon*3+3]
    if not all([lett in 'ACTG' for lett in codon2] ):       
      if not silent:        printerr('count_changes WARNING skipping codon2 n.'+str(i_codon+1)+' : '+codon2, 1)
      continue
    if ('-' in codon and codon!='---') or  ('-' in codon2 and codon2!='---') : raise Exception, "count_changes ERROR the sequences must be aligned by codon, i.e. the gaps must be in groups of three"
    if '-' in codon or '-' in codon2: 
      #skipping gap position
      continue

    if codon!=codon2:
      #change! (codon)
      if transl(codon)!=transl(codon2):          
        if "*" in transl(codon) +transl(codon2): is_nonsyn=2 #nonsense mutation
        else:                                    is_nonsyn=1
      else:                                      is_nonsyn=False
      for i_within_codon in range(3):
        if codon[i_within_codon]!=codon2[i_within_codon]:
          #change! (nt)
          i_cds=i_codon*3+i_within_codon
          nt=cds[i_cds]
          is_cpg= (   nt == 'G' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='C' )  or (i_cds!=0 and cds[i_cds-1] =='C' )    )     ) or \
            (   nt == 'C' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='G' )  or (i_cds!=0 and cds[i_cds-1] =='G' )    )     )        # # G and the next or previous is C OR #C and the

          if not is_nonsyn:    syn+=1
          elif is_nonsyn==2:   nonsense+=1
          else:                nonsyn+=1                
          if is_cpg: 
            if not is_nonsyn:    cpg_syn+=1
            elif is_nonsyn==2:   cpg_nonsense+=1
            else:                cpg_nonsyn+=1                

  if split_nonsense:    return [nonsyn, syn, nonsense,  cpg_nonsyn, cpg_syn, cpg_nonsense]
  else:                 return [nonsyn+nonsense, syn, cpg_nonsyn+cpg_nonsense, cpg_syn]


def count_unique_changes(cds, other_cds_list, silent=True, split_nonsense=False):
  """ This function is an extension of count_changes, but counts changes in several sequences at the time, and filters out those common to more than a sequence. The non-uniqness is computed just based on the sequences and no tree is inspected. 
  returns [non_syn,  syn, cpg_non_syn, cpg_syn ]  
 if nonsense==True, nonsense (stop) mutations are differentiated from nonsyn mutations. the function instead returns [ nonSyn, Syn, NonSense,  CpG_nonSyn, CpG_syn, CpG_nonsense ] 

  """
  cds=replace(upper(cds), 'U', 'T')

  position_to_change={}  ### hash keeping track of changes we already saw   k: position (0based) -> list of other nts observed in any of other_cds
  
  syn=0    ;  nonsyn=0   ;   nonsense=0
  cpg_syn=0 ; cpg_nonsyn=0 ; cpg_nonsense=0
  
  for cds2 in other_cds_list:
    cds2=replace(upper(cds2), 'U', 'T')

    if len(cds)%3!=0 or len(cds2)%3!=0: raise Exception, "count_unique_changes ERROR the sequences must be composed of codons (length multiple of 3)"
    if len(cds)!=len(cds2):             raise Exception, "count_changes ERROR the sequences do not have the same length"
    for i_codon in range(len(cds)/3):
      ## cycling codons
      codon=cds[i_codon*3:i_codon*3+3]
      if not all([lett in 'ACTG' for lett in codon] ):       
        if not silent:        printerr('count_unique_changes WARNING skipping codon n.'+str(i_codon+1)+' : '+codon, 1)
        continue
      codon2=cds2[i_codon*3:i_codon*3+3]
      if not all([lett in 'ACTG' for lett in codon2] ):       
        if not silent:        printerr('count_unique_changes WARNING skipping codon2 n.'+str(i_codon+1)+' : '+codon2, 1)
        continue
      if ('-' in codon and codon!='---') or  ('-' in codon2 and codon2!='---') : raise Exception, "count_unique_changes ERROR the sequences must be aligned by codon, i.e. the gaps must be in groups of three"
      if '-' in codon or '-' in codon2: 
        #skipping gap position
        continue

      if codon!=codon2:
        #change! (codon)
        if transl(codon)!=transl(codon2):          
          if "*" in transl(codon) +transl(codon2): is_nonsyn=2 #nonsense mutation
          else:                                    is_nonsyn=1
        else:                                      is_nonsyn=False
        for i_within_codon in range(3):
          if codon[i_within_codon]!=codon2[i_within_codon]:
            #change! (nt)
            i_cds=i_codon*3+i_within_codon
            nt=cds[i_cds];    nt2=cds2[i_cds]
            #### dtermining if the change is uniq
            if position_to_change.has_key( i_cds ) and nt2 in position_to_change[i_cds]: 
              continue  #skipping non-uniq
            #is uniq
            if not position_to_change.has_key( i_cds ): position_to_change[i_cds]=[]
            position_to_change[i_cds].append(nt2)
                        
            is_cpg= (   nt == 'G' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='C' )  or (i_cds!=0 and cds[i_cds-1] =='C' )    )     ) or \
              (   nt == 'C' and (  (  i_cds+1<len(cds)  and   cds[i_cds+1] =='G' )  or (i_cds!=0 and cds[i_cds-1] =='G' )    )     )        # # G and the next or previous is C OR #C and the
            if not is_nonsyn:    syn+=1
            elif is_nonsyn==2:   nonsense+=1
            else:                nonsyn+=1                
            if is_cpg: 
              if not is_nonsyn:    cpg_syn+=1
              elif is_nonsyn==2:   cpg_nonsense+=1
              else:                cpg_nonsyn+=1                
              
  if split_nonsense:    return [nonsyn, syn, nonsense,  cpg_nonsyn, cpg_syn, cpg_nonsense]
  else:                 return [nonsyn+nonsense, syn, cpg_nonsyn+cpg_nonsense, cpg_syn]


def add_ancestral_states(t, rst_file, add=True):
  """ This function must be called after a codeml analysis with RateAncestor=1 has been run using a certain tree t. this function parses the rst_file and populates the non-leaf node and adds them a .sequence attribute. A .index attribute is also added to all nodes (useful particularly for non-leaf nodes).
  If you don't want to modify the tree t, specify add=False. in this case, a hash would be returned with k: node -> value:  sequence 
  """
  try: PhyloTree()
  except: 
    from ete2 import PhyloTree
  fh=open(rst_file, 'r')
  line=fh.readline()
  while line and line != "tree with node labels for Rod Page's TreeView\n":    line=fh.readline()
  if not line: raise IOError, "add_ancestral_states ERROR can't find line \"tree with node labels for Rod Page's TreeView\" in file: "+rst_file
  line=fh.readline() #now in next line: tree with nodes names like:       (((((1_Homo_sapiens, 3_Pan_troglodytes) 11 , 2_Gorilla_gorilla) 10 , 4_Pongo_pygmaeus) 9 , 5_Macaca_mulatta) 8 , 6_Callithrix_jacchus) 7 ;
  tree_text=replace(line, ') ', '):')  #making it read like it was the distance
  t2=PhyloTree(tree_text)
  for node in t2.traverse():
    if node.is_leaf():      
      node.index= int(node.name.split('_')[0])
      node.name= join( node.name.split('_')[1:], '_' )
    else: 
      node.name= str(int(node.dist))
      node.index=    int(node.dist)

  while line and line != "List of extant and reconstructed sequences\n":    line=fh.readline()  
  line=fh.readline()  ;   line=fh.readline()  ;   line=fh.readline()  ;   line=fh.readline()   
  #now in line of first node
  while line.strip():
    if line.startswith('node'): #ancestral node, reconstructed
      index=int(line.split()[1][1:])
      #print 'searching index '+str(index)
      (t2&(str(index))).sequence=   join(line.split()[2:], '')
      #print 'adding to node '+str(index)+' '+(t2&(str(index))).seq
    line=fh.readline()   
  fh.close()

  if add==True:
    for node_t1, node_t2 in mapping_trees(t, t2):
      node_t1.index=node_t2.index
      if not node_t1.is_leaf() :     node_t1.sequence=node_t2.sequence #and hasattr(node_t2, 'sequence')
  else:
    o={}
    for node_t1, node_t2 in mapping_trees(t, t2):
      if not node_t1.is_leaf():         o[node_t1]=node_t2.sequence    
    return o

def mapping_trees(t1, t2):
  """ Given two identical trees, it maps the nodes of one into the other, exploiting the leaves under each node -- a uniq property of trees (right?).  The algorithm is not very efficient, it uses almost brute force"""
  hash_two2one={}
  for node1 in t1.traverse():
    id1=   join(sorted(node1.get_leaf_names()), '&?')
    for node2 in t2.traverse():
      if not hash_two2one.has_key(node2):
        id2=   join(sorted(node2.get_leaf_names()), '&?')
        if id2==id1:        hash_two2one[node2]=node1
  out=[]
  for node2 in hash_two2one:
    out.append( [hash_two2one[node2], node2] )
  return out

def color_scale(percent, c1='#000000', c2='#FFFFFF' ):
  """ This function is useful when using colors to represent some value. Given two colors and a normalized value between 0 and 1, it returns a color between c1 and c2. Extreme case: value0 --> returns c1, value1 -> returns c2. the colors must be provided as RGB hex-string, e.g.   #008800   OR  008800 ; it is returned in the form #008800 """
  def dec2hex(n):
      """return the hexadecimal string representation of integer n"""
      return "%X" % n
  def hex2dec(s):
      """return the integer value of a hexadecimal string s"""
      return int(s, 16)
  c1=c1.strip('#');      r1, g1, b1=     hex2dec( c1[0:2] ), hex2dec( c1[2:4] ), hex2dec( c1[4:] ), 
  c2=c2.strip('#');      r2, g2, b2=     hex2dec( c2[0:2] ), hex2dec( c2[2:4] ), hex2dec( c2[4:] ), 
  r3= r1+ (r2-r1)*float(percent)
  g3= g1+ (g2-g1)*float(percent)
  b3= b1+ (b2-b1)*float(percent)
  hr3=dec2hex(r3).rjust(2, '0')
  hg3=dec2hex(g3).rjust(2, '0')
  hb3=dec2hex(b3).rjust(2, '0')
  return '#'+hr3+hg3+hb3

def color_scale_midpoint(percent, c1, c2, mid="#FFFFFF"):
  """ same concept of color_scale, but uses a midpoint color. so if the value is below 0.5, the midcolor between c1 and mid is returned, otherwise the midcolor betwee c2 and mid is returned.   """
  def dec2hex(n):
      """return the hexadecimal string representation of integer n"""
      return "%X" % n
  def hex2dec(s):
      """return the integer value of a hexadecimal string s"""
      return int(s, 16)
  midpoint = 0.5
  c2 = c2.strip('#');      r1, g1, b1 =      hex2dec( c2[0:2] ), hex2dec( c2[2:4] ), hex2dec( c2[4:] ), 
  c1 = c1.strip('#');      r2, g2, b2 =      hex2dec( c1[0:2] ), hex2dec( c1[2:4] ), hex2dec( c1[4:] ),
  midcolor =   mid.strip('#'); rm, gm, bm =  hex2dec( midcolor[0:2] ), hex2dec( midcolor[2:4] ), hex2dec( midcolor[4:] ),
  
  if percent >= midpoint:
    r3 = rm+ (r2-rm)*float(percent - midpoint) * 2
    g3 = gm+ (g2-gm)*float(percent - midpoint) * 2
    b3 = bm+ (b2-bm)*float(percent - midpoint) * 2  
  elif percent < midpoint:
    r3 = rm+ (r1-rm)*float(percent) * 2
    g3 = gm+ (g1-gm)*float(percent) * 2
    b3 = bm+ (b1-bm)*float(percent) * 2
    
  hr3 = dec2hex(r3).rjust(2, '0')
  hg3 = dec2hex(g3).rjust(2, '0')
  hb3 = dec2hex(b3).rjust(2, '0')
  return '#'+hr3+hg3+hb3


def LRT(lnL_M0,lnL_M1,np_M0, np_M1   ):
  """ returns the p-value of a LRT test using approximate chi2.
  inputs are: lnL_M0,lnL_M1,np_M0, np_M1     where M0 stands for null mode, M1 for alternative model,    lnL is the logarithm of the likelihood, and np the number of parameters
"""
  try: chi2
  except NameError: from scipy.stats import chi2
  D= 2*(lnL_M1 - lnL_M0 )
  ddf=  np_M1 - np_M0
  p_value = 1 - chi2.cdf(D, ddf)    
  return p_value

def tagged_tree2codeml(t, format=9, tag="tag"):
  """This is to get newick representation of a tree that can be input to codeml. One or more species are tagged with numeric tags.
  To define these, you need to add features to the tree before running this function:
  node= t&"whatever"
  node.add_feature('tag', '1')
  
  then you can use this function to get something codeml friendly:
  print tagged_tree2codeml(  t  )
  
  --> 
  Normally the ete format is 9, which is, without branch support or distances. The tag name is "tag"
  """
  s=t.write(format=format, features=[tag])
  def repl_function(m):
    return "#"+m.group(1)
  return re.sub(r"\[&&NHX:"+tag+"=(\d+)\]", repl_function, s)

def align_coding_sequences(titles_seqs_list, protein_alignment=False):
  """Aligns a set of coding sequences according to their peptide sequences. Normally, mafft is invoked to align residues. Otherwise, the alignment of their translation can be provided (alignemtn instance) as argument. 
  titles_seqs_list is a list of elements like [title, seq]
Returns an alignment instance. NOTE: it requires  temp_folder to be defined if protein_alignment is not provided --> the realign function is used """
  if not protein_alignment:
    protein_alignment=alignment()
    for t, s in titles_seqs_list:
      protein_alignment.add(t, transl(s))   
    protein_alignment.realign()

  if sorted(protein_alignment.titles()) != sorted([t for t, s in titles_seqs_list ]): 
    for i in  sorted(protein_alignment.titles()):
      print i
    print "--------------------------------"
    for i in   sorted([t for t, s in titles_seqs_list ]):
      print i
    raise Exception, "align_coding_sequences ERROR the titles provided and those in the protein_alignment do not correspond!"

  cds_alignment=alignment()
  for t, s in titles_seqs_list:
    gaps=0; aligned_cds=''
    for p in range(protein_alignment.length()):
      if protein_alignment.seq_of(t)[p]=='-':
        aligned_cds+='---'
        gaps+=1
      else:
        aligned_cds+=   s[(p-gaps)*3:(p-gaps)*3+3]
    cds_alignment.add(t, aligned_cds)
  return cds_alignment


def function_help(f):
  """ returns a ipython style doc for a certain function, describing arguments, defaults, and the __doc__ string"""
  import inspect 
  
  h=str(f.__name__)+'('
  args, varargs, varkw, defaults = inspect.getargspec(f)
  if defaults is None: n_defaults=0
  else:                n_defaults=len(defaults)
  for i in range(len(args) - n_defaults  ): #non defaulted args
    h+=args[i]+', '  
  for i in range(n_defaults  ): #non defaulted args
    value=defaults[i]
    if type(value)==str: value='"'+value+'"'
    else: value=str(value)
    h+=args[i+len(args) - n_defaults]+'=' +value+', '
  h=h[:-2]    +')\n\n'
  h+=str(f.__doc__)+'\n'  
  return h

  
def interactive_mode(vars=None, message="welcome to the shell" ):
  """ To open an interactive shell inside a python script. Usage: interactive_mode()() ; double parenthesis is because this returns a pointer to a function. """
  #prompt_message = "Welcome!  Useful: G is the graph, DB, C"
  prompt_message = message
  try:
      from IPython.Shell import IPShellEmbed
      ipshell = IPShellEmbed(argv=[''],banner=prompt_message,exit_msg="Goodbye")
      return  ipshell
  except ImportError:
      if vars is None:  vars=globals()
      import code
      import rlcompleter
      import readline
      readline.parse_and_bind("tab: complete")
      # calling this with globals ensures we can see the environment
      print prompt_message
      shell = code.InteractiveConsole(vars)
      return shell.interact


def load_chromosome_lengths(chromosome_length_file, max_chars=0, exception_raised=Exception):
  """Utility to load chromosome lenghts from a fastalength output file and also set it as a MMlib variable; also performing controls on the file """
  global chromosome_lengths; chromosome_lengths={}
  for line in open(chromosome_length_file, 'r'):     
    fasta_identifier = line.split()[1]
    length=int(line.split()[0])
    if chromosome_lengths.has_key(fasta_identifier): 
      bash('rm '+chromosome_length_file)
      raise exception_raised, "ERROR the target file has a duplicate fasta identifier! ("+line.split()[1]+') Please modify it and rerun. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if length==0: 
      bash('rm '+chromosome_length_file)
      raise exception_raised, "ERROR the target file has a length zero entry! ("+line.split()[1]+') Please modify it and rerun. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if is_number(fasta_identifier) and fasta_identifier[0]=='0':
      bash('rm '+chromosome_length_file)
      raise exception_raised, "ERROR the target file has a numeric fasta identifier starting with zero!  ("+line.split()[1]+') This would cause an unexpected blast behavior. Please modify this or these ids and rerun. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if ':subseq(' in fasta_identifier: 
      bash('rm '+chromosome_length_file)
      raise  exception_raised, "ERROR with fasta header: "+fasta_identifier+' ; this was generated by fastasubseq and will cause unexpected behavior of this program, since it is using fastasubseq itself to cut sequences. Please clean the titles in your target file from ":subseq(" tags. Note: remove the .index and *.fa.n* blast formatting files after changing the target file '
    if max_chars and   len(fasta_identifier)>max_chars: 
      bash('rm '+chromosome_length_file)
      raise  exception_raised, "ERROR with fasta header: "+fasta_identifier+' is too long. The maximum length for a fasta identifier (first word of the title) is '+str(max_chars)+' characters. Please clean the titles in your target file. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'
    if '/' in fasta_identifier:
      bash('rm '+chromosome_length_file)
      raise  exception_raised, "ERROR with fasta header: "+fasta_identifier+' has forbidden character: "/" \nPlease clean the titles in your target file. Note: remove the .index and *.fa.n* blast formatting files after changing the target file'

    chromosome_lengths[fasta_identifier]=length
  set_MMlib_var('chromosome_lengths', chromosome_lengths)

def get_chromosome_lengths():
  global chromosome_lengths
  return chromosome_lengths

def RNAplot(seq, ss, fileout, label='', options=''):
  """Produce a structure plot with RNAplot into fileout. 
  Standard postscript file is converted to specified format (from extension) if not .ps  (imagemagik suite required) 
  Label, if provided, adds a text label below the structure
"""
  global temp_folder
  bbash( 'cd {temp}; echo ">title\n{seq}\n{ss}" | RNAplot {opts} '.format(temp=temp_folder, seq=seq, ss=ss, opts=options) )
  expected_file=temp_folder+'title_ss.ps'
  output_extension = fileout.split('.')[-1]
  if output_extension=='ps' and not label:
    bbash('mv {ps} {out}'.format(ps=expected_file, out=fileout))
  else:
    labelbit='' if not label else '-label "{lab}"'.format(lab=label)
    bbash( 'montage -geometry +1+1 -density 150 {labelbit} {ps} {out}'.format(ps=expected_file, out=fileout, labelbit=labelbit) )

def RNAfold(seq, constraints=None, img=None, options='', title=None, rnaplot_options=''):
  """Runs RNAfold to predict the secondary structure of this sequence. 
  Optionally can produce an image (multiple extensions accepted);
  It can accept fold constraints; possible characters:   (.)    ;
 Returns:  [secondary structure, free energy]
  Secondary structure reported will have the length of the input sequence (without gaps, if any)"""
  seq=nogap(seq)
  if constraints: 
    accepted_chars={'|':0, 'x':0, '<':0, '>':0, '(':0, ')':0, '.':0}
    if not len(constraints)==len(seq): 
      raise Exception, 'RNAfold ERROR lenght of sequence and of constraints must be the same!  Seq.length: {s}  Constraints: {c}'.format(s=len(seq), c=len(constraints))
    if not all([c in accepted_chars  for c in constraints] ): 
      raise Exception, 'RNAfold ERROR illegal characters in constraints. Provided: {c}'.format(c=constraints)

  add='' if not constraints else '\n'+constraints
  if constraints: options+=' -C '
  rnafold_out=bbash( 'echo ">title\n{seq}{add}" | RNAfold --noPS {opts} '.format(temp=temp_folder, seq=seq, add=add, opts=options) )
  free_energy=    float(rnafold_out.split('(')[-1][:-1] )
  ss=  rnafold_out.split('\n')[-1].split( ' (')[0]
  label='{tit}E= {e}'.format(tit=title+'\n' if title else '', e=free_energy)
  if not img is None: RNAplot(seq, ss, fileout=img, label=label, options=rnaplot_options)
  return [free_energy, ss]

