#!/usr/bin/python -u
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
from selenoprofiles_3 import p2ghit, blasthit, exoneratehit, genewisehit

try:  
  from PyQt4 import QtCore, QtGui
  from ete2 import *
except: 
  printerr('ERROR PyQt4 and ete2 must be installed to use this program!', 1)
  raise 


help_msg="""This program serves to draw the distribution of selenoprofiles results on a species tree. It uses ete2.

usage:

selenoprofiles_tree_drawer.py    alignment1 alignment2 .. alignmentN   -t tree_topology.newick  [options]

Alignments must be .ali files as produced by selenoprofiles_join_alignments. In every alignment, only the sequences with a title matching that of a typical selenoprofiles output are considered. The species to which they belong is derived from the title.
Every such species must be present in the tree, or an exception is raised and the program crashes (unless option -g is active).  Alteratively, the tree file can be omitted. In this case the output is drawed on an unstructured tree with only leafs, with the names of the species found in the alignments. The tree can contain also species which do not contain any result. In this case, those nodes normally are not drawn (unless option -e is active)

### Options:
-m                        the input tree has masked species names; e.g. Xenopus (Silurana) tropicalis is coded as Xenopus_{ch40}Silurana{ch41}_tropicalis
-e                        empty nodes (without any assigned result) are also drawed
-o                        outgroup species. Without arguments, the midpoint species is computed and used.
-g                        results which are not found in the tree are tolerated, and are just ignored (the program does not crash).  a warning is printed for each result not found
-C                        circular tree
-common                   convert species names from scientific to common (only for species for which it's hardcoded -- see common_names hash in the script code)

-a                        suppress normal output; instead of images, colored numbers are used to summarize the number of hits for each family
-c                        Only if option -a is active, it compress the number of boxes summing up all those whose labels are not among some allowed ones, defined by option -explicit_labels arg1,arg2,... 
# only if option -a is NOT active:
-W                        gene brick width in pixels (default: 150)
-H                        gene brick height in pixels (default: 40)
-T                        don't display text on gene bricks (chromosome names and positions)
-I                        don't display intron markers 
-no_id                    don't display numeric id on the left
-f                        display colored vertical lines for features indicated in the fasta headers. Argument to this option must be of the form  feat1,feat2,feat3 ; each of these will be searched in the fasta headers, looking for something like (in case of feat1):  " feat1:1,33,41 " ; the numbers are understood as positions along the aminoacid sequence of this gene, in which this feature is present. Normally the features are displayed as black lines, but one can define each color adding an RGB color after the feature name, e.g  feat1#33FF44,feat2#EE0000
-add                      provide an tab delimited file to display additional features for the results. This file must have an entry per line, like protein_idTABtext (or type2)
-sp_add                   add annotation for each species. A tab delimited file must be provided with species TAB value. The values will be appended right next to the species name (aligned). The special value COLOR#xxxxx can be provided to change the color of the species name to xxxxxx (RGB, hex)

-out                      save image in the output file provided as argument (allowed extensions: pdf, png, svg) -- it does not open ete interactive environment.
-img_w                    image width, used for -out option
-img_h                    image height, used for -out option

The following options control the size of text drawn. If a text doesn't not fit its dedicated space, its point size is reduced until it does.
-tsize                    size of text for titles
-ssize                    size of text for species 
-bsize                    size of body text
-nsize                    size of numbers in the colored boxes (when option -a is active)

-prompt                   open interactive prompt at the end
-F                        print family names in each species line, not just as header of their columns 
-temp  +                  temporary folder. A folder with random name is created here, used and deleted at the end of the computation.
-print_opt                print currently active options
-v                        verbose mode
-h OR --help              print this help and exit

"""

def set_selenoprofiles_tree_drawer_var(varname, value):
  globals()[varname]=value

command_line_synonyms={}

def_opt= {'temp':'/users/rg/mmariotti/temp', 'common':0, 'no_id':0, 'add':0, 'sp_add':0,
't':'/users/rg/mmariotti/Genomes/.tree', 'prompt':0,
'm':0,
's':'', 'f':'', 'out':'',
'a':0, 'C':0, 'c':0, 'explicit_labels':'selenocysteine,cysteine,homologue',
'I':0, 
'F':0, 'e':0, 'g':0, 'T':0,'f':0,
'o':0,
'O':'', 
'tsize': 20, 'ssize': 20, 'bsize': 8, 'nsize':14, 'margin_boxes':10,
'W':150, 'H':40,
'v':0, 'Q':0, 'ali':0, 'abs_pos':0, 'intron_color':'#FFFFFF',
'img_h':None, 'img_w':None, 
'*':[]
}


#########################################################
###### start main program function
                                                      
global label_to_color; label_to_color={'selenocysteine':'#89C900', #green 
'cysteine': '#F40000', #red
'arginine': '#F15BA6', #pink
'threonine': '#8B1B8D', #dark purple
'pseudo':    '#545454', #  gray
'uga_containing': '#EEA347',
'unaligned':'#979796','gapped':'#979796', #grey
'homologue':'#E1D100', #yellow 
'unknown':'#979796', #light grey
'serine':'#D09E5F',
'readthrough':'#0E6DBB', #blue
'tRNA':'#C5D51B',
'glycine':'#FFB20B', #orange
'leucine':'#834D1A', #brown
'tRNA?':'#7A850D',}


global ordered_labels; ordered_labels=['selenocysteine', 'cysteine', 'arginine', 'threonine', 'uga_containing', 'unaligned', 'gapped',
'alanine', 'asparagine', 'aspartic_acid', 'glutamic_acid', 'glutamine', 'glycine', 'histidine', 'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'tryptophan', 'tyrosine', 'valine', 
'homologue','pseudo', 'readthrough', 'tRNA', 'tRNA?']

common_names= {'Ornithorhynchus_anatinus': 'Platypus', 'scientific_name': 'common name', 'Oryctolagus_cuniculus': 'Rabbit', 'Dipodomys_ordii': 'Kangaroo_rat', 'Sus_scrofa': 'Pig', 'Pongo_pygmaeus': 'Orangutan', 'Sorex_araneus': 'Shrew', 'Myotis_lucifugus': 'Microbat', 'Oryzias_latipes': 'Medaka', 'Monodelphis_domestica': 'Opossum', 'Anolis_carolinensis': 'Lizard', 'Callorhinchus_milii': 'Elephant_shark', 'Microcebus_murinus': 'Mouse_lemur', 'Lama_pacos': 'Llama', 'Taeniopygia_guttata': 'Finch', 'Pan_troglodytes': 'Chimpanzee', 'Gasterosteus_aculeatus': 'Stickleback', 'Homo_sapiens': 'Human', 'Tupaia_belangeri': 'Tree_shrew', 'Dasypus_novemcinctus': 'Armadillo', 'Macaca_mulatta': 'Macaque', 'Otolemur_garnettii': 'Galago', 'Spermophilus_tridecemlineatus': 'Squirrel', 'Rattus_norvegicus': 'Rat', 'Macropus_eugenii': 'Wallaby', 'Xenopus_{ch40}Silurana{ch41}_tropicalis': 'Frog', 'Gallus_gallus': 'Chicken', 'Erinaceus_europaeus': 'Hedgehog', 'Gorilla_gorilla': 'Gorilla', 'Mus_musculus': 'Mouse', 'Choloepus_hoffmanni': 'Sloth', 'Tursiops_truncatus': 'Dolphin', 'Takifugu_rubripes': 'Fugu', 'Felis_catus': 'Cat', 'Callithrix_jacchus': 'Marmoset', 'Bos_taurus': 'Cow', 'Equus_caballus': 'Horse', 'Canis_lupus_familiaris': 'Dog', 'Pteropus_vampyrus': 'Macrobat', 'Danio_rerio': 'Zebrafish', 'Procavia_capensis': 'Hyrax', 'Loxodonta_africana': 'Elephant', 'Cavia_porcellus': 'Guinea_pig', 'Tetraodon_nigroviridis': 'Pufferfish', 'Tarsius_syrichta': 'Tarsier', 'Branchiostoma_floridae':'Lancelet', 'Drosophila_melanogaster':'Fruit_fly', 'Anopheles_gambiae':'Malaria_mosquito', 'Aedes_aegypti':'Yellow_fever_mosquito', 'Acyrthosiphon_pisum':'Pea_aphid', 'Pediculus_humanus':'Louse', 'Tribolium_castaneum':'Beetle', 'Apis_mellifera': 'Honey_bee', 'Nasonia_vitripennis':'Jewel_wasp'}

global warnings_reduced_text; warnings_reduced_text={} #keys: categories

def reduce_font_if_necessary(simpleTextItem, width=-1, height=-1, category='unknown'):
  psize=99 #not used, this is just for the first time last condition of while is checked
  reducing_of_points=0
  while (width>0 and simpleTextItem.boundingRect().width() > width) or (height>0 and simpleTextItem.boundingRect().height()> height ) and psize>0:
    font=simpleTextItem.font()
    psize= font.pointSize()
    font.setPointSize(psize-1)
    simpleTextItem.setFont(font)
    reducing_of_points+=1
  if reducing_of_points:  
    if not warnings_reduced_text.has_key(category): warnings_reduced_text[category]=0
    warnings_reduced_text[category]+=1
    #printerr('Reduced point size for text "'+simpleTextItem.text()+'" in category: '+category, 1)
  
class NumberedBoxFace(faces.Face):
  """This face displays numbers in colored boxes. Initialise it with a list like [ [number, '#hexcodeforcolor'], [], [number, '#hexcodeforcolor']  ...  ]
  Each element of this list is a two-element list with the number to write and the corresponding color in hex code. Empty elements are display like white boxes (nothing is drawn in this position) this is useful to keep boxes in different lines still aligned 
 """
  def __init__(self, input_list):
      faces.Face.__init__(self)
      self.type = "item"
      self.item = None
      self.list=input_list
        
  def update_items(self):
    self.item = QtGui.QGraphicsRectItem(0, 0, numbered_box_width*len(self.list), numbered_box_height) #parent item: bkg
    #self.item.setBrush(QtGui.QBrush(QtGui.QColor(self.gene.color())))
    for boxdata_index in range(len(self.list)):
      boxdata=self.list[boxdata_index]
      if boxdata:
        number, hexcolor = boxdata        
        colored_box=     QtGui.QGraphicsRectItem(numbered_box_width*boxdata_index, 0, numbered_box_width, numbered_box_height)
        colored_box.setBrush(QtGui.QBrush(QtGui.QColor( hexcolor )))
        text_in_box=     QtGui.QGraphicsSimpleTextItem() 
        #printerr("***"+str(number), 1)
        text_in_box.setText(str(number))	
        font=text_in_box.font();  font.setPointSize(opt['nsize']);     text_in_box.setFont(font) #setting default text size         
        reduce_font_if_necessary(text_in_box, numbered_box_width-1, numbered_box_height, category='numbered box' )
        text_in_box.setZValue(1) #default is 0, this is to be sure it is on top of colored box  
        text_in_box.setPos(boxdata_index*numbered_box_width+1, 0)        
        colored_box.setParentItem(self.item)
        text_in_box.setParentItem(self.item)

  def _width(self):        return self.item.rect().width()
  def _height(self):       return self.item.rect().height()


class GeneFace(faces.Face):
    """ Colored rectangle for a gene, the color represent the label, the width and position how it spans the profile, the black or white lines the introns"""
    def __init__(self, limited_p2ghit_instance):
      faces.Face.__init__(self)
      self.type = "item"
      self.item = None
      self.gene= limited_p2ghit_instance
        
    def update_items(self):
      offset_for_additional_here=offset_for_additional*len(self.gene.additional_features)      
      self.item = QtGui.QGraphicsRectItem(0, 0, gene_brick_width+offset_for_id+offset_for_additional_here, gene_brick_height) #parent item: bkg
      self.item.setPen(QtGui.QPen(QtCore.Qt.NoPen)) #no black border      
      if not opt['no_id']:
        ### rectangle for prediction id
        chrom_id_rect = QtGui.QGraphicsRectItem(0, gene_brick_height/6, offset_for_id, gene_brick_height*2/3) #parent item: bkg
        chrom_id_rect.setParentItem(self.item)
        ## text for prediction id
        obj=QtGui.QGraphicsSimpleTextItem() 
        obj.setText(self.gene.id)
        reduce_font_if_necessary(obj, offset_for_id, gene_brick_height*2/3, category='prediction id')
        obj.setPos( (offset_for_id-obj.boundingRect().width())/2 , (gene_brick_height-obj.boundingRect().height())/2 ) #centering the text
        obj.setParentItem(chrom_id_rect)

      for index, x in enumerate( self.gene.additional_features ):
        #write(self.gene.full_id+' drawing secis '+' '+str(offset_for_additional), 1)
        secis_rect = QtGui.QGraphicsRectItem(gene_brick_width+offset_for_id+offset_for_additional*index, gene_brick_height/8, offset_for_additional-1, gene_brick_height*6/8) #parent item: bkg
        secis_rect.setParentItem(self.item)
        secis_rect.setBrush(QtGui.QBrush(QtGui.QColor( x.color()    )))
        #pen=QtGui.QPen(); pen.width=3
        #pen.setColor(QtGui.QColor("#FFFFFF"))
        #secis_rect.setPen(QtGui.QPen(color="#FFFFFF", width=2) )
        obj=QtGui.QGraphicsSimpleTextItem()         
        obj.setText(x.text)
        obj.setBrush( QtGui.QBrush(QtGui.QColor( "#FFFFFF"  ) ))
        reduce_font_if_necessary(obj, offset_for_additional_here, gene_brick_height*6/8, category='secis description')
        obj.setParentItem(secis_rect)
        obj.setPos(  gene_brick_width+offset_for_id+offset_for_additional*index+ (offset_for_additional-obj.boundingRect().width())/2 , (gene_brick_height-obj.boundingRect().height())/2 ) #centering the text

      #drawing line to represent the 100 coverage for profile
      line_for_full_coverage= QtGui.QGraphicsLineItem(offset_for_id, 0, offset_for_id+gene_brick_width, 0)
      line_for_full_coverage.setParentItem(self.item)
      ## rectangle for gene brick
      gene_brick=QtGui.QGraphicsRectItem(offset_for_id+self.gene.relative_boundaries()[0]*gene_brick_width, 0, (self.gene.relative_boundaries()[1]-self.gene.relative_boundaries()[0])*gene_brick_width, gene_brick_height) 
      gene_brick.setBrush(QtGui.QBrush(QtGui.QColor(self.gene.color())))
      gene_brick.setParentItem(self.item)

      #drawing lines for introns
      if not opt['I'] and len(self.gene.exons)>1:
        cds_so_far=0
        tot_cds=self.gene.length()
        for exon_index in range(len(self.gene.exons[:-1])):
          start, end = self.gene.exons[exon_index]
          cds_so_far+= end-start+1
          aa_so_far = cds_so_far/3.0
          if self.gene.strand=='+':     intron_length= self.gene.exons[exon_index+1][0] - end -1   #######
          elif self.gene.strand=='-':   intron_length= start - self.gene.exons[exon_index+1][1] -1   #######
          x_intron  =  self.gene.relative_position_in_ali_of(aa_so_far) * gene_brick_width  +  offset_for_id
          line_intron= QtGui.QGraphicsLineItem( x_intron, 1, x_intron, gene_brick_height-2)
          
          color=opt['intron_color'] # '#FFFFFF' #white
          if intron_length<=5: color='#EE0000' #red for frameshifts
          line_intron.setPen(QtGui.QPen(QtGui.QColor(color)))
          line_intron.setZValue(1)
          line_intron.setParentItem(gene_brick)          

      if opt['f']:
        for feature_name in self.gene.graphical_features: 
          for aa_position in self.gene.graphical_features[feature_name]:
            x_feature= self.gene.relative_position_in_ali_of(aa_position) * gene_brick_width  +  offset_for_id
            line_feature = QtGui.QGraphicsLineItem( x_feature, 1, x_feature, gene_brick_height-2)
            if graphical_features_colors.has_key(feature_name):
              line_feature.setPen(QtGui.QPen(QtGui.QColor( graphical_features_colors[feature_name]  )))
            line_feature.setZValue(1.1)
            line_feature.setParentItem(gene_brick)
            
      if not opt['T']:
      ## text for chromosome
        obj=QtGui.QGraphicsSimpleTextItem() 
        obj.setPos(offset_for_id+1, 0)
        obj.setText(self.gene.chromosome)
        font=obj.font();  font.setPointSize(opt['bsize']);     obj.setFont(font) #setting default text size
        reduce_font_if_necessary(obj, gene_brick_width, category='text for chromosome')
        obj.setZValue(2)
        obj.setParentItem(gene_brick)
      ## text for positions
        obj=QtGui.QGraphicsSimpleTextItem() 
        obj.setPos(offset_for_id+1,gene_brick_height/2 +1 )
        obj.setText(join([str(i) for i in self.gene.boundaries()], self.gene.strand)) # e.g. 1+101 (positive strand) 45-80 (negative strand
        font=obj.font();  font.setPointSize(opt['bsize']);     obj.setFont(font) #setting default text size      
        reduce_font_if_necessary(obj, gene_brick_width, category='positions text')
        obj.setZValue(2.1)
        obj.setParentItem(gene_brick)

    def _width(self):        return self.item.rect().width()
    def _height(self):       return self.item.rect().height()

def is_selenoprofiles_title( title ):
  """ Returns True if it is a selenoprofiles2 title, False if not """
  if 'chromosome:' in title and 'strand:' in title and  'positions:' in title and title.split()[0].count('.') in [4, 2]: return True
  return False

class gene_attribute(object):
  """ """
  def __init__(self, text=None, type=None, color=None):
    self.type=type
    self.text=text
    self.data={}
    self.custom_color=color
  def color(self):
    if self.custom_color: return self.custom_color
    else: return "#000000"  

class limited_p2ghit(gene):
  """This class is analog to p2ghit, but lacks some of its data. """
  def load_from_header(self, header):
    gene.load_from_header(self, header)
    if 'species:' in header:        
      tline=header.split('species:')[1]
      if tline.startswith('"'): species_name= tline[1:].split('"')[0]
      else:                     species_name= tline.split()[0]
    elif len(self.id.split('.'))>=5:     species_name= unmask_characters(replace_chars(self.id.split('.')[3], '_', ' '))       
    else: species_name='None'
    self.species= species(  species_name  )
    #self.program= header.split('prediction_program:')[1].split()[0]
    self.label  =         self.id.split('.')[2]
    self.profile_name=    self.id.split('.')[0]
    if self.id.split('.')>=5:  self.target_name = self.id.split('.')[-1]
    else:                      self.target_name = base_filename(self.target)
    self.full_id=self.id
    self.id=self.id.split('.')[1]
    self.relative_boundaries_data=[]  #filled when the relative_boundaries function is called
    self.graphical_features={}
    for feature_name in graphical_features_names:
      if feature_name+':' in header:
        tt=header.split(feature_name+':')[1]
        if tt and not tt[0] in ' \n': 
          self.graphical_features[feature_name]=[float(n) for n in tt.split()[0].split(',')]
    self.additional_features=[]
    
  def summary(self):    return gene.summary(self, other_fields=['program', 'label', 'profile_name'])
    
  def relative_boundaries(self):
    """Return a list of two float numbers, indicating the boundaries in percentage of where the prediction is spanning the profile. max: 0.000001, 1.0 """
    if not self.relative_boundaries_data:
      try:      self.profile; assert self.profile
      except:   raise Exception, "relative_boundaries ERROR the .profile attribute must be defined to call this method! it failed on this object: "+str(self)
      self.relative_boundaries_data=  [self.relative_position_in_ali_of(1), self.relative_position_in_ali_of(len(nogap(self.seq)))]
    return self.relative_boundaries_data

  def relative_position_in_ali_of(self, position):
    """ Returns the position in the alignment corresponding to position pos in the aminacid sequence of self. Analog to position_in_ali of class alignment, but specific for this class. input is 1 based, output is a float, max=1.0 ... like relative boundaries"""
    pos_seq=0
    for p, aa in enumerate(self.seq):
      if aa!='-':        pos_seq+=1
      if pos_seq>=position:        return float(p+1)/len(self.seq)

  def color(self):   
    """Return the hex representation of the color with which this gene will be drawn, depending on the label of this object. """
    if label_to_color.has_key(self.label): return label_to_color[self.label]
    else: return label_to_color['unknown']

  def target_name(self):
    """Return the filename of the target for this prediction, removing the extension """
    return join(base_filename(self.target).split('.')[:-1], '.')

  def sequence_identity_with_profile(self, dont_count_gaps=0):
    """This maps the prediction back to the profile and returns the average sequence identity with the sequences of the profile. dont_count_gaps is flag for the following behavior:                                                                                                                                                                             
    0 -> terminal gaps are not counted                                                                                                                                           
    1 -> gaps are not counted                                                                                                                                                    
    2 -> everything is counted (resulting score will be much lower than the other methods for uncomplete hits)                                                                   
  """
    dont_count_gaps=int(dont_count_gaps)
    dcg= bool(dont_count_gaps==1)
    dctg=bool(dont_count_gaps==0)
    return self.alignment_with_profile(dont_shrink=True, title='UnMaTcHaBlE').average_sequence_identity_of(title='UnMaTcHaBlE', dont_count_gaps=dcg, dont_count_terminal_gaps=dctg)

  def alignment_with_profile(self, profile_ali='', dont_shrink=False, title=''):
    """ This functions is designed to be equivalent to the one with the same name in the p2ghit class of selenoprofiles, but to be working for limited_p2ghit
    """
    if not title: title=self.output_id()
    a=self.profile.copy()
    a.add(title, self.seq)
    a.remove_empty_columns()
    return a

  def output_id(self):
    return self.profile_name+'.'+self.id+'.'+self.label+'.'+self.species.name+'.'+self.target_name

graphical_features_names=[]
graphical_features_colors={}
gene_brick_width=150
gene_brick_height=40
offset_for_id=25
offset_for_additional=25

def main():
#########################################################
############ loading options
  global opt; opt=command_line(def_opt, help_msg, '*', synonyms=command_line_synonyms )
  global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)

  global gene_brick_width;  gene_brick_width=opt['W']
  global gene_brick_height; gene_brick_height=opt['H']
  global numbered_box_width;  numbered_box_width=30 
  global numbered_box_height; numbered_box_height=40

  global offset_for_id
  if opt['no_id']: offset_for_id=0
  global offset_for_additional

  global explicit_labels; explicit_labels=opt['explicit_labels'].split(',')
  global graphical_features_names
  global graphical_features_colors
  if opt['f']:
    for piece in opt['f'].split(','):
      splt=piece.split('#')
      gf_name=splt[0]
      if len(splt)>1:        graphical_features_colors[gf_name]= '#'+splt[1]
      graphical_features_names.append( gf_name )

  #checking input
  global tree_input_file;
  if opt['t']: 
    tree_input_file=opt['t']
    check_file_presence(tree_input_file, 'tree_input_file')
    t=PhyloTree(tree_input_file)
    if opt['m']:
      for node in t.traverse():  node.name =  unmask_characters(replace(node.name, '_', ' '))
    for node in t:
      if node.is_leaf():
        node.is_used=0
        node.columns={} #indexed with family name
  else:
    print 'WARNING no tree was provided in input: an unstructured tree with all species encountered in the input alignment will be used. To derive the phylogenetic tree of species included in ncbi taxonomy, visit: https://github.com/jhcepas/ncbi_taxonomy '
    t=PhyloTree()

  tree_style = TreeStyle()  
  if opt['C']:   
    tree_style.mode='c'
    #tree_style.scale *= 10
    tree_style.allow_face_overlap = True
  else:          tree_style.mode='r'
  tree_style.branch_vertical_margin = 12
  tree_style.draw_aligned_faces_as_table = True
  tree_style.aligned_table_style = 1
  tree_style.show_leaf_name = False

  global p2ghit_by_name; p2ghit_by_name={}
  global   families_order; families_order=[]
  global labels_seen_for_family_hash;         labels_seen_for_family_hash={}
  for ali_file in opt['*']:
    if ':' in ali_file and '-' in ali_file.split(':')[1]:
      ali_file_position_range=   [int(i) for i in ali_file.split(':')[1].split('-')]  # [start, stop]
      ali_file=ali_file.split(':')[0]
    
    print "Input alignment: "+ali_file
    ali=alignment(ali_file) 

    profile_ali=alignment()
    for title in ali.titles():
      if not is_selenoprofiles_title(title):
        profile_ali.add( title, ali.seq_of(title) )

    if not profile_ali.nseq():      profile_ali.add( 'puppet', "X"*profile_ali.length() )
    #profile_ali.remove_useless_gaps()

    any_title_was_included=False
    for title in ali.titles():
      #print "Title: "+title
      seq=ali.seq_of(title)
      seq_no_terminal_gaps=seq 
      while seq_no_terminal_gaps.startswith('-'):   seq_no_terminal_gaps=seq_no_terminal_gaps[1:]
      while seq_no_terminal_gaps[-1]=='-':          seq_no_terminal_gaps=seq_no_terminal_gaps[:-1]

      
      if is_selenoprofiles_title( title ):
        any_title_was_included=True
        x=limited_p2ghit()
        x.load_from_header(title)
        x.seq=seq
        x.profile=profile_ali
        family=x.profile_name
        profile_ali.name = family ## will be the same for each selenoprofiles title in this alignment . nonetheless we don't want to know it before here otherwise we'd need to make assumptions

        species_name = x.species.name
        
        p2ghit_by_name[  x.full_id  ] =x
        if not labels_seen_for_family_hash.has_key(family):               labels_seen_for_family_hash[family]={}      
        if not labels_seen_for_family_hash[family].has_key(x.label):      
          labels_seen_for_family_hash[family][x.label]=1
          
          if opt['a'] and opt['c'] and not  x.label in explicit_labels:   labels_seen_for_family_hash[family]['others']=1

        if not opt['t']: #no tree was specified. Let's build a puppet tree with the organisms names
          try: 
            node=t&species_name
          except:   
            node=t.add_child(name=species_name)
            node.is_used=0
            node.columns={} #indexed with family name
            node=t&species_name    
        try:
          node=t&species_name
          node.is_used=1

          if not node.is_leaf(): raise Exception, "ERROR species: "+species_name+" is not a leaf in the input tree!"
          if not node.columns.has_key(family): node.columns[family]=[]
          node.columns[family].append( x )


        except ValueError: 
#          print species_name
#          raise
          if not opt['g']:        raise Exception, "ERROR can't find a node in the tree for species: "+species_name
          else:                   print "ignoring species not found: "+species_name

    if any_title_was_included:   
      families_order.append(profile_ali.name)
      title_face=faces.TextFace(family, fsize=opt['tsize'])
      title_face.hz_align = 1
      tree_style.aligned_header.add_face(title_face, column=len(families_order))
    else:
      printerr("WARNING no single title was included for alignment: " +ali_file, 1)

  if opt['add']:
    for line in open(opt['add']):
      stripped=line.strip()
      if stripped:
        splt=stripped.split('\t')
        p2gname=splt[0]
        text=replace (   splt[1]  , '\\n', '\n')
        if len(splt)>2 and splt[2]:           color=splt[2]
        else: color=None

        p2ghit_by_name[p2gname].additional_features.append( gene_attribute(text=text, color=color)  )


  if opt['sp_add']:
    for line in open(opt['sp_add']):
      stripped=line.strip()
      if stripped:
        splt=stripped.split('\t')
        species=splt[0].strip()
        try: 
          node = t&species
        except:
          try:
            species=unmask_characters( replace(species, '_', ' ') )
            node = t&species
          except:
            raise Exception, "ERROR can't find a node in the tree for species in -sp_add file: "+str([species])
          
          value=splt[1]
          if value.startswith('COLOR'):
            node.species_coloring='#'+value.split('#')[1]
          else:
            if not hasattr(node, 'species_attributes'):              node.species_attributes=[]
            node.species_attributes.append(value)


        
  
      


  #parse tree and prune useless nodes!
  if not opt['e']: 
    nodes_to_keep=[t]
    for node in t:
      if node.is_used: 
        nodes_to_keep.append(node)

    #print [i.name for i in nodes_to_keep]
    #print t
    #prune_tree(t, nodes_to_keep)
    #print t
    t.prune(nodes_to_keep)
    
  #setting outgroup
  if opt['o']:
    if opt['o'] in [1, True]:
      t.set_outgroup(t.get_midpoint_outgroup())
    else:
      outgroup=str(opt['o'])
      matches = t.search_nodes(name = outgroup)
      #print "*********"+str(matches)
      if len(matches) > 0: t.set_outgroup(matches[0].name)
      else: raise Exception, "ERROR the species "+outgroup+' was not found in the tree.' +" Maybe it is because this species had no results and was pruned out: try with option -e"*int(not opt['e'])
        
  t.ladderize()
  #print t
  if not opt['out']: t.show(mylayout, tree_style) 
  else:              t.render(opt['out'], layout=mylayout, tree_style=tree_style, h=opt['img_h'], w=opt['img_w'])

  for category in warnings_reduced_text:
    printerr('WARNING '+str(warnings_reduced_text[category])+' '+category+' text(s) were reduced in point size to fit the dedicated space', 1)

  if opt['prompt']: interactive_mode(message='Tree is loaded into variable: t')()

  ###############

def mylayout(node):

  ## Phylogeny layout
  if node.is_leaf():
    ## Setting the leaf color name

    #node.img_style['bgcolor']='#DDDDDD'
    column_index=0

    ## modify to have a customized species name printed
    main_species=node.name
    main_species_in_filenames=replace_chars(mask_characters(main_species), ' ', '_')
    if opt['common']:
      main_species = common_names.setdefault(main_species_in_filenames, main_species_in_filenames)

    if hasattr(node, 'species_coloring'): fgcolor=node.species_coloring
    else:                                 fgcolor="#000000"
    
    nameFace=faces.TextFace(main_species, fgcolor=fgcolor, fsize=opt['ssize'])         #species name!
    nameFace.vt_align=1
    faces.add_face_to_node(nameFace, node, column = column_index, position="branch-right")
    column_index+=1

    if hasattr(node, 'species_attributes'):
      for value in node.species_attributes:
        valueFace=  faces.TextFace(value, fgcolor="#000000", fsize=opt['ssize'])
        valueFace.vt_align=1
        faces.add_face_to_node(valueFace, node, column = column_index, position="branch-right")    
      column_index+=1        


    for family in families_order:
      
      if opt['a']:
      ####################### ABSTRACT MODE
        if node.columns.has_key(family):
          if opt['F']:
            familyNameFace=faces.TextFace(family, fsize=opt['bsize'], fgcolor="#000000")
            familyNameFace.margin_left = 20;             familyNameFace.margin_right = 20; 
            faces.add_face_to_node(familyNameFace, node, column = column_index, aligned = True)

          count_per_label={}
          for gene_index in  range(len(node.columns[family])):
            x=node.columns[family][gene_index]
            if not count_per_label.has_key(x.label): count_per_label[x.label]=0
            count_per_label[x.label]+=1

          #create rectangle with width of: box_width * len(labels_seen_for_family_hash[family].keys())

          list_to_build_numbered_box=[]
          if not opt['c']:
            for label in ordered_labels:  
              if labels_seen_for_family_hash[family].has_key(label):
                if not count_per_label.has_key(label):   list_to_build_numbered_box.append( [] )
                else:                                    list_to_build_numbered_box.append( [count_per_label[label], label_to_color.setdefault(label, label_to_color['unknown']) ] )
          else: #condensating boxes in max n columns. possible labels are defined by explicit_labels option

            count_per_label['others']=0
            for label in ordered_labels:  
              if count_per_label.has_key(label) and not label in explicit_labels:                 count_per_label['others']+=count_per_label[label]
            for label in explicit_labels+['others']:
              if labels_seen_for_family_hash[family].has_key(label):
                if not count_per_label.has_key(label):   list_to_build_numbered_box.append( [] )
                elif not count_per_label[label]:         list_to_build_numbered_box.append( [] ) 
                else:                                    list_to_build_numbered_box.append( [count_per_label[label], label_to_color.setdefault(label, label_to_color['unknown']) ] )
          numbered_boxes_face=NumberedBoxFace(list_to_build_numbered_box)
          numbered_boxes_face.margin_right = opt['margin_boxes']
          numbered_boxes_face.hz_align=1
          numbered_boxes_face.rotable = False
          faces.add_face_to_node(numbered_boxes_face, node, column = column_index, aligned = True)                        
          
        else:
          if opt['F']:
            emptyFamilyNameFace=faces.TextFace(family, fgcolor="#FFFFFF")
            emptyFamilyNameFace.margin_left=4;       emptyFamilyNameFace.margin_right=2; 
            faces.add_face_to_node(emptyFamilyNameFace, node, column = column_index, aligned = True)

      else:
      ####################### NORMAL MODE      
        if node.columns.has_key(family):
          if opt['F']:
            familyNameFace=faces.TextFace(family, fgcolor="#000000")
            familyNameFace.margin_left=2;       familyNameFace.margin_right=2;             
            faces.add_face_to_node(familyNameFace, node, column = column_index, aligned = True)

          list_of_genes=node.columns[family]
          try:            list_of_genes.sort(key=lambda x:ordered_labels.index(x.label)) # this list will host the ordered list of genes to draw in this species. The ordered is determined by the appearances of labels in ordered_labels
          except: 
            for g in list_of_genes:  
              if not g.label in ordered_labels:
                print " WARNING label "+g.label + " was not found among the known ones. "
                
         
          for gene_index in  range(len(list_of_genes)):
            x=list_of_genes[gene_index]
            gene_face=GeneFace(x)
            gene_face.margin_left=5;       gene_face.margin_right=5; 
            faces.add_face_to_node(gene_face, node, column = column_index, aligned = True)
            
            
#            if opt['ali']:
              
            
          #add separator?
          
        else:
          if opt['F']:
            emptyFamilyNameFace=faces.TextFace(family, fgcolor="#FFFFFF")
            emptyFamilyNameFace.margin_left=2;       emptyFamilyNameFace.margin_right=2; 
            faces.add_face_to_node(emptyFamilyNameFace, node, column = column_index, aligned = True)

      column_index+=1





#######################################################################################################################################
def prune_tree(t, nodes_to_keep):
    """ Fixed (and faster) prunning algorithm. Use this until I fix
    the problem within the main ETE branch.

    'nodes_to_keep' must be the list of node instances that you want
    to keep in the final tree. All nodes must be leaves, if not, they
    are automatically converted into leaves by removing their
    children.

    So far, this function is quite verbose. Printing slows down a bit
    the process, but you can follow the progress...
    """ 
    #print "Getting tree path..."
    # Converts to set to speed up searches
    if type(nodes_to_keep) == set:
        to_keep=nodes_to_keep
    else:
        to_keep=set(nodes_to_keep)

    #print "Checking that all nodes are leaves..."
    not_leaves = [n for n in nodes_to_keep if not n.is_leaf()]
    if len(not_leaves)>0:
    #    print "\nFixing", len(not_leaves), "non-leaf nodes..."
        # Converts all internal species nodes into leaves by removing all
        # their sub-species or strains
        for nl in not_leaves: 
            for c in nl.get_children():
                c.detach()
    to_detach = []
    
    #print "Selecting unused nodes"
    counter = 0
    for node in t.traverse("postorder"):
    #    print "\r", counter,
        counter +=1
        for c in node.children:
            if c in to_keep:
                to_keep.add(node)
                break
        if node not in to_keep:
            to_detach.append(node)
            for c in node.children:
                to_detach.remove(c)
    #print "\nDetaching", len(to_detach), "nodes"
    counter = 0
    for node in to_detach:
    #    print "\r", counter,
        counter +=1
        node.detach()
    #print "\nFixing", len(to_keep), "orphan nodes"
    counter = 0
    for node in to_keep:
    #    print "\r", counter,
        counter +=1
        if len(node.children) == 1:
            node.delete()
            

    if len(t.children)==1:
      try:
        a=t.children[0]
        a.delete()

      except: pass

    
    ############################            
            
    return t
 ##########################3
 



def close_program():
  if opt['debug']: raw_input('check temp folder:'+temp_folder)
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
    
else:
  global opt;
  opt=get_MMlib_var('opt')

