#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/08/05 14:46:46

import os
from re import match, sub
from optparse import OptionParser



__version__ = "0.01b"
__title__   = "fatiscan_concatenator v%s" % __version__



def main():
    '''
    main
    '''
    opts = get_options()
    dico = fatiParser (opts)
    dicoAnnot = getAnnot (opts.annot)
    if opts.write == 'R' or opts.write == 'both':
        Writer(dico, opts.out, dicoAnnot, opts.lvl)
    if opts.write == 'py' or opts.write == 'both':
        from cPickle import dump
        out = open (_check_out (opts.out + '.dic'), 'w')
        dump (dico, out)
        out.close()
        
def _print_status(indir, col):
    '''
    if verbose print status of parsing
    '''
    os.system('clear')
    print '\n\nParsing files from fatiscan...'
    print '\t in: '+indir
    print 'parsing file: ' + col

def _check_out (out):
    '''
    check that we are not overwriting
    '''
    if os.path.exists(out):
        todo = ''
        while todo != 'R' and todo != 'C' and todo != 'Q':
            todo = raw_input('File %s already exists, [R]eplace it add \
[C]hange output name to %s_bis ([Q]uit): ' % (out, out))
        if todo == 'C':
            out = out + '_bis'
        elif todo == 'Q':
            exit()
    return out

def getAnnot(annot):
    '''
    Parses .annot file and return a dictionary, just to get
    description and level of annotations...
    Should have those columns: Level, id, term, Type, #Seqs,
    score, sequences
    '''
    dicoAnnot = {}
    for line in open(annot):
        if line.startswith('Level'):
            continue
        dicoAnnot [line.strip().split()[1]] = \
                  (int(line.strip().split('\t')[0]),line.strip().split('\t')[2])
    return dicoAnnot


def Writer(dico, out, dicoAnnot, lvl):
    '''
    write file for R script...
    soon deprecated... hope
    '''
    dicPercent = {}
    for s in dico[dico.keys()[0]].keys():
        if s == 'save' : continue
        dicPercent[s] = {
            'reds' : [],
            'blues': []
            }
    # get list of available levels:
    levels = list (set(zip(*dicoAnnot.values())[0]))
    lvl =  levels if lvl == 'all' else map (int, lvl.split(','))
    if len (set (lvl) & set (levels)) != len (lvl):
        exit('ERROR: level %s not available\n' %\
             (str (lvl)))
    ng = open(_check_out (out + '_' + str (lvl[0]) + '.ng') ,'w')
    np = open(_check_out (out + '_' + str (lvl[0]) + '.np') ,'w')
    ng.write('name\t'          +\
             sub('\tsave','',\
                 '\t'.join (sorted (dico[dico.keys()[0]].keys()))) \
             + '\n')
    np.write ('name\tdescr\tnumGenes\t' + \
              sub('\tsave', '', \
                  '\t'.join (sorted (dico[dico.keys()[0]].keys()))) \
              + '\n')
    ngs = ''
    for l in map (lambda x: x[1], \
                  sorted (map (lambda x: \
                               (dicoAnnot[x][1].lower(), x), dico.keys()))):
        if not dicoAnnot[l][0] in lvl:
            continue
        ngline = []
        npline = []
        if not dico[l]['save']: continue
        for c in sorted(dico[l].keys()):
            if c == 'save': continue
            k = 'bet'
            if dico[l][c]['bet']['up'] and dico[l][c].has_key('min'):
                k = 'min'
            elif dico[l][c].has_key('max'):
                k = 'max'
            if dico[l][c]['bet']['apv'] <= 0.05:
                ngline.append(str (\
                    int (\
                    abs (\
                    dico[l][c][k]['up'] - \
                            (dico[l][c][k]['p2']/\
                             (dico[l][c][k]['p1'] + dico[l][c][k]['p2']\
                              )))*100)\
                    ))
                if dico[l][c][k]['up']:
                    dicPercent[c]['reds'].append (float (ngline[-1]))
                else: dicPercent[c]['blues'].append(float(ngline[-1]))
                npline.append (str (dico[l][c]['bet']['up']*2 + \
                                    (dico[l][c]['bet']['apv']<0.001)+1))
                ngs = str (dico[l][c][k]['p1']+dico[l][c][k]['p2'])
            else:
                ngline.append('')
                npline.append('0')
        if not ngline.count('') > (len (ngline)-1):
            np.write (l + '\t' + dicoAnnot[l][1] + '\t' + ngs + '\t' \
                      + '\t'.join (npline) + '\n')
            ng.write(l+'\t'+'\t'.join(ngline)+'\n')
    ng.close()
    np.close()

def fatiParser(opts):
    '''
    parse fatiscan outfile, returns a big dictionary
    '''
    cols = []
    dico = {}
    for infile in filter (lambda x: not x.startswith('.'), \
                          sorted (os.listdir(opts.indir))):
        if not opts.pattern in infile:
            continue
        apvs = {}
        first = {}
        last  = {}
        col = infile
        cols.append(col)
        if opts.verbose:
            _print_status(opts.indir, col)
        fatiout = (open (opts.indir + infile)).readlines()
        for i in range(0, len(fatiout)):
            ks = []
            if (match('#', fatiout[int(i)]))\
                   or (match('^\n',fatiout[i]))\
                   or (match('^ $',fatiout[i])):
                continue
            line  = fatiout[i].strip().split('\t')
            if int(line[3]) + int(line[4]) < 6: continue
            apv = float(line[8])
            np  = float(line[0])
            try:
                up = (float(line[3])/float(line[4])) > \
                     (float(line[5])/float(line[6]))
            except ZeroDivisionError:
                continue
            name = line[1]
            if not first.has_key(name):
                first[name] = 30
                last[name] = 0
                apvs[name] = [apv]
            else:
                apvs[name].append(apv)
            if apv <= min (apvs[name]): ks.append('bet')
            if apv <= 0.05:
                if up and np < first[name] :
                    first[name] = np
                    ks.append('min')
                if up == False and np > last[name]:
                    last[name] = np
                    ks.append('max')
            if not dico.has_key(name):
                dico[name] = {'save' : False}
            if not dico[name].has_key(col):
                dico[name][col] = {}
            if apv <= 0.05 and not '.dw' in infile:
                dico[name]['save'] = True
            for k in ks:
                dico[name][col][k] = {}
                dico[name][col][k]['np'  ] = np + 1              
                dico[name][col][k]['lim' ] = float(line[2])
                dico[name][col][k]['p1'  ] = float(line[3])
                dico[name][col][k]['n1'  ] = float(line[4])
                dico[name][col][k]['p2'  ] = float(line[5])
                dico[name][col][k]['n2'  ] = float(line[6])
                dico[name][col][k]['apv' ] = apv
                dico[name][col][k]['up'  ] = up
                list1 = line[9].split(',')
                list2 = line[10].split(',')
                if list1[-1] == ''    : list1.pop()
                if list2[-1] == ''    : list2.pop()
                if list1[-1] == 'null': list1.pop()
                if list2[-1] == 'null': list2.pop()
                dico[name][col][k]['l1'  ] = list1
                dico[name][col][k]['l2'  ] = list2
                dico[name][col][k]['tot' ] = float(line[3])+float(line[5])
    return dico


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
        Concatenate outfile from fatiscan either short/long, into big
python dictionary or 2 tables for R script, or both
./concat_fatiscan.py -i /home/francisco/project/functional_anaysis_vs_evolution/v_56/Mammals/2_fatiOut/GO_2-8/ -o ici -a /home/francisco/project/functional_anaysis_vs_evolution/v_56/Mammals/funcDB/biol_proc_detail.txt --format both
.                                                                           .
********************************************                                      
TODO:                                                                                     
      many things....                                                                      
********************************************                                            
"""
        )
    parser.add_option('-i', dest='indir', metavar="PATH", \
                      help='path to directory containing outfiles of fatiscan.')
    parser.add_option('-o', dest='out', metavar="PATH", \
                      help='path to output file.')
    parser.add_option('--format', dest='write', default='R', \
                      help=\
                      '''[%default] outfile format, "R" to generate 2 matrix
                      readable by little R script, "py" to generate python
                      dictionary, "both".
                      ''')
    parser.add_option('-a', dest='annot', metavar="PATH", \
                      help=\
                      '''path to annotations file. Should have those columns:
                      Level, id, term, Type, #Seqs, score, sequences.
                      ''')
    parser.add_option('-v', dest='verbose', action='store_true', \
                      default=False, help=\
                      '''[%default] Verbose output.
                      ''')
    parser.add_option('-l', dest='lvl', default='all', metavar='LIST', help=\
                      '''[%default] level(s) of annotations you want to include
                      in outfiles. (e.g.: "-l 2,3")
                      '''
                      )
    parser.add_option('-p', dest='pattern', metavar="PATH", \
                      default='', help=
                      '''common pattern of fatiscan outfiles.       
                      e.g.: fati.out, will match with lala_fati.out but also
                      with lala_fati.out2etc.''')
    opts = parser.parse_args()[0]
    if not opts.indir or not opts.annot:
        exit(parser.print_help())
    return opts

if __name__ == "__main__":
    exit(main())
