#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/08/17 16:46:02

# easy_install fisher/ if not also in extra stats package
try:
    from extra_stats.fisher import pvalue
except ImportError:
    from fisher import pvalue

# in my extra_stats package
import sys
sys.path.append('/home/francisco/toolbox/utils/')

from extra_stats.fdr import bh_qvalues
from numpy import log

from optparse import OptionParser
from bisect import bisect_left

__version__ = "0.10"
__title__   = "gene set tool kit v%s" % __version__


class Gene_set:
    '''
    Fatiscan with upper case, it is an object.
    from gene_set import Gene_set
    infile = '/home/francisco/project/functional_anaysis_vs_evolution/v_56/Mammals/0_dataset/Homo_sapiens.dN_val'
    annot  = '/home/francisco/project/functional_anaysis_vs_evolution/v_56/Mammals/funcDB/biol_proc_2-8.annot'
    gaga = Gene_set (infile, annot)
    lala = gaga.run_gsea()
    '''
    def __init__(self, infile, annot, partitions=30, use_order=True):
        '''
        init function, what is done when object is called.
        '''
        # get gene list and corresponding values
        self.infile = infile
        self.use_order = use_order
        self.genes, self.values, self.order = self._parse_infile()

        self.annot = self._parse_annot(annot)
        # sort genes in annot by their values...
        # useful to know which gene in list1 or list2
        self.annot = self._order_genes_in_annot()

        self.gsea_dic = {}

    def _parse_infile(self):
        '''
        parse in file in format:
        geneID_1 <tab> value1
        geneID_2 <tab> value2
        ...
        genes should be ordered by value
        returns genes, values and order of values
        '''
        genes, values = zip (*sorted ((i.strip().split('\t') \
                              for i in open(self.infile)), \
                             key=lambda x: float(x[1])))
        values = map (float, values)
        if self.use_order:
            order = map (values.index, values)
        else:
            order = values[:]
        return genes, dict (zip (genes, values)), order

    def _parse_annot(self, annot):
        '''
        parse annotations file in format:
        annotationsA <tab> geneID_1
        annotationsA <tab> geneID_2
        ...
        speed notes: * iterator on for for
                     * dico
        '''
        dico = {}
        for gene, annot in (i.strip().split('\t') for i in open (annot)):
            # only store genes that we have in our list
            if self.values.has_key(gene):
                dico.setdefault (annot, []).append (gene)
        return dico

    def _order_genes_in_annot(self):
        '''
        order genes in annot dict by their values
        '''
        dico = {}
        for annot in self.annot.iterkeys():
            dico[annot] = sorted (self.annot[annot], \
                                  key=lambda x: self.values[x])
        return dico

    def run_fatigo (self, list1, list2):
        '''
        computes an enrichment test between two lists of genes
        '''
        

    def run_gsea (self, partitions=30):
        '''
        run gsea needs python fisher, and fdr from extra stats
        speed notes: * making annot genes and order local does not
                       significantly speed up process.
                     * putting external method with bissect part neither
                     * no more ideas...
        '''
        pvalues     = []
        pos         = []
        total_len   = len (self.genes)
        # create local functions to skip call... faster
        append_pos  = pos.append
        append_pv   = pvalues.append
        order_index = self.order.index
        iter_annot  = self.annot.iteritems
        # intialize dict
        dico = dict (((k, []) for k in self.annot))
        # define part size. order[-1] == max(order)
        rank = float (len (self.order)-1)/partitions
        # define cutoff value for each partition
        dico['thresh'] = (bisect_left (self.order, rank * (part + 1)) \
                          for part in xrange(partitions))
        # start fishers
        for part in xrange(partitions):
            try:
                genes1     = set (self.genes[:order_index (dico['thresh'].next())])
            except ValueError:
                continue
            len_genes1 = len (genes1)
            len_genes2 = total_len - len_genes1
            local_def  = genes1.intersection
            for annot, annot_genes in iter_annot():
                append_pos ((part, annot))
                #p1 = len (annot_genes  & genes1)
                p1 = len (local_def ( set (annot_genes)))
                p2 = len (annot_genes) - p1
                n1 = len_genes1  - p1
                n2 = len_genes2  - p2
                dico[annot].append ({'p1': p1, 'n1': n1, 'p2': p2, 'n2': n2})
                append_pv (pvalue (p1, n1, p2, n2).two_tail)
        # compute adjustment of pvalues
        qvalues = iter (bh_qvalues (pvalues))
        pvalues = iter (pvalues)
        for part, annot in pos:
            dico[annot][part]['apv'] = qvalues.next()
            dico[annot][part]['pv' ] = pvalues.next()
        # store this in Gene_set
        self.gsea_dic = dico

    def write_gsea (self, outfile, max_apv=1, all_parts=False):
        '''
        write to file, or pickle
        '''
        def _get_string(dico, annot):
            '''
            get string from gsea_dic current value
            '''
            string = []
            string.append (dico['p1' ])
            string.append (dico['n1' ])
            string.append (dico['p2' ])
            string.append (dico['n2' ])
            odd = 1
            try:
                odd = log ((float (dico['p1'])\
                            /dico['n1'])\
                           /(float (dico['p2'])\
                             /dico['n2']))
            except ZeroDivisionError:
                if dico['p2'] == 0:
                    odd = float ('inf')
                elif dico['n1'] == 0:
                    print >> stderr, "WARNING: Empty partition: " \
                          + str(part)
                    odd = float ('-inf')
                elif dico['n2'] == 0:
                    print >> stderr, "WARNING: Empty partition: " \
                          + str(part)
                    odd = float ('inf')
            string.append (odd)
            string.append (dico ['pv' ])
            string.append (dico ['apv'])
            string.append (', '.join (list (annot[:dico['p1']])))
            string.append (', '.join (list (annot[dico['p1']:])))
            return map (str, string)

        if self.gsea_dic == {}:
            print >> stderr, 'ERROR: you do not have run GSEA yet...'
            raise Exception
        cols = ['#term', 'part', 'term_size', \
                'list1_positives', 'list1_negatives', \
                'list2_positives', 'list2_negatives', \
                'odds_ratio_log', 'pvalue', 'adj_pvalue', \
                'list1_positive_ids', 'list2_positive_ids']
        out = open (outfile, 'w')
        out.write('\t'.join(cols)+'\n')
        if all_parts:
            for annot in filter (lambda x: x!= 'thresh', self.gsea_dic):
                for part in xrange(30):
                    if self.gsea_dic[annot][part]['apv'] > max_apv:
                        continue
                    string = _get_string(self.gsea_dic[annot][part], \
                                         self.annot[annot])
                    out.write ('\t'.join ([annot] + [str(part)] +string) + '\n')
        else:
            for annot in filter (lambda x: x!= 'thresh', self.gsea_dic):
                part = min (self.gsea_dic[annot], key=lambda x: x['apv'])
                if max_apv < part['apv']:
                    continue
                string = _get_string(part, self.annot[annot])
                part = self.gsea_dic[annot].index(part)
                out.write ('\t'.join ([annot] + [str(part)] + string) + '\n')
        out.close()

def main ():
    '''
    for direct command line call
    '''
    opts = get_options()
    gene_set = Gene_set(opts.infile, opts.annot)
    gene_set.run_gsea(partitions = opts.partitions)
    if opts.pickle:
        from cPickle import dump
        dump (open (outfile, 'w'), self)
    else:
        gene_set.write_gsea(opts.outfile, all_parts=opts.all_parts, \
                            max_apv=float(opts.max_apv))

def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
Gene set enrichment analysis                                                                                                
.                                                                                .
********************************************                                      
"""
        )

    parser.add_option('-i', dest='infile', metavar="PATH", \
                      help='''path to input file with a ranked list, in format:                             
                      geneID_1 <tab> value1                                                        
                      geneID_2 <tab> value2                                                         
                      ...
                      ''')
    parser.add_option('-a', dest='annot', metavar="PATH", \
                      help='''path to annotations file in format:                                           
                      annotationsA <tab> geneID_1                                                       
                      annotationsA <tab> geneID_2                                                   
                      ...
                      ''')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='path to output file tab separated file')
    parser.add_option('-R', '--use_rank', action='store_true', \
                      dest='use_order', default=False, \
                      help=\
                      '''[%default] Use rank of genes in stead of provided value to determine thresh value for delimiting different partitions.''')
    parser.add_option('-p', metavar="INT", dest='partitions', default=30, \
                      help='''[%default] Number of partitions.''')
    parser.add_option('--max_apv', metavar="FLOAT", dest='max_apv', default=1, \
                      help='''[%default] Only write to outfile results with adjusted pvalue higher than specified value.''')
    parser.add_option('--long', dest='all_parts', action='store_true', default=False, \
                      help='''[%default] Write results for all partitions.''')
    parser.add_option('--pickle', action='store_true', \
                      dest='pickle', default=False, \
                      help='[%default] Store results in python dict (cPickle) format.')
    parser.add_option('--verbose', action='store_true', \
                      dest='verb', default=False, \
                      help=\
                      '[%default] Talk a bit... ')
    opts = parser.parse_args()[0]
    if not opts.infile or not opts.annot or not opts.outfile:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
